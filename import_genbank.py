#!/usr/bin/env python
"""Import a GenBank results file into the antiSMASH database"""
from __future__ import print_function
import sys
import hashlib
from Bio import Entrez, SeqIO
import psycopg2
import psycopg2.extensions
psycopg2.extensions.register_type(psycopg2.extensions.UNICODE)
psycopg2.extensions.register_type(psycopg2.extensions.UNICODEARRAY)

DB_CONNECTION = "host='localhost' port=15432 user='postgres' password='secret' dbname='antismash'"
Entrez.email = "kblin@biosustain.dtu.dk"


def main():
    '''Run the import'''
    if len(sys.argv) < 2:
        print('Usage: {} <gbk file>'.format(sys.argv[0]), file=sys.stderr)
        sys.exit(2)

    connection = psycopg2.connect(DB_CONNECTION)

    rec = SeqIO.read(sys.argv[1], 'genbank')
    with connection:
        with connection.cursor() as cursor:
            load_record(rec, cursor)

    connection.close()


def load_record(rec, cur):
    '''Load a record into the database using the cursor'''
    genome_id = get_or_create_genome(rec, cur)
    print("genome_id: {}".format(genome_id))
    seq_id = get_or_create_dna_sequence(rec, cur, genome_id)
    print("seq_id: {}".format(seq_id))

    FEATURE_HANDLERS = {
        'CDS': handle_cds,
        'CDS_motif': handle_cds_motif,
        'aSDomain': handle_asdomain,
        'cluster': handle_cluster
    }

    for feature in rec.features:
        if feature.type not in FEATURE_HANDLERS:
            continue
        handler = FEATURE_HANDLERS[feature.type]
        handler(rec, cur, seq_id, feature)



def get_or_create_dna_sequence(rec, cur, genome_id):
    '''Fetch existing dna_sequence entry or create a new one'''
    params = {}
    params['seq'] = str(rec.seq)
    params['md5sum'] = hashlib.md5(params['seq']).hexdigest()
    params['accession'] = rec.annotations['accessions'][0]
    params['version'] = rec.annotations['sequence_version']

    cur.execute("SELECT sequence_id FROM antismash.dna_sequences WHERE md5 = %s", (params['md5sum'],))
    ret = cur.fetchone()
    if ret is None:
        params['genome'] = genome_id
        cur.execute("INSERT INTO antismash.dna_sequences (dna, md5, acc, version, genome)"
                    "VALUES (%(seq)s, %(md5sum)s, %(accession)s, %(version)s, %(genome)s) RETURNING sequence_id;", params)
        ret = cur.fetchone()

    return ret[0]


def get_or_create_genome(rec, cur):
    '''Fetch existing genome entry or create a new one'''
    taxid = get_taxid(rec)
    cur.execute("SELECT genome_id FROM antismash.genomes WHERE taxon = %s", (taxid,))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.genomes (taxon) VALUES (%s) RETURNING genome_id;", (taxid,))
        ret = cur.fetchone()

    return ret[0]


def get_taxid(rec):
    '''Extract the taxid from a record'''
    for feature in rec.features:
        if feature.type != 'source':
            continue
        if 'db_xref' not in feature.qualifiers:
            return 0
        for entry in feature.qualifiers['db_xref']:
            if entry.startswith('taxon:'):
                return int(entry[6:])


def get_or_create_locus(cur, seq_id, feature):
    '''Get or create a new locus tag'''
    params = {}
    params['start_pos'] = feature.location.nofuzzy_start
    params['end_pos'] = feature.location.nofuzzy_end
    if feature.location.strand == 1:
        params['strand'] = '+'
    elif feature.location.strand == -1:
        params['strand'] = '-'
    else:
        params['strand'] = '0'
    params['sequence'] = seq_id
    cur.execute("SELECT locus_id FROM antismash.loci WHERE sequence = %(sequence)s AND start_pos = %(start_pos)s AND end_pos = %(end_pos)s AND strand = %(strand)s",
                params)
    ret = cur.fetchone()

    if ret is None:
        cur.execute("INSERT INTO antismash.loci (start_pos, end_pos, strand, sequence) VALUES "
                    "(%(start_pos)s, %(end_pos)s, %(strand)s, %(sequence)s) RETURNING locus_id", params)
        ret = cur.fetchone()

    return ret[0]


def handle_cds(rec, cur, seq_id, feature):
    '''Handle CDS features'''
    pass


def handle_cds_motif(rec, cur, seq_id, feature):
    '''Handle CDS_motif features'''
    pass


def handle_asdomain(rec, cur, seq_id, feature):
    '''Handle aSDomain features'''
    pass


def handle_cluster(rec, cur, seq_id, feature):
    '''Handle cluster features'''
    params = {}
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)
    params['evidence'] = 'prediction'
    print("locus_id: {}".format(params['locus_id']))
    cur.execute("SELECT bgc_id FROM antismash.biosynthetic_gene_clusters WHERE locus = %(locus_id)s", params)
    ret = cur.fetchone()
    if ret is None:
        cur.execute("""
INSERT INTO antismash.biosynthetic_gene_clusters (locus, evidence)
SELECT val.locus, f.evidence_id FROM (
    VALUES (%(locus_id)s, %(evidence)s) ) val (locus, evidence)
LEFT JOIN antismash.evidences f ON val.evidence = f.name
RETURNING bgc_id
""", params)
        ret = cur.fetchone()
    params['bgc_id'] = ret[0]

    print("bgc_id: {}".format(params['bgc_id']))

    for product in feature.qualifiers['product'][0].split('-'):
        nx_create_rel_clusters_types(cur, params, product)

    print(feature.qualifiers)


def nx_create_rel_clusters_types(cur, params, product):
    '''create relation table to bgc_types'''
    cur.execute("""
SELECT * FROM antismash.rel_clusters_types WHERE bgc_id = %s AND
    bgc_type_id = (SELECT bgc_type_id FROM antismash.bgc_types WHERE term = %s)""",
                (params['bgc_id'], product))
    ret = cur.fetchone()
    if ret is None:
        print("creating link for type {}".format(product))
        cur.execute("""
INSERT INTO antismash.rel_clusters_types (bgc_id, bgc_type_id)
SELECT val.bgc_id, f.bgc_type_id FROM ( VALUES (%s, %s) ) val (bgc_id, bgc_type)
LEFT JOIN antismash.bgc_types f ON val.bgc_type = f.term
""", (params['bgc_id'], product))


def get_lineage(taxid):
    '''Get the full lineage for a taxid from Entrez'''
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    lineage = {}
    for entry in records[0]['LineageEx']:
        if entry['Rank'] == 'no rank':
            continue
        lineage[entry['Rank']] = entry['ScientificName']

    return lineage


if __name__ == "__main__":
    main()
