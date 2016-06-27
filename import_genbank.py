#!/usr/bin/env python
"""Import a GenBank results file into the antiSMASH database"""
from __future__ import print_function
import hashlib
import sys
import re
from Bio import Entrez, SeqIO
import psycopg2
import psycopg2.extensions
psycopg2.extensions.register_type(psycopg2.extensions.UNICODE)
psycopg2.extensions.register_type(psycopg2.extensions.UNICODEARRAY)

DB_CONNECTION = "host='localhost' port=15432 user='postgres' password='secret' dbname='antismash'"
Entrez.email = "kblin@biosustain.dtu.dk"
SMCOG_PATTERN = re.compile(r"smCOG: (SMCOG[\d]{4}):[\w',/\s\(\)-]+\(Score: ([\d.e-]+); E-value: ([\d.e-]+)\);")


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
    params = {}
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT gene_id FROM antismash.genes WHERE locus = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        params['locus_tag'] = feature.qualifiers['locus_tag'][0]
        params['func_class'] = get_functional_class(cur, feature)
        params['evidence'] = 'prediction'
        params['smcog_hit_id'] = None
        params['translation'] = get_translation(feature)
        try:
            smcog_name, smcog_score, smcog_evalue = parse_smcog(feature)
            smcog_score = float(smcog_score)
            smcog_evalue = float(smcog_evalue)
            smcog_id = get_smcog_id(cur, smcog_name)
            cur.execute("INSERT INTO antismash.smcog_hits (score, evalue, smcog_class) VALUES (%s, %s, %s) RETURNING smcog_hit_id", (smcog_score, smcog_evalue, smcog_id))
            params['smcog_hit_id'] = cur.fetchone()[0]
        except ValueError:
            pass

        cur.execute("""
INSERT INTO antismash.genes (functional_class, locus_tag, locus, smcog_hit, translation, evidence) VALUES
( (SELECT functional_class_id FROM antismash.functional_classes WHERE name = %(func_class)s),
  %(locus_tag)s, %(locus_id)s, %(smcog_hit_id)s, %(translation)s,
  (SELECT evidence_id FROM antismash.evidences WHERE name = %(evidence)s) )
""", params)


def get_translation(feature):
    '''Get the protein translation from the feature'''
    if 'translation' not in feature.qualifiers:
        return None

    return feature.qualifiers['translation'][0]


def get_smcog_id(cur, smcog):
    '''Get the smcog_id given the smCOG identifier in a feature'''
    cur.execute("SELECT smcog_id FROM antismash.smcogs WHERE name = %s", (smcog,))
    ret = cur.fetchone()
    if ret is None:
        return None
    return ret[0]


def get_functional_class(cur, feature):
    '''Get the functional class from a CDS feature'''
    func_class = 'other'
    try:
        smcog, _, _ = parse_smcog(feature)
        cur.execute("""
SELECT name FROM antismash.functional_classes WHERE functional_class_id =
    (SELECT functional_class FROM antismash.smcogs WHERE name = %s)""", (smcog,))
        ret = cur.fetchone()
        if ret is not None:
            func_class = ret[0]
    except ValueError:
        pass

    if 'sec_met' not in feature.qualifiers:
        return func_class

    for entry in feature.qualifiers['sec_met']:
        if entry.startswith('Kind: '):
            func_class = entry[6:]
            if func_class == 'biosynthetic':
                func_class = 'bgc_seed'
                return func_class
            break

    return func_class


def parse_smcog(feature):
    '''Parse the smCOG feature qualifier'''
    if 'note' not in feature.qualifiers:
        raise ValueError('No note qualifier in {}'.format(feature))

    for entry in feature.qualifiers['note']:
        if not entry.startswith('smCOG:'):
            continue
        match = SMCOG_PATTERN.search(entry)
        if match is None:
            print(entry)
            raise ValueError('Failed to parse smCOG line {!r}'.format(entry))
        return match.groups()

    raise ValueError('No smcog qualifier in {}'.format(feature))


def handle_cds_motif(rec, cur, seq_id, feature):
    '''Handle CDS_motif features'''
    pass


def handle_asdomain(rec, cur, seq_id, feature):
    '''Handle aSDomain features'''
    params = {}
    params['pks_signature'] = None
    params['minowa'] = None
    params['nrps_predictor'] = None
    params['stachelhaus'] = None
    params['consensus'] = None
    params['kr_activity'] = None
    params['kr_stereochemistry'] = None
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT as_domain_id FROM antismash.as_domains WHERE locus = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        params['name'] = feature.qualifiers['domain'][0] if 'domain' in feature.qualifiers else None
        params['database'] = feature.qualifiers['database'][0] if 'database' in feature.qualifiers else None
        params['score'] = float(feature.qualifiers['score'][0]) if 'score' in feature.qualifiers else None
        params['evalue'] = float(feature.qualifiers['evalue'][0]) if 'evalue' in feature.qualifiers else None
        params['translation'] = feature.qualifiers['translation'][0] if 'translation' in feature.qualifiers else None
        params['locus_tag'] = feature.qualifiers['locus_tag'][0] if 'locus_tag' in feature.qualifiers else None
        params['detection'] = feature.qualifiers['detection'][0] if 'detection' in feature.qualifiers else None

        parse_specificity(feature, params)

        cur.execute("""
INSERT INTO antismash.as_domains (
    name,
    database,
    detection,
    score,
    evalue,
    translation,
    pks_signature,
    minowa,
    nrps_predictor,
    stachelhaus,
    consensus,
    kr_activity,
    kr_stereochemistry,
    locus,
    gene
) VALUES (
    %(name)s,
    %(database)s,
    %(detection)s,
    %(score)s,
    %(evalue)s,
    %(translation)s,
    %(pks_signature)s,
    %(minowa)s,
    %(nrps_predictor)s,
    %(stachelhaus)s,
    %(consensus)s,
    %(kr_activity)s,
    %(kr_stereochemistry)s,
    %(locus_id)s,
    (SELECT gene_id FROM antismash.genes WHERE locus_tag = %(locus_tag)s)
) RETURNING as_domain_id""", params)
        as_domain_id = cur.fetchone()[0]
        if params['consensus'] is not None:
            monomer_id = get_or_create_monomer(cur, params['consensus'])
            cur.execute("SELECT as_domain_id FROM antismash.rel_as_domains_monomers WHERE as_domain_id = %s AND monomer_id = %s", (as_domain_id, monomer_id))
            ret = cur.fetchone()
            if ret is None:
                cur.execute("INSERT INTO antismash.rel_as_domains_monomers (as_domain_id, monomer_id) VALUES (%s, %s)", (as_domain_id, monomer_id))


def get_or_create_monomer(cur, name):
    '''get the monomer_id for a monomer by name, create monomer entry if needed'''
    cur.execute("SELECT monomer_id FROM antismash.monomers WHERE name = %s", (name,))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.monomers (name) VALUES (%s) RETURNING monomer_id", (name,))
        ret = cur.fetchone()

    return ret[0]


def parse_specificity(feature, params):
    '''parse the feature's specificity entries'''
    if 'specificity' in feature.qualifiers:
        for spec in feature.qualifiers['specificity']:
            if spec.startswith('KR activity: '):
                params['kr_activity'] = False if spec.endswith('inactive') else True
                continue
            if spec.startswith('KR stereochemistry: '):
                params['kr_stereochemistry'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('NRPSpredictor2 SVM: '):
                params['nrps_predictor'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('Stachelhaus code: '):
                params['stachelhaus'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('Minowa: '):
                params['minowa'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('PKS signature: '):
                params['pks_signature'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('consensus: '):
                params['consensus'] = spec.split(':')[-1].strip()
                continue


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


    for product in feature.qualifiers['product'][0].split('-'):
        nx_create_rel_clusters_types(cur, params, product)


def nx_create_rel_clusters_types(cur, params, product):
    '''create relation table to bgc_types'''
    cur.execute("""
SELECT * FROM antismash.rel_clusters_types WHERE bgc_id = %s AND
    bgc_type_id = (SELECT bgc_type_id FROM antismash.bgc_types WHERE term = %s)""",
                (params['bgc_id'], product))
    ret = cur.fetchone()
    if ret is None:
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
