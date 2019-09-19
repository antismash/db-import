#!/usr/bin/env python3
"""Import a GenBank results file into the antiSMASH database."""
from argparse import ArgumentParser
from collections import defaultdict
import hashlib
import os
import sys
import re
import urllib

import antismash
from Bio import Entrez
import psycopg2
import psycopg2.extensions
psycopg2.extensions.register_type(psycopg2.extensions.UNICODE)
psycopg2.extensions.register_type(psycopg2.extensions.UNICODEARRAY)

DB_CONNECTION = "host='localhost' port=5432 user='postgres' password='secret' dbname='antismash'"
Entrez.email = "kblin@biosustain.dtu.dk"
# the black list contains accessions that contain duplicate or invalid locus tags
BLACKLIST = [
    'NC_017443',
    'NC_017447',
    'NZ_AONY01000056',
    'NZ_AOOL01000046',
    'NZ_AOOR01000035',
    'NZ_CP024070',
    'NZ_KE136343',
    'NZ_AOOE01000052',
    'NZ_CP022430',
]
REPORTED_TYPES = set()

MINIMAL = False
VISIBILITY = 'public'


def main():
    """Run the import."""
    global MINIMAL
    global VISIBILITY

    parser = ArgumentParser()
    parser.add_argument('--minimal', action='store_true', default=False,
                        help="Set when importing results of a minimal/fast mode antiSMASH run")
    parser.add_argument('--visibility', choices=['public'], default='public',
                        help="Set the record visibility (default: %(default)s)")
    parser.add_argument('filename')
    args = parser.parse_args()

    if args.minimal:
        MINIMAL = True
    VISIBILITY = args.visibility

    connection = psycopg2.connect(DB_CONNECTION)

    recs = antismash.common.secmet.Record.from_genbank(args.filename)
    with connection:
        with connection.cursor() as cursor:
            assembly_id = None
            for rec in recs:
                if rec.name in BLACKLIST:
                    print('Skipping blacklisted record {!r}'.format(rec.name), file=sys.stderr)
                    continue
                if not assembly_id:
                    assembly_id = get_assembly_id(rec)
                load_record(rec, cursor, assembly_id)
            if assembly_id:
                assembly_id = assembly_id.split('.')[0]
                input_basename = os.path.basename(args.filename)
                if input_basename.endswith('.final.gbk'):
                    input_basename = input_basename[:-10]
                cursor.execute("INSERT INTO antismash.filenames (assembly_id, base_filename) VALUES (%s, %s)", (assembly_id, input_basename))

    connection.close()


def load_record(rec, cur, assembly_id):
    """Load a record into the database using the cursor."""
    genome_id = get_or_create_genome(rec, cur, assembly_id)
    print("genome_id: {}".format(genome_id))
    seq_id = get_or_create_dna_sequence(rec, cur, genome_id)
    print("seq_id: {}".format(seq_id))

    FEATURE_HANDLERS = {
        'CDS_motif': handle_cds_motif,
        'aSDomain': handle_asdomain,
        'gene': handle_gene,
        'misc_feature': handle_misc_feature,
        'PFAM_domain': handle_pfamdomain,
    }

    for region in sorted(rec.get_regions()):
        handle_region(rec, cur, seq_id, region)
        for cds in region.cds_children:
            handle_cds(rec, cur, seq_id, cds)

    for feature in rec.get_all_features():
        if feature.type in ["region", "CDS"]:
            continue
        if feature.type not in FEATURE_HANDLERS:
            if feature.type not in REPORTED_TYPES:
                print("Skipping unknown feature type", feature.type, file=sys.stderr)
                REPORTED_TYPES.add(feature.type)
            continue
        handler = FEATURE_HANDLERS[feature.type]
        handler(rec, cur, seq_id, feature)


def get_or_create_dna_sequence(rec, cur, genome_id):
    """Fetch existing dna_sequence entry or create a new one."""
    params = {}
    params['seq'] = str(rec.seq)
    params['md5sum'] = hashlib.md5(params['seq'].encode('utf-8')).hexdigest()
    params['accession'] = rec.annotations['accessions'][0]
    params['version'] = rec.annotations.get('sequence_version', '0')
    params['definition'] = rec.description

    cur.execute("SELECT sequence_id FROM antismash.dna_sequences WHERE md5 = %s", (params['md5sum'],))
    ret = cur.fetchone()
    if ret is None:
        params['genome_id'] = genome_id
        cur.execute("INSERT INTO antismash.dna_sequences (dna, md5, acc, version, definition, genome_id)"
                    "VALUES (%(seq)s, %(md5sum)s, %(accession)s, %(version)s, %(definition)s, %(genome_id)s)"
                    "RETURNING sequence_id;", params)
        ret = cur.fetchone()

    return ret[0]


def get_or_create_genome(rec, cur, assembly_id):
    """Fetch existing genome entry or create a new one."""
    try:
        taxid = get_or_create_tax_id(cur, get_taxid(rec), get_strain(rec))
    except psycopg2.ProgrammingError:
        print(rec)
        raise
    cur.execute("SELECT genome_id FROM antismash.genomes WHERE tax_id = %s AND assembly_id = %s", (taxid, assembly_id))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.genomes (tax_id, assembly_id) VALUES (%s, %s) RETURNING genome_id;", (taxid, assembly_id))
        ret = cur.fetchone()

    return ret[0]


def get_taxid(rec):
    """Extract the taxid from a record."""
    for feature in rec.get_misc_feature_by_type("source"):
        if feature.type != 'source':
            continue
        refs = feature.get_qualifier("db_xref")
        if refs is None:
            return 0
        for entry in refs:
            if entry.startswith('taxon:'):
                return int(entry[6:])


def get_strain(rec):
    """Extract the strain from a record."""
    for feature in rec.get_misc_feature_by_type("source"):
        strain = feature.get_qualifier("strain")
        if strain:
            return strain[0]
        serovar = feature.get_qualifier("serovar")
        if serovar:
            return serovar[0]

    return None


def get_assembly_id(rec):
    """Extract the NCBI assembly ID from a record."""
    for ref in rec.dbxrefs:
        if not ref.startswith('Assembly:'):
            continue
        return ref[9:]

    return None


def get_or_create_locus(cur, seq_id, feature):
    """Get or create a new locus tag."""
    params = {}
    params['start_pos'] = int(feature.location.start)
    params['end_pos'] = int(feature.location.end)
    if feature.location.strand == 1:
        params['strand'] = '+'
    elif feature.location.strand == -1:
        params['strand'] = '-'
    else:
        params['strand'] = '0'
    params['sequence_id'] = seq_id
    cur.execute("SELECT locus_id FROM antismash.loci WHERE sequence_id = %(sequence_id)s AND "
                "start_pos = %(start_pos)s AND end_pos = %(end_pos)s AND strand = %(strand)s",
                params)
    ret = cur.fetchone()

    if ret is None:
        cur.execute("INSERT INTO antismash.loci (start_pos, end_pos, strand, sequence_id) VALUES "
                    "(%(start_pos)s, %(end_pos)s, %(strand)s, %(sequence_id)s) RETURNING locus_id", params)
        ret = cur.fetchone()

    return ret[0]


def handle_gene(rec, cur, seq_id, feature):
    """Handle gene features."""
    params = {}
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT gene_id FROM antismash.genes WHERE locus_id = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        params['locus_tag'] = feature.locus_tag
        cur.execute("INSERT INTO antismash.genes (locus_tag, locus_id) VALUES (%(locus_tag)s, %(locus_id)s)", params)


def handle_cds(rec, cur, seq_id, feature):
    """Handle CDS features."""
    params = {}
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT cds_id FROM antismash.cdss WHERE locus_id = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        params['locus_tag'] = feature.locus_tag
        params['name'] = feature.gene
        params['product'] = feature.product
        params['protein_id'] = feature.protein_id
        params['func_class'] = str(feature.gene_function)
        params['evidence'] = 'prediction'
        params['translation'] = feature.translation
        cur.execute("""
INSERT INTO antismash.cdss (
    functional_class_id,
    locus_tag,
    name,
    product,
    protein_id,
    locus_id,
    translation,
    evidence_id
) VALUES
( (SELECT functional_class_id FROM antismash.functional_classes WHERE name = %(func_class)s),
  %(locus_tag)s, %(name)s, %(product)s, %(protein_id)s, %(locus_id)s, %(translation)s,
  (SELECT evidence_id FROM antismash.evidences WHERE name = %(evidence)s) ) RETURNING cds_id
""", params)
        ret = cur.fetchone()
    cds_id = ret[0]
    for gene_function in feature.gene_functions.get_by_tool("smcogs"):
        create_smcog_hit(cur, gene_function, cds_id)
        print("     Skipping all  but first SMCOG"); break  # TODO
    create_profile_hits(cur, feature, cds_id)


def handle_misc_feature(rec, cur, seq_id, feature):
    """Handle (a subset of) misc_feature features."""

    # only want to handle misc_features created by antiSMASH
    if not feature.created_by_antismash:
        return

    is_tta_codon = False
    for entry in feature.notes:
        if entry.startswith("tta leucine codon"):
            is_tta_codon = True
            break

    if not is_tta_codon:
        return

    params = {}
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT tta_codon_id FROM antismash.tta_codons WHERE locus_id = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.tta_codons (locus_id) VALUES ( %(locus_id)s )", params)


def create_smcog_hit(cur, gene_function, cds_id):
    """Create an smCOG hit entry."""
    assert gene_function.tool == "smcogs"
    smcog_name = gene_function.product
    smcog_score = 0  # TODO
    smcog_evalue = 100.0  # TODO
    smcog_id = get_smcog_id(cur, smcog_name)
    cur.execute("SELECT cds_id, smcog_id FROM antismash.smcog_hits WHERE "
                "smcog_id = %s AND cds_id = %s", (smcog_id, cds_id))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.smcog_hits (score, evalue, smcog_id, cds_id)"
                    "VALUES (%s, %s, %s, %s)", (smcog_score, smcog_evalue, smcog_id, cds_id))


def create_profile_hits(cur, feature, cds_id):
    """Create profile hit entries for a feature."""
    detected_domains = parse_domains_detected(feature)
    for domain in detected_domains:
        domain['cds_id'] = cds_id
        cur.execute("""
SELECT cds_id FROM antismash.profile_hits WHERE
    cds_id = %(cds_id)s AND
    name = %(name)s AND
    evalue = %(evalue)s AND
    bitscore = %(bitscore)s""", domain)
        ret = cur.fetchone()
        if ret is None:
            try:
                cur.execute("""
INSERT INTO antismash.profile_hits (cds_id, name, evalue, bitscore, seeds)
    VALUES (%(cds_id)s, %(name)s, %(evalue)s, %(bitscore)s, %(seeds)s)""", domain)
            except psycopg2.IntegrityError:
                print(feature)
                print(domain)
                raise


def parse_domains_detected(feature):
    """Parse detected domains."""
    domains = []
    if not feature.sec_met:
        return domains

    for domain in feature.sec_met.domains:
        dom = {}
        dom['name'] = domain.name
        dom['evalue'] = domain.evalue
        dom['bitscore'] = domain.bitscore
        dom['seeds'] = domain.nseeds
        domains.append(dom)

    return domains


def get_translation(feature):
    """Get the protein translation from the feature."""
    if 'translation' not in feature.qualifiers:
        return None

    return feature.qualifiers['translation'][0]


def get_smcog_id(cur, smcog):
    """Get the smcog_id given the smCOG identifier in a feature."""
    cur.execute("SELECT smcog_id FROM antismash.smcogs WHERE name = %s", (smcog,))
    ret = cur.fetchone()
    if ret is None:
        return None
    return ret[0]


def handle_cds_motif(rec, cur, seq_id, feature):
    """Handle CDS_motif features."""
    if not isinstance(feature, antismash.common.secmet.Prepeptide):
        # This is a CDS_motif from the nrpspks module, ignore
        # We can find all info we need in the aSDomain record
        return

    params = defaultdict(lambda: None)
    for region in rec.get_regions():
        if region.overlaps_with(feature):
            params['bgc_id'] = region.get_region_number()
            break
    assert "bgc_id" in params
    params['locus_tag'] = feature.locus_tag
    parse_ripp_core(feature, params)
    if params['peptide_sequence'] is None:
        return

    cur.execute("SELECT compound_id FROM antismash.compounds WHERE "
                "peptide_sequence = %(peptide_sequence)s AND locus_tag = %(locus_tag)s", params)
    ret = cur.fetchone()
    if ret is None:
        cur.execute("""
INSERT INTO antismash.compounds (
    peptide_sequence,
    molecular_weight,
    monoisotopic_mass,
    alternative_weights,
    bridges,
    class,
    score,
    locus_tag
) VALUES (
    %(peptide_sequence)s,
    %(molecular_weight)s,
    %(monoisotopic_mass)s,
    %(alternative_weights)s,
    %(bridges)s,
    %(class)s,
    %(score)s,
    %(locus_tag)s
) RETURNING compound_id""", params)
        ret = cur.fetchone()

    compound_id = ret[0]

    cur.execute("SELECT bgc_id FROM antismash.rel_clusters_compounds WHERE "
                "bgc_id = %s AND compound_id = %s", (params['bgc_id'], compound_id))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.rel_clusters_compounds (bgc_id, compound_id)"
                    "VALUES (%s, %s)", (params['bgc_id'], compound_id))


def parse_ripp_core(feature, params):
    """Parse RIPP core features."""
    params['monoisotopic_mass'] = feature.monoisotopic_mass
    params['molecular_weight'] = feature.molecular_weight
    params['alternative_weights'] = str(feature.alternative_weights)
    params['class'] = feature.peptide_class
    params['score'] = feature.score
    params['peptide_sequence'] = feature.core
    if feature.detailed_information:
        bridges = feature.detailed_information.to_biopython_qualifiers().get("number_of_bridges")
        if bridges is not None:
            params["bridges"] = bridges


def get_bgc_id_from_overlap(cur, seq_id, feature):
    """Query for bgc_ids that contain the feature."""
    start_pos = int(feature.location.start)
    end_pos = int(feature.location.end)
    cur.execute("""
SELECT bgc.bgc_id FROM antismash.biosynthetic_gene_clusters bgc JOIN antismash.loci l USING (locus_id)
    WHERE l.sequence_id = %s AND int4range(l.start_pos, l.end_pos) @> int4range(%s, %s)""",
                (seq_id, start_pos, end_pos))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError('No bgc found overlapping {}'.format(feature))

    return ret[0]


def handle_asdomain(rec, cur, seq_id, feature):
    """Handle aSDomain features."""
    params = {}
    params['pks_signature'] = None
    params['minowa'] = None
    params['nrps_predictor'] = None
    params['stachelhaus'] = None
    params['phmm'] = None
    params['predicat'] = None
    params['pid'] = None
    params['snn_score'] = None
    params['consensus'] = None
    params['kr_activity'] = None
    params['kr_stereochemistry'] = None
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT as_domain_id FROM antismash.as_domains WHERE locus_id = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        params['score'] = feature.score
        params['evalue'] = feature.evalue
        params['translation'] = feature.translation
        params['locus_tag'] = feature.locus_tag
        params['detection'] = feature.detection
        params['as_domain_profile_id'] = get_as_domain_profile_id(cur, feature.domain_subtype or feature.domain)

        parse_specificity(feature, params)

        try:
            cur.execute("""
INSERT INTO antismash.as_domains (
    detection,
    score,
    evalue,
    translation,
    pks_signature,
    minowa,
    nrps_predictor,
    stachelhaus,
    phmm,
    predicat,
    pid,
    snn_score,
    consensus,
    kr_activity,
    kr_stereochemistry,
    as_domain_profile_id,
    locus_id,
    cds_id
) VALUES (
    %(detection)s,
    %(score)s,
    %(evalue)s,
    %(translation)s,
    %(pks_signature)s,
    %(minowa)s,
    %(nrps_predictor)s,
    %(stachelhaus)s,
    %(phmm)s,
    %(predicat)s,
    %(pid)s,
    %(snn_score)s,
    %(consensus)s,
    %(kr_activity)s,
    %(kr_stereochemistry)s,
    %(as_domain_profile_id)s,
    %(locus_id)s,
    (SELECT cds_id FROM antismash.cdss WHERE locus_tag = %(locus_tag)s)
) RETURNING as_domain_id""", params)
        except psycopg2.ProgrammingError:
            print("error fetching cds_id for locus_tag", params['locus_tag'])
            raise
        as_domain_id = cur.fetchone()[0]
        if params['consensus'] is None:
            return

        for monomer in params['consensus'].split('|'):
            monomer_id = get_or_create_monomer(cur, monomer)
            cur.execute("SELECT as_domain_id FROM antismash.rel_as_domains_monomers WHERE "
                        "as_domain_id = %s AND monomer_id = %s", (as_domain_id, monomer_id))
            ret = cur.fetchone()
            if ret is None:
                cur.execute("INSERT INTO antismash.rel_as_domains_monomers (as_domain_id, monomer_id) VALUES "
                            "(%s, %s)", (as_domain_id, monomer_id))


def handle_pfamdomain(rec, cur, seq_id, feature):
    """Handle PFAM_domain features."""
    params = {}

    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT pfam_domain_id FROM antismash.pfam_domains WHERE locus_id = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        params['score'] = feature.score
        params['evalue'] = feature.evalue
        params['translation'] = feature.translation
        params['locus_tag'] = feature.locus_tag
        params['detection'] = feature.detection
        params['database'] = feature.database
        params['pfam_id'] = get_pfam_id(cur, feature.identifier)

        cur.execute("""
INSERT INTO antismash.pfam_domains (
    database,
    detection,
    score,
    evalue,
    translation,
    pfam_id,
    locus_id,
    cds_id
) VALUES (
    %(database)s,
    %(detection)s,
    %(score)s,
    %(evalue)s,
    %(translation)s,
    %(pfam_id)s,
    %(locus_id)s,
    (SELECT cds_id FROM antismash.cdss WHERE locus_tag = %(locus_tag)s)
)""", params)


def get_pfam_id(cur, identifier):
    """Get the pfam_id for a domain by db_xref string."""
    assert identifier
    cur.execute("SELECT pfam_id FROM antismash.pfams WHERE pfam_id = %s", (identifier,))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("Invalid pfam_id {!r}".format(identifier))

    return ret[0]


def get_as_domain_profile_id(cur, name):
    """Get the as_domain_profile_id for a domain by name."""
    if name is None:
        return None

    cur.execute("SELECT as_domain_profile_id FROM antismash.as_domain_profiles WHERE name = %s", (name, ))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError('Invalid asDomain name {!r}'.format(name))

    return ret[0]


def get_or_create_monomer(cur, name):
    """Get the monomer_id for a monomer by name, create monomer entry if needed."""
    cur.execute("SELECT monomer_id FROM antismash.monomers WHERE name = %s", (name,))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.monomers (name) VALUES (%s) RETURNING monomer_id", (name,))
        ret = cur.fetchone()

    return ret[0]


def parse_specificity(feature, params):
    """Parse the feature's specificity entries."""
    for prediction in feature.specificity:
        method, pred = prediction.split(": ", 1)
        if method == 'KR activity':
            params['kr_activity'] = not pred.endswith('inactive')
        elif method == 'KR stereochemistry':
            params['kr_stereochemistry'] = pred
        elif method == 'NRPSpredictor':
            params['nrps_predictor'] = pred
        elif method == 'Stachelhaus':
            params['stachelhaus'] = pred
        elif method == 'Minowa':
            params['minowa'] = pred
        elif method == 'PKS signature':
            params['pks_signature'] = pred
        elif method == 'consensus':
            params['consensus'] = pred
        else:
            raise ValueError("unknown method: %s")


def handle_region(rec, cur, seq_id, feature):
    """Handle cluster features."""
    params = defaultdict(lambda: None)
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)
    params['evidence'] = 'prediction'
    print("locus_id: {}".format(params['locus_id']))
    params['contig_edge'] = feature.contig_edge
    params['minimal'] = MINIMAL
    params['visibility'] = VISIBILITY

    params["cluster_number"] = feature.get_region_number()

    cur.execute("SELECT bgc_id FROM antismash.biosynthetic_gene_clusters WHERE locus_id = %(locus_id)s", params)
    ret = cur.fetchone()
    if ret is None:
        cur.execute("""
INSERT INTO antismash.biosynthetic_gene_clusters (cluster_number, locus_id, evidence_id, visibility_id, contig_edge, minimal)
SELECT val.cluster_number::int4, val.locus_id, f.evidence_id, v.visibility_id, val.contig_edge, val.minimal FROM (
    VALUES (%(cluster_number)s, %(locus_id)s, %(evidence)s, %(visibility)s, %(contig_edge)s, %(minimal)s) ) val
    (cluster_number, locus_id, evidence, visibility, contig_edge, minimal)
LEFT JOIN antismash.evidences f ON val.evidence = f.name
LEFT JOIN antismash.visibilities v ON val.visibility = v.name
RETURNING bgc_id
""", params)
        ret = cur.fetchone()
    params['bgc_id'] = ret[0]
    assert params['bgc_id']

    for product in feature.products:
        product = product.lower()
        try:
            nx_create_rel_clusters_types(cur, params, product)
        except psycopg2.IntegrityError:
            raise RuntimeError("no definition in schema for product type: %s", product)

    print("         skipping clusterblast parts")  # TODO
#    store_clusterblast(cur, feature, 'clusterblast', params['bgc_id'])
#    store_clusterblast(cur, feature, 'knownclusterblast', params['bgc_id'])
#    store_clusterblast(cur, feature, 'subclusterblast', params['bgc_id'])

    flat_monomers = []
    nested_monomers = []
    for cds in feature.cds_children:
        cds_monomers = []
        for module in cds.modules:
            assert not module.is_multigene_module()
            monomers = [monomer for _, monomer in module.monomers]
            cds_monomers.extend(monomers)
            flat_monomers.append(monomers)
        if cds_monomers:
            nested_monomers.append(cds_monomers)

    if not flat_monomers:
        # We don't have a compound prediction at this stage
        return

    monomer_str = " + ".join("(%s)" % " - ".join(cds_monomers) for cds_monomers in nested_monomers)

    cur.execute("SELECT compound_id FROM antismash.compounds WHERE peptide_sequence = %s", (monomer_str,))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.compounds (peptide_sequence) VALUES (%s) RETURNING compound_id",
                    (monomer_str,))
        ret = cur.fetchone()
    compound_id = ret[0]

    cur.execute("SELECT bgc_id FROM antismash.rel_clusters_compounds WHERE bgc_id = %s AND compound_id = %s",
                (params['bgc_id'], compound_id))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.rel_clusters_compounds (bgc_id, compound_id) VALUES (%s, %s)",
                    (params['bgc_id'], compound_id))

    for i, monomers in enumerate(flat_monomers):
        position = i + 1
        for monomer in monomers:
            monomer_id = get_or_create_monomer(cur, monomer)
            cur.execute("SELECT position FROM antismash.rel_compounds_monomers WHERE "
                        "compound_id = %s AND monomer_id = %s AND position = %s",
                        (compound_id, monomer_id, position))
            ret = cur.fetchone()
            if ret is None:
                cur.execute("""
INSERT INTO antismash.rel_compounds_monomers (compound_id, monomer_id, position)
    VALUES (%s, %s, %s)""", (compound_id, monomer_id, position))


def store_clusterblast(cur, feature, algorithm, bgc_id):
    """Store XClusterBlast results in the database."""
    if algorithm not in feature.qualifiers:
        return

    cur.execute("SELECT algorithm_id FROM antismash.clusterblast_algorithms WHERE name = %s", (algorithm, ))
    ret = cur.fetchone()
    if ret is None:
        print('Did not find algorithm_id for {!r}!'.format(algorithm))
        return

    algorithm_id = ret[0]
    for entry in feature.qualifiers[algorithm]:
        params = parse_clusterblast_line(entry)
        params['algorithm_id'] = algorithm_id
        params['bgc_id'] = bgc_id
        cur.execute("SELECT clusterblast_hit_id FROM antismash.clusterblast_hits WHERE "
                    "bgc_id = %(bgc_id)s AND algorithm_id = %(algorithm_id)s AND acc = %(acc)s",
                    params)
        ret = cur.fetchone()
        if ret is None:
            cur.execute("""
INSERT INTO antismash.clusterblast_hits (rank, acc, description, similarity, algorithm_id, bgc_id) VALUES
    (%(rank)s, %(acc)s, %(description)s, %(similarity)s, %(algorithm_id)s, %(bgc_id)s)""", params)


def parse_clusterblast_line(line):
    """Parse a clusterblast result line."""
    pattern = r'''([\d]+)\. ([\w]+)\s+([\w,.:;='/&\#\+\* \(\)\[\]-]+)\(([\d]+)%'''
    m = re.search(pattern, line)
    if m is None:
        raise ValueError(line)

    rank = int(m.group(1))
    acc = m.group(2)
    desc = m.group(3)
    similarity = int(m.group(4))

    # Ugh, this is ugly, but only the MIBiG results have underscores that cause random
    # mid-word line breaks that Biopython converts into spaces.
    if acc.startswith('BGC'):
        desc = desc.replace(' ', '')
    desc = desc.replace('_', ' ').strip()
    if desc.endswith('biosynthetic gene cluster'):
        desc = desc[:-len('biosynthetic gene cluster')].strip()
    res = {
        'rank': rank,
        'acc': acc,
        'description': desc,
        'similarity': similarity
    }
    return res


def nx_create_rel_clusters_types(cur, params, product):
    """Create relation table to bgc_types."""
    assert params.get("bgc_id")
    assert product
    cur.execute("""
SELECT * FROM antismash.rel_clusters_types WHERE bgc_id = %s AND
    bgc_type_id = (SELECT bgc_type_id FROM antismash.bgc_types WHERE term = %s)""",
                (params['bgc_id'], product))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("""
INSERT INTO antismash.rel_clusters_types (bgc_id, bgc_type_id)
SELECT val.bgc_id, f.bgc_type_id FROM ( VALUES (%s, %s) ) val (bgc_id, bgc_type)
LEFT JOIN antismash.bgc_types f ON val.bgc_type = f.term""",
                    (params['bgc_id'], product))


def get_or_create_tax_id(cur, taxid, strain):
    """Get the tax_id or create a new one."""
    if taxid == 0:
        return 0
    cur.execute("SELECT tax_id FROM antismash.taxa WHERE tax_id = %s", (taxid, ))
    ret = cur.fetchone()
    if ret is None:
        lineage = get_lineage(taxid)
        lineage['tax_id'] = taxid
        lineage['strain'] = strain
        try:
            cur.execute("""
INSERT INTO antismash.taxa (tax_id, superkingdom, phylum, class, taxonomic_order, family, genus, species, strain) VALUES
    (%(tax_id)s, %(superkingdom)s, %(phylum)s, %(class)s, %(order)s, %(family)s, %(genus)s, %(species)s, %(strain)s)""",
                        lineage)
        except KeyError:
            print('Error inserting {!r}'.format(lineage))
            raise
    return taxid


def get_lineage(taxid):
    """Get the full lineage for a taxid from Entrez."""
    api_key = os.environ.get('ASDBI_ENTREZ_API_KEY', '')
    extra_params = {}
    if api_key:
        extra_params['api_key'] = api_key
    retries = 5
    while retries > 0:
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml", **extra_params)
            break
        except urllib.error.HTTPError as err:
            retries -= 1
            print("Got http error:", err, retries, "left")

    records = Entrez.read(handle)
    lineage = defaultdict(lambda: 'Unclassified')
    for entry in records[0]['LineageEx']:
        if entry['Rank'] == 'no rank':
            continue
        lineage[entry['Rank']] = entry['ScientificName'].split(' ')[-1]

    if 'species' not in lineage and 'ScientificName' in records[0] and len(records[0]['ScientificName'].split()) == 2:
        lineage['species'] = records[0]['ScientificName'].split()[-1]
    return lineage


if __name__ == "__main__":
    main()
