#!/usr/bin/env python
"""Import a GenBank results file into the antiSMASH database."""
from __future__ import print_function
from collections import defaultdict
import hashlib
import os
import sys
import typing as t
import re
from Bio import Entrez, SeqIO
import urllib
import psycopg2
import psycopg2.extensions
psycopg2.extensions.register_type(psycopg2.extensions.UNICODE)
psycopg2.extensions.register_type(psycopg2.extensions.UNICODEARRAY)

DB_CONNECTION = "host='localhost' port=5432 user='postgres' password='secret' dbname='antismash'"
Entrez.email = "kblin@biosustain.dtu.dk"
SMCOG_PATTERN = re.compile(r"smCOG: (SMCOG[\d]{4}):[\w'`:,/\s\(\)\[\]-]+\(Score: ([\d.e-]+); E-value: ([\d.e-]+)\);")
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


def main():
    """Run the import."""
    if len(sys.argv) < 2:
        print('Usage: {} <gbk file>'.format(sys.argv[0]), file=sys.stderr)
        sys.exit(2)

    connection = psycopg2.connect(DB_CONNECTION)

    recs = SeqIO.parse(sys.argv[1], 'genbank')
    with connection:
        with connection.cursor() as cursor:
            assembly_id = None
            for rec in recs:
                if rec.name in BLACKLIST:
                    print('Skipping blacklisted record {!r}'.format(rec.name), file=sys.stderr)
                    continue
                if not assembly_id:
                    assembly_id = get_assembly_id(rec)
                load_record(rec, cursor)
            if assembly_id:
                assembly_id = assembly_id.split('.')[0]
                input_basename = os.path.basename(sys.argv[1])
                if input_basename.endswith('.final.gbk'):
                    input_basename = input_basename[:-10]
                cursor.execute("INSERT INTO antismash.filenames (assembly_id, base_filename) VALUES (%s, %s)", (assembly_id, input_basename))

    connection.close()


def load_record(rec, cur):
    """Load a record into the database using the cursor."""
    genome_id = get_or_create_genome(rec, cur)
    print("genome_id: {}".format(genome_id))
    seq_id = get_or_create_dna_sequence(rec, cur, genome_id)
    print("seq_id: {}".format(seq_id))

    FEATURE_HANDLERS = {
        'CDS': handle_cds,
        'CDS_motif': handle_cds_motif,
        'aSDomain': handle_asdomain,
        'cluster': handle_cluster,
        'gene': handle_gene,
        'misc_feature': handle_misc_feature,
        'PFAM_domain': handle_pfamdomain,
    }

    for feature in rec.features:
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

    cur.execute("SELECT sequence_id FROM antismash.dna_sequences WHERE md5 = %s", (params['md5sum'],))
    ret = cur.fetchone()
    if ret is None:
        params['genome_id'] = genome_id
        cur.execute("INSERT INTO antismash.dna_sequences (dna, md5, acc, version, genome_id)"
                    "VALUES (%(seq)s, %(md5sum)s, %(accession)s, %(version)s, %(genome_id)s)"
                    "RETURNING sequence_id;", params)
        ret = cur.fetchone()

    return ret[0]


def get_or_create_genome(rec, cur):
    """Fetch existing genome entry or create a new one."""
    try:
        taxid = get_or_create_tax_id(cur, get_taxid(rec), get_strain(rec))
    except psycopg2.ProgrammingError:
        print(rec)
        raise
    assembly_id = get_assembly_id(rec)
    cur.execute("SELECT genome_id FROM antismash.genomes WHERE tax_id = %s AND assembly_id = %s", (taxid, assembly_id))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("INSERT INTO antismash.genomes (tax_id, assembly_id) VALUES (%s, %s) RETURNING genome_id;", (taxid, assembly_id))
        ret = cur.fetchone()

    return ret[0]


def get_taxid(rec):
    """Extract the taxid from a record."""
    for feature in rec.features:
        if feature.type != 'source':
            continue
        if 'db_xref' not in feature.qualifiers:
            return 0
        for entry in feature.qualifiers['db_xref']:
            if entry.startswith('taxon:'):
                return int(entry[6:])


def get_strain(rec):
    """Extract the strain from a record."""
    for feature in rec.features:
        if feature.type != 'source':
            continue
        if 'strain' in feature.qualifiers:
            return feature.qualifiers['strain'][0]
        if 'serovar' in feature.qualifiers:
            return feature.qualifiers['serovar'][0]

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
    params['start_pos'] = feature.location.nofuzzy_start
    params['end_pos'] = feature.location.nofuzzy_end
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
        params['locus_tag'] = feature.qualifiers['locus_tag'][0]
        cur.execute("INSERT INTO antismash.genes (locus_tag, locus_id) VALUES (%(locus_tag)s, %(locus_id)s)", params)


def handle_cds(rec, cur, seq_id, feature):
    """Handle CDS features."""
    params = {}
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)

    cur.execute("SELECT cds_id FROM antismash.cdss WHERE locus_id = %s", (params['locus_id'],))
    ret = cur.fetchone()
    if ret is None:
        params['locus_tag'] = feature.qualifiers['locus_tag'][0]
        params['name'] = feature.qualifiers.get('gene', [None])[0]
        params['product'] = feature.qualifiers.get('product', [None])[0]
        params['protein_id'] = feature.qualifiers.get('protein_id', [None])[0]
        params['func_class'] = get_functional_class(cur, feature)
        params['evidence'] = 'prediction'
        params['translation'] = get_translation(feature)
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
    create_smcog_hit(cur, feature, cds_id)
    create_profile_hits(cur, feature, cds_id)
    create_terpene_cyclisations(cur, feature, cds_id)


def handle_misc_feature(rec, cur, seq_id, feature):
    """Handle (a subset of) misc_feature features."""

    # only want to handle misc_features created by antiSMASH
    if 'tool' not in feature.qualifiers or 'antiSMASH'not in feature.qualifiers['tool']:
        return

    is_tta_codon = False
    for entry in feature.qualifiers['note']:
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


def create_smcog_hit(cur, feature, cds_id):
    """Create an smCOG hit entry."""
    try:
        smcog_name, smcog_score, smcog_evalue = parse_smcog(feature)
        smcog_score = float(smcog_score)
        smcog_evalue = float(smcog_evalue)
        smcog_id = get_smcog_id(cur, smcog_name)
        cur.execute("SELECT cds_id, smcog_id FROM antismash.smcog_hits WHERE "
                    "smcog_id = %s AND cds_id = %s", (smcog_id, cds_id))
        ret = cur.fetchone()
        if ret is None:
            cur.execute("INSERT INTO antismash.smcog_hits (score, evalue, smcog_id, cds_id)"
                        "VALUES (%s, %s, %s, %s)", (smcog_score, smcog_evalue, smcog_id, cds_id))
    except ValueError as e:
        # no smcog qualifier is an expected error, don't log that
        err_msg = str(e)
        if not (err_msg.startswith('No smcog qualifier') or
                err_msg.startswith('No note qualifier')):
            print(e, file=sys.stderr)


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


def create_terpene_cyclisations(cur, feature, cds_id):
    """Create terpene cyclisations entry"""
    detected = parse_terpene_cyclisation(feature)
    if not detected:
        return

    detected['cds_id'] = cds_id
    cur.execute("""
SELECT cds_id FROM antismash.terpene_cyclisations LEFT JOIN antismash.terpenes USING (terpene_id) WHERE
    cds_id = %(cds_id)s AND
    from_carbon = %(from_carbon)s AND
    to_carbon = %(to_carbon)s AND
    name = %(synthase_type)s""", detected)
    ret = cur.fetchone()

    if ret is None:
        cur.execute("""
INSERT INTO antismash.terpene_cyclisations (terpene_id, cds_id, from_carbon, to_carbon) VALUES (
    (SELECT terpene_id FROM antismash.terpenes WHERE name = %(synthase_type)s),
    %(cds_id)s,
    %(from_carbon)s,
    %(to_carbon)s
)""", detected)


def parse_terpene_cyclisation(feature):
    """Parse terpene cyclisation patterns."""
    detected = {}

    if not 'note' in feature.qualifiers:
        return detected

    for note in feature.qualifiers['note']:
        if not note.startswith('Cyclization pattern: '):
            continue
        if note.endswith('no prediction'):
            return detected

        parts = note[21:].split('-')
        if len(parts) != 3:
            return detected

        try:
            detected['from_carbon'] = int(parts[0])
            detected['to_carbon'] = int(parts[1])
            detected['synthase_type'] = parts[2]
        except ValueError:
            return {}
        break

    return detected


def parse_domains_detected(feature):
    """Parse detected domains."""
    pattern = re.compile(r'([\w-]+) \(E-value: ([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?), '
                         r'bitscore: ([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?), seeds: (\d+)')
    domains = []
    if 'sec_met' not in feature.qualifiers:
        return domains

    domain_line = None
    for secmet_line in feature.qualifiers['sec_met']:
        prefix = 'Domains detected: '
        if secmet_line.startswith(prefix):
            domain_line = secmet_line[len(prefix):]
            break

    if domain_line is None:
        return domains

    for domain in domain_line.split(';'):
        match = pattern.search(domain)
        if not match:
            print("no match found for", domain, file=sys.stderr)
            print("line:", domain_line, file=sys.stderr)
            continue
        dom = {}
        dom['name'] = match.group(1)
        dom['evalue'] = match.group(2)
        dom['bitscore'] = match.group(4)
        dom['seeds'] = match.group(6)
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


def get_functional_class(cur, feature):
    """Get the functional class from a CDS feature."""
    func_class = 'other'
    try:
        smcog, _, _ = parse_smcog(feature)
        cur.execute("""
SELECT name FROM antismash.functional_classes WHERE functional_class_id =
    (SELECT functional_class_id FROM antismash.smcogs WHERE name = %s)""", (smcog,))
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
    """Parse the smCOG feature qualifier."""
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
    """Handle CDS_motif features."""
    if 'aSTool' in feature.qualifiers:
        # This is a CDS_motif from the nrpspks module, ignore
        # We can find all info we need in the aSDomain record
        return

    if 'note' not in feature.qualifiers:
        return

    params = defaultdict(lambda: None)
    try:
        params['bgc_id'] = get_bgc_id_from_overlap(cur, seq_id, feature)
    except ValueError:
        return
    params['locus_tag'] = feature.qualifiers['locus_tag'][0]
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
    for note in feature.qualifiers['note']:
        if note.startswith('monoisotopic mass:'):
            params['monoisotopic_mass'] = note.split(':')[-1].strip()
            continue
        if note.startswith('molecular weight:'):
            params['molecular_weight'] = note.split(':')[-1].strip()
            continue
        if note.startswith('alternative weights:'):
            params['alternative_weights'] = note.split(':')[-1].strip()
            continue
        if note.startswith('number of bridges:'):
            params['bridges'] = note.split(':')[-1].strip()
            continue
        if note.startswith('predicted core seq:'):
            params['peptide_sequence'] = note.split(':')[-1].strip()
            continue
        if note.startswith('predicted class:'):
            params['class'] = note.split(':')[-1].strip()
            continue
        if note.startswith('score:'):
            params['score'] = note.split(':')[-1].strip()
            continue


def get_bgc_id_from_overlap(cur, seq_id, feature):
    """Query for bgc_ids that contain the feature."""
    start_pos = feature.location.nofuzzy_start
    end_pos = feature.location.nofuzzy_end
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
        params['score'] = float(feature.qualifiers['score'][0]) if 'score' in feature.qualifiers else None
        params['evalue'] = float(feature.qualifiers['evalue'][0]) if 'evalue' in feature.qualifiers else None
        params['translation'] = feature.qualifiers['translation'][0] if 'translation' in feature.qualifiers else None
        params['locus_tag'] = feature.qualifiers['locus_tag'][0] if 'locus_tag' in feature.qualifiers else None
        params['detection'] = feature.qualifiers['detection'][0] if 'detection' in feature.qualifiers else None
        params['as_domain_profile_id'] = get_as_domain_profile_id(cur, feature.qualifiers.get('domain', [None])[0])

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
        params['score'] = float(feature.qualifiers['score'][0]) if 'score' in feature.qualifiers else None
        params['evalue'] = float(feature.qualifiers['evalue'][0]) if 'evalue' in feature.qualifiers else None
        params['translation'] = feature.qualifiers['translation'][0] if 'translation' in feature.qualifiers else None
        params['locus_tag'] = feature.qualifiers['locus_tag'][0] if 'locus_tag' in feature.qualifiers else None
        params['detection'] = feature.qualifiers['detection'][0] if 'detection' in feature.qualifiers else None
        params['database'] = feature.qualifiers['database'][0] if 'database' in feature.qualifiers else None
        params['pfam_id'] = get_pfam_id(cur, feature.qualifiers.get('db_xref', [None])[0])

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


def get_pfam_id(cur, xref):
    """Get the pfam_id for a domain by db_xref string."""
    if xref is None:
        return None

    if not xref.startswith("PFAM: "):
        raise ValueError("Invalid db_xref trying to parse pfam_id: {!r}".format(xref))

    pfam = xref.split()[-1].strip()
    cur.execute("SELECT pfam_id FROM antismash.pfams WHERE pfam_id = %s", (pfam,))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("Invalid pfam_id {!r}".format(pfam))

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
    if 'specificity' in feature.qualifiers:
        for spec in feature.qualifiers['specificity']:
            if spec.startswith('KR activity: '):
                params['kr_activity'] = False if spec.endswith('inactive') else True
                continue
            if spec.startswith('KR stereochemistry: '):
                params['kr_stereochemistry'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('NRPSpredictor3 SVM: '):
                params['nrps_predictor'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('Stachelhaus code: '):
                params['stachelhaus'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('Minowa: '):
                params['minowa'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('pHMM: '):
                params['phmm'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('PrediCAT '):
                params['predicat'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('PID to NN: '):
                params['pid'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('SNN score: '):
                params['snn_score'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('PKS signature: '):
                params['pks_signature'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('consensus: '):
                params['consensus'] = spec.split(':')[-1].strip()
                continue
            if spec.startswith('SANDPUMA ensemble: '):
                params['consensus'] = spec.split(':')[-1].strip()
                continue



def handle_cluster(rec, cur, seq_id, feature):
    """Handle cluster features."""
    params = defaultdict(lambda: None)
    params['locus_id'] = get_or_create_locus(cur, seq_id, feature)
    params['evidence'] = 'prediction'
    print("locus_id: {}".format(params['locus_id']))

    for note in feature.qualifiers['note']:
        if note.startswith('Cluster number: '):
            params['cluster_number'] = note.split(':')[-1].strip()

    cur.execute("SELECT bgc_id FROM antismash.biosynthetic_gene_clusters WHERE locus_id = %(locus_id)s", params)
    ret = cur.fetchone()
    if ret is None:
        cur.execute("""
INSERT INTO antismash.biosynthetic_gene_clusters (cluster_number, locus_id, evidence_id)
SELECT val.cluster_number::int4, val.locus_id, f.evidence_id FROM (
    VALUES (%(cluster_number)s, %(locus_id)s, %(evidence)s) ) val (cluster_number, locus_id, evidence)
LEFT JOIN antismash.evidences f ON val.evidence = f.name
RETURNING bgc_id
""", params)
        ret = cur.fetchone()
    params['bgc_id'] = ret[0]

    for product in feature.qualifiers['product'][0].split('-'):
        product = product.lower().replace(' ', '')
        try:
            nx_create_rel_clusters_types(cur, params, product)
        except psycopg2.IntegrityError:
            print("Failed to insert product type", product)
            raise

    store_clusterblast(cur, feature, 'clusterblast', params['bgc_id'])
    store_clusterblast(cur, feature, 'knownclusterblast', params['bgc_id'])
    store_clusterblast(cur, feature, 'subclusterblast', params['bgc_id'])

    monomer_list, monomer_str = parse_monomers(feature)
    if monomer_str is None:
        # We don't have a compound prediction at this stage
        return

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

    for i, monomers in enumerate(monomer_list):
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


def parse_monomers(feature):
    """Parse a feature's monomoers string."""
    for note in feature.qualifiers['note']:
        if note.startswith('Monomers prediction: '):
            monomers_str = note.split(':')[-1].strip()
            if monomers_str == '':
                return [], None
            monomers_condensed = monomers_str.replace('(', '').replace(')', '').replace(' ', '').replace('+', '-')
            monomers = monomers_condensed.split('-')
            return [monomer.split('|') for monomer in monomers], monomers_str
    return [], None


def nx_create_rel_clusters_types(cur, params, product):
    """Create relation table to bgc_types."""
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


def get_or_create_tax_id(cur, taxid, strain):
    """Get the tax_id or create a new one."""
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

