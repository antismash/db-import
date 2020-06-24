#!/usr/bin/env python3
"""Import a GenBank results file into the antiSMASH database."""
from argparse import ArgumentParser
from collections import defaultdict
import hashlib
import os
import sys
import time
import urllib

# pylint: disable=line-too-long,missing-docstring

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

DEFAULT_AS_OPTIONS = antismash.config.build_config(["--minimal"], modules=antismash.main.get_all_modules())


class ExistingRecordError(ValueError):
    pass

class MissingAssemblyIdError(ValueError):
    pass

class RecordData:
    def __init__(self, cursor, record, record_id, assembly_id, module_results, record_no):
        self.cursor = cursor
        self.record = record
        self.record_id = record_id
        assert record_id
        self.assembly_id = assembly_id
        self.module_results = module_results
        self.record_no = record_no

        self._current_region = None
        self._current_region_id = None
        self.feature_mapping = {}

    @property
    def current_region(self):
        assert self._current_region
        return self._current_region

    @current_region.setter
    def current_region(self, region):
        assert isinstance(region, antismash.common.secmet.Region)
        self._current_region = region
        self._current_region_id = self.feature_mapping[region]

    @property
    def current_region_id(self):
        assert self._current_region_id
        return self._current_region_id

    def insert(self, statement, values):
        self.cursor.execute(statement, values)
        if "RETURNING" in statement:
            return self.cursor.fetchone()[0]
        return None


def get_return_id(cur):
    ret = cur.fetchone()
    if not ret:
        raise ValueError("no id to return")
    return ret[0]


def main(filename, db_connection):
    """Run the import."""
    connection = psycopg2.connect(db_connection)
    connection.autocommit = False

    results = antismash.common.serialiser.AntismashResults.from_file(filename)
    with connection.cursor() as cursor:
        try:
            assembly_id = get_assembly_id(results.records[0])
            if not assembly_id:
                short_name = os.path.basename(filename)
                id_parts = short_name.split("_")
                if id_parts[0] not in ("GCF", "GCA"):
                    raise MissingAssemblyIdError()
                assembly_id = "_".join(id_parts[:2])

            print("assembly_id:", assembly_id, end="\t")
            if assembly_id:
                input_basename = os.path.basename(filename)
                if input_basename.endswith('.final.gbk'):
                    input_basename = input_basename[:-10]
                cursor.execute("SELECT (assembly_id) FROM antismash.filenames WHERE base_filename = %s AND assembly_id = %s", (input_basename, assembly_id))
                if cursor.fetchone() is not None:
                    print("skipping previously processed file/assembly")
                    raise ExistingRecordError()
                cursor.execute("INSERT INTO antismash.filenames (assembly_id, base_filename) VALUES (%s, %s)", (assembly_id, input_basename))
            record_no = 0
            for rec, module_results in zip(results.records, results.results):
                record_no += 1
                if rec.name in BLACKLIST:
                    print('Skipping blacklisted record {!r}'.format(rec.name), file=sys.stderr)
                    continue
                prepare_record(rec, module_results)
                load_record(rec, module_results, cursor, assembly_id, record_no)
            connection.commit()
            print("changes committed")
        except ExistingRecordError:
            connection.rollback()
            print("no changes committed")
        except Exception:
            connection.rollback()
            print("no changes committed")
            raise
    connection.close()


def prepare_record(record, module_results):
    record.strip_antismash_annotations()
    if antismash.detection.hmm_detection.__name__ not in module_results:
        return
    # reannotate in the order that antismash does, for correct regions etc
    antismash.main.run_detection(record, DEFAULT_AS_OPTIONS, module_results)

    def regen(raw, module):
        assert raw
        regenerated = module.regenerate_previous_results(raw, record, DEFAULT_AS_OPTIONS)
        assert regenerated is not None, "%s results failed to generate for %s" % (module.__name__, record.id)
        regenerated.add_to_record(record)
        return regenerated

    for module in antismash.main.get_analysis_modules():
        if module.__name__ in module_results:
            module_results[module.__name__] = regen(module_results[module.__name__], module)

    for val in module_results.values():
        assert not isinstance(val, dict)


def load_record(rec, module_results, cur, assembly_id, record_no):
    """Load a record into the database using the cursor."""
    if not rec.get_regions():
        return
    genome_id = get_or_create_genome(rec, cur, assembly_id)
    print("genome_id: {}".format(genome_id))
    try:
        seq_id = get_or_create_dna_sequence(rec, cur, genome_id, record_no)
    except ExistingRecordError:
        print("skipping existing record:", rec.id)
        raise
    print("seq_id: {}".format(seq_id))

    data = RecordData(cur, rec, seq_id, assembly_id, module_results, record_no)

    for region in sorted(rec.get_regions()):
        handle_region(data, seq_id, region)
        handle_ripps(data)

    add_tta_codons(data)

    for pfam in rec.get_pfam_domains():
        handle_pfamdomain(data, pfam)
    for gene in rec.get_genes():
        handle_gene(data, gene)


def get_or_create_dna_sequence(rec, cur, genome_id, record_no):
    """Fetch existing dna_sequence entry or create a new one."""
    # record level entry
    params = {}
    params['seq'] = str(rec.seq)
    params['md5sum'] = hashlib.md5(params['seq'].encode('utf-8')).hexdigest()
    params['accession'] = rec.annotations['accessions'][0]
    params['version'] = rec.annotations.get('sequence_version', '0')
    params['definition'] = rec.description
    params['genome_id'] = genome_id
    params['record_number'] = record_no

    cur.execute("INSERT INTO antismash.dna_sequences (dna, md5, accession, version, definition, genome_id, record_number)"
                "VALUES (%(seq)s, %(md5sum)s, %(accession)s, %(version)s, %(definition)s, %(genome_id)s, %(record_number)s)"
                "RETURNING accession;", params)
    return cur.fetchone()[0]


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
        cur.execute("INSERT INTO antismash.genomes (tax_id, assembly_id) VALUES (%s, %s) RETURNING genome_id, tax_id, assembly_id;", (taxid, assembly_id))
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
    return None


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


def handle_gene(data, gene):
    """Handle gene features."""
    region_id = None
    for region in data.record.get_regions():
        if gene.is_contained_by(region):
            region_id = data.feature_mapping[region]
            break
    if not region_id:
        return

    params = {
        'location': str(gene.location),
        'locus_tag': gene.locus_tag,
        'region_id': region_id,
    }

    data.insert("""
INSERT INTO antismash.genes (locus_tag, location, region_id)
VALUES (%(locus_tag)s, %(location)s, %(region_id)s)""", params)


def handle_cds(data, region_id, cds):
    """Handle CDS features."""
    params = {}
    params['location'] = str(cds.location)

    params['locus_tag'] = cds.locus_tag
    params['name'] = cds.gene
    params['product'] = cds.product
    params['protein_id'] = cds.protein_id
    params['func_class'] = str(cds.gene_function)
    params['translation'] = cds.translation
    params['region_id'] = region_id
    cds_id = data.insert("""
INSERT INTO antismash.cdss (
    functional_class_id,
    locus_tag,
    name,
    product,
    protein_id,
    location,
    translation,
    region_id
) VALUES
( (SELECT functional_class_id FROM antismash.functional_classes WHERE name = %(func_class)s),
  %(locus_tag)s, %(name)s, %(product)s, %(protein_id)s, %(location)s, %(translation)s, %(region_id)s)
RETURNING cds_id
""", params)
    assert cds_id
    data.feature_mapping[cds] = cds_id

    genefunctions = data.module_results[antismash.detection.genefunctions.__name__]._tools

    all_smcog_results = genefunctions.get("smcogs")
    hit = all_smcog_results.best_hits.get(cds.get_name())
    if hit:
        create_smcog_hit(data.cursor, hit, cds_id)

    all_resfam_results = genefunctions.get("resist")
    if all_resfam_results:
        hit = all_resfam_results.best_hits.get(cds.get_name())
        if hit:
            create_resfam_hit(data.cursor, hit, cds_id)

    create_profile_hits(data, cds)


def add_tta_codons(data):
    tta_results = data.module_results.get(antismash.modules.tta.__name__)
    if not tta_results:
        return
    for feature in tta_results.features:
        data.insert("INSERT INTO antismash.tta_codons (location, seq_id) VALUES (%s, %s)", (str(feature.location), data.record_id))


def create_resfam_hit(cursor, hit, cds_id):
    cursor.execute("SELECT resfam_id from antismash.resfams WHERE name = %s", (hit.hit_id,))
    ret = cursor.fetchone()
    if not ret:
        raise ValueError("unknown resfam ID: %s" % hit.hit_id)
    resfam_id = ret[0]
    cursor.execute("INSERT INTO antismash.resfam_domains (score, evalue, resfam_id, cds_id) VALUES (%s, %s, %s, %s)", (hit.bitscore, hit.evalue, resfam_id, cds_id))


def create_smcog_hit(cur, hit, cds_id):
    """Create an smCOG hit entry."""
    smcog_name = hit.hit_id.split(":", 1)[0]
    smcog_score = hit.bitscore
    smcog_evalue = hit.evalue
    smcog_id = get_smcog_id(cur, smcog_name)
    cur.execute("INSERT INTO antismash.smcog_hits (score, evalue, smcog_id, cds_id)"
                "VALUES (%s, %s, %s, %s)", (smcog_score, smcog_evalue, smcog_id, cds_id))


def create_profile_hits(data, cds):
    """Create profile hit entries for a feature."""
    detected_domains = parse_domains_detected(cds)
    for domain in detected_domains:
        domain['cds_id'] = data.feature_mapping[cds]
        try:
            data.insert("""
INSERT INTO antismash.profile_hits (cds_id, name, evalue, bitscore, seeds)
VALUES (%(cds_id)s, %(name)s, %(evalue)s, %(bitscore)s, %(seeds)s)""", domain)
        except psycopg2.IntegrityError:
            print(cds)
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


def get_smcog_id(cur, smcog):
    """Get the smcog_id given the smCOG identifier in a feature."""
    cur.execute("SELECT smcog_id FROM antismash.smcogs WHERE name = %s", (smcog,))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("could not find matching smCOG entry for: %s" % smcog)
    return ret[0]


def handle_ripps(data):
    """Handle RiPP predictions."""
    for protocluster in data.current_region.get_unique_protoclusters():
        if not protocluster.product.endswith("peptide"):
            continue
        for motif in data.record.get_cds_motifs():
            if motif.overlaps_with(protocluster) and isinstance(motif, antismash.common.secmet.Prepeptide):
                handle_ripp(data, protocluster, motif)


def handle_ripp(data, protocluster, motif):
    """Handle a single RiPP prediction."""
    params = defaultdict(lambda: None)
    params['protocluster_id'] = data.feature_mapping[protocluster]
    params['locus_tag'] = motif.locus_tag
    params['cds_id'] = data.feature_mapping[data.record.get_cds_by_name(motif.locus_tag)]
    parse_ripp_core(motif, params)
    if params['peptide_sequence'] is None:
        print("skipping ripp without core in", protocluster)
        return
    print("inserting ripp with core:", params['peptide_sequence'])
    assert params["bridges"] is None or isinstance(params["bridges"], int), params
    data.insert("""
INSERT INTO antismash.ripps (
    protocluster_id,
    peptide_sequence,
    molecular_weight,
    monoisotopic_mass,
    alternative_weights,
    bridges,
    class,
    subclass,
    score,
    locus_tag,
    cds_id
) VALUES (
    %(protocluster_id)s,
    %(peptide_sequence)s,
    %(molecular_weight)s,
    %(monoisotopic_mass)s,
    %(alternative_weights)s,
    %(bridges)s,
    %(class)s,
    %(subclass)s,
    %(score)s,
    %(locus_tag)s,
    %(cds_id)s
)""", params)


def parse_ripp_core(feature, params):
    """Parse RIPP core features."""
    params['monoisotopic_mass'] = feature.monoisotopic_mass
    params['molecular_weight'] = feature.molecular_weight
    params['alternative_weights'] = str(feature.alternative_weights)
    params['class'] = feature.peptide_class
    params['subclass'] = feature.peptide_subclass
    params['score'] = feature.score
    params['peptide_sequence'] = feature.core
    if feature.detailed_information:
        # go for a generic approach, since the detail qualifiers vary
        bridges = feature.detailed_information.to_biopython_qualifiers().get("number_of_bridges", [None])[0]
        if bridges is not None:
            params["bridges"] = int(bridges)


def handle_asdomain(data, domain, module_id, function_id, follows):
    """Handle aSDomain features."""
    params = {}
    params['pks_signature'] = None
    params['minowa'] = None
    params['nrps_predictor'] = None
    params['stachelhaus'] = None
    params['consensus'] = None
    params['kr_activity'] = None
    params['kr_stereochemistry'] = None
    params['location'] = str(domain.location)

    params['score'] = domain.score
    params['evalue'] = domain.evalue
    params['translation'] = domain.translation
    params['locus_tag'] = domain.locus_tag
    params['cds_id'] = data.feature_mapping[data.record.get_cds_by_name(domain.locus_tag)]
    params['detection'] = domain.detection
    params['as_domain_profile_id'] = get_as_domain_profile_id(data.cursor, domain.domain_subtype or domain.domain)
    params['function_id'] = function_id
    params['module_id'] = module_id
    params['follows'] = follows

    parse_specificity(domain, params)

    try:
        as_domain_id = data.insert("""
INSERT INTO antismash.as_domains (
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
    as_domain_profile_id,
    location,
    cds_id,
    module_id,
    function_id,
    follows
) VALUES (
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
    %(as_domain_profile_id)s,
    %(location)s,
    %(cds_id)s,
    %(module_id)s,
    %(function_id)s,
    %(follows)s
) RETURNING as_domain_id""", params)
    except psycopg2.ProgrammingError:
        print("error fetching cds_id for locus_tag", params['locus_tag'])
        raise
    data.feature_mapping[domain] = as_domain_id
    print(as_domain_id, domain, "follows", follows)
    if params['consensus'] is None:
        return as_domain_id

    for substrate in params['consensus'].split('|'):
        substrate_id = get_substrate(data.cursor, substrate)
        assert substrate_id is not None
        data.insert("INSERT INTO antismash.rel_as_domains_substrates (as_domain_id, substrate_id) VALUES "
                    "(%s, %s)", (as_domain_id, substrate_id))
    return as_domain_id


def get_go_id(cursor, go_identifier):
    cursor.execute("SELECT go_id FROM antismash.gene_ontologies WHERE identifier = %s", (go_identifier,))
    ret = cursor.fetchone()
    if ret is None:
        raise ValueError("unknown GO identifier %s" % go_identifier)
    return ret[0]


def handle_pfamdomain(data, feature):
    """Handle PFAM_domain features."""
    params = {}

    params['location'] = str(feature.location)

    params['score'] = feature.score
    params['evalue'] = feature.evalue
    params['translation'] = feature.translation
    params['locus_tag'] = feature.locus_tag
    params['detection'] = feature.detection
    params['database'] = feature.database
    params['pfam_id'] = get_pfam_id(data.cursor, feature.identifier)
    params['location'] = str(feature.location)
    params['cds_id'] = data.feature_mapping[data.record.get_cds_by_name(feature.locus_tag)]

    pfam_id = data.insert("""
INSERT INTO antismash.pfam_domains (
    database,
    detection,
    score,
    evalue,
    translation,
    pfam_id,
    location,
    cds_id
) VALUES (
    %(database)s,
    %(detection)s,
    %(score)s,
    %(evalue)s,
    %(translation)s,
    %(pfam_id)s,
    %(location)s,
    %(cds_id)s
)
RETURNING pfam_domain_id""", params)

    if not feature.gene_ontologies:
        return

    statement = """
INSERT INTO antismash.pfam_go_entries (pfam_domain_id, go_id)
VALUES (%s, %s)"""
    for go_identifier in feature.gene_ontologies.ids:
        go_id = get_go_id(data.cursor, go_identifier)
        data.insert(statement, (pfam_id, go_id))


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


def get_substrate(cur, name):
    """Get the substrate_id for a substrate by name."""
    cur.execute("SELECT substrate_id FROM antismash.substrates WHERE name = %s", (name.lower(),))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("missing substrate in substrates table:", name)
    return ret[0]


def get_or_create_monomer(cur, name, substrate_id, modified):
    """Get the monomer_id for a monomer by name, create monomer entry if needed."""
    cur.execute("SELECT monomer_id FROM antismash.monomers WHERE name = %s AND substrate_id = %s", (name.lower(), substrate_id))
    ret = cur.fetchone()
    if ret:
        return ret[0]
    cur.execute("SELECT description FROM antismash.substrates WHERE substrate_id = %s", (substrate_id,))
    ret = cur.fetchone()
    desc = ("modified " if modified else "") + ret[0]
    print("inserting new monomer:", name, " -> ", desc)
    cur.execute("INSERT INTO antismash.monomers (substrate_id, name, description) VALUES (%s, %s, %s) RETURNING monomer_id", (substrate_id, name.lower(), desc))
    return cur.fetchone()[0]


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


def handle_region(data, sequence_id, region):
    """Handle cluster features."""
    assert region
    params = defaultdict(lambda: None)
    params['contig_edge'] = region.contig_edge
    params['location'] = str(region.location)
    params['sequence_id'] = sequence_id
    params["region_number"] = region.get_region_number()

    region_id = data.insert("""
INSERT INTO antismash.regions (region_number, accession, location, contig_edge)
VALUES (%(region_number)s, %(sequence_id)s, %(location)s, %(contig_edge)s)
RETURNING region_id""", params)
    params['region_id'] = region_id
    data.feature_mapping[region] = region_id
    data.current_region = region

    for product in region.products:
        product = product.lower()
        try:
            nx_create_rel_regions_types(data.cursor, params, product)
        except psycopg2.IntegrityError:
            raise RuntimeError("no definition in schema for product type: %s" % product)

    for cds in region.cds_children:
        handle_cds(data, region_id, cds)

    for candidate in region.candidate_clusters:
        handle_candidate(data, candidate)
    for protocluster in region.get_unique_protoclusters():
        handle_protocluster(data, protocluster)
        if protocluster.product == "T2PKS":
            handle_t2pks(data, protocluster)
    link_proto_to_candidates(data)

    clusterblast_results = data.module_results.get(antismash.modules.clusterblast.__name__)

    if clusterblast_results:
        if clusterblast_results.general:
            store_clusterblast(data, clusterblast_results.general.region_results[region.get_region_number()-1], 'clusterblast')
        if clusterblast_results.knowncluster:
            store_clusterblast(data, clusterblast_results.knowncluster.region_results[region.get_region_number()-1], 'knownclusterblast')
        if clusterblast_results.subcluster:
            store_clusterblast(data, clusterblast_results.subcluster.region_results[region.get_region_number()-1], 'subclusterblast')

    handle_region_nrpspks(data)

    return


def handle_region_nrpspks(data):
    modules = []
    domain_results = data.module_results[antismash.detection.nrps_pks_domains.__name__]
    for cds in data.current_region.cds_children:
        if not cds.nrps_pks:
            continue
        cds_domain_results = domain_results.cds_results[cds]
        for module, raw_module in zip(cds.modules, cds_domain_results.modules):
            modules.append(module)
            assert not module.is_multigene_module()  # this will have to be handled when these are detected
            handle_module(data, raw_module, cds_domain_results, module)

    if not modules:
        return

    nrps_results = data.module_results.get(antismash.modules.nrps_pks.__name__)
    if not nrps_results:
        return

    for candidate in data.current_region.candidate_clusters:
        for module in modules:
            if module.is_contained_by(candidate):
                data.insert("INSERT INTO antismash.rel_candidates_modules (candidate_id, module_id) VALUES (%s, %s)",
                            (data.feature_mapping[candidate], data.feature_mapping[module]))


def handle_module(data, raw_module, domain_results, secmet_module):
    statement = """
INSERT INTO antismash.modules (location, type, complete, iterative, region_id, trans_at)
VALUES (%(location)s, %(type)s, %(complete)s, %(iterative)s, %(region_id)s, %(trans_at)s)
RETURNING module_id"""
    assert raw_module.is_trans_at() in [True, False]
    values = {
        "location": str(secmet_module.location),
        "type": secmet_module.type,
        "complete": secmet_module.is_complete(),
        "iterative": secmet_module.is_iterative(),
        "region_id": data.current_region_id,
        "trans_at": raw_module.is_trans_at(),
    }

    module_id = data.insert(statement, values)
    data.feature_mapping[secmet_module] = module_id

    if not hasattr(handle_module, "_domain_function_mapping"):
        data.cursor.execute("SELECT * FROM antismash.module_domain_functions")
        handle_module._domain_function_mapping = {func: func_id for func_id, func in data.cursor.fetchall()}

    singles = {
        "starter": raw_module._starter,
        "loader": raw_module._loader,
        "carrier_protein": raw_module._carrier_protein,
        "finalisation": raw_module._end,
    }
    modification_domains = {domain_results.domain_features[component.domain] for component in raw_module._modifications}

    previous = None
    for domain in sorted(secmet_module.domains, key=lambda x: x.protein_location.start):
        function = "other"
        for label, component in singles.items():
            if component and domain is domain_results.domain_features[component.domain]:
                function = label
                break
        else:
            if domain in modification_domains:
                function = "modification"
        previous = handle_asdomain(data, domain, module_id, handle_module._domain_function_mapping[function], follows=previous)

    if secmet_module.is_complete():
        for substrate, monomer in secmet_module.get_substrate_monomer_pairs():
            substrate_id = get_substrate(data.cursor, substrate)
            modified = substrate != monomer
            monomer_id = get_or_create_monomer(data.cursor, monomer, substrate_id, modified)
            statement = """
INSERT INTO antismash.rel_modules_monomers (module_id, substrate, monomer)
VALUES (%s, %s, %s)"""
            data.insert(statement, (module_id, substrate_id, monomer_id))


def get_product_id(cur, product):
    cur.execute("SELECT bgc_type_id FROM antismash.bgc_types WHERE term = %s", (product.lower(),))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("missing product type from products table:", product)
    return ret[0]


def handle_candidate(data, candidate):
    polymer = None
    nrps_results = data.module_results.get(antismash.modules.nrps_pks.__name__)
    if nrps_results:
        for cand_result in nrps_results.region_predictions[data.current_region.get_region_number()]:
            if cand_result.candidate_cluster_number == candidate.get_candidate_cluster_number():
                polymer = cand_result.polymer
                break

    candidate_id = data.insert("""
INSERT INTO antismash.candidates (region_id, location, smiles, polymer)
VALUES (%s, %s, %s, %s)
RETURNING candidate_id""",
                               (data.current_region_id, str(candidate.location), candidate.smiles_structure, polymer))
    data.feature_mapping[candidate] = candidate_id

    for product in candidate.products:
        data.cursor.execute("""
INSERT INTO antismash.rel_candidates_types (candidate_id, bgc_type_id)
VALUES (%s, %s)""", (candidate_id, get_product_id(data.cursor, product)))


def handle_protocluster(data, protocluster):
    product_id = get_product_id(data.cursor, protocluster.product)
    protocluster_id = data.insert("""
INSERT INTO antismash.protoclusters (region_id, location, bgc_type_id)
VALUES (%s, %s, %s)
RETURNING protocluster_id""", (data.feature_mapping[data.current_region], str(protocluster.location), product_id))
    data.feature_mapping[protocluster] = protocluster_id


def link_proto_to_candidates(data):
    for candidate in data.current_region.candidate_clusters:
        cand_id = data.feature_mapping[candidate]
        for proto in candidate.protoclusters:
            proto_id = data.feature_mapping[proto]
            data.insert("""
INSERT INTO antismash.rel_candidates_protoclusters (candidate_id, protocluster_id)
VALUES (%s, %s)""", (cand_id, proto_id))


def store_clusterblast(data, results, algorithm):
    """Store XClusterBlast results in the database."""
    data.cursor.execute("SELECT algorithm_id FROM antismash.clusterblast_algorithms WHERE name = %s", (algorithm, ))
    ret = data.cursor.fetchone()
    if ret is None:
        raise ValueError('Did not find algorithm_id for {!r}!'.format(algorithm))
    algorithm_id = ret[0]

    assert isinstance(results, antismash.modules.clusterblast.results.RegionResult), type(results)
    # limit to the number of drawn hits, normally set with DEFAULT_AS_OPTIONS.cb_nclusters
    for i, hit in enumerate(results.ranking[:len(results.svg_builder.hits)]):
        ref_cluster, _ = hit
        params = {
            "algorithm_id": algorithm_id,
            "region_id": data.current_region_id,
            "rank": i + 1,
            "acc": ref_cluster.accession + (("_" + ref_cluster.cluster_label) if algorithm != "knownclusterblast" else ""),
            "description": ref_cluster.description,
            "similarity": results.svg_builder.hits[i].similarity,
        }
        data.insert("""
INSERT INTO antismash.clusterblast_hits (rank, region_id, acc, description, similarity, algorithm_id)
VALUES (%(rank)s, %(region_id)s, %(acc)s, %(description)s, %(similarity)s, %(algorithm_id)s)
RETURNING clusterblast_hit_id
    """, params)


def nx_create_rel_regions_types(cur, params, product):
    """Create relation table to bgc_types."""
    assert params.get("region_id")
    assert product
    cur.execute("""
SELECT * FROM antismash.rel_regions_types WHERE region_id = %s AND
    bgc_type_id = (SELECT bgc_type_id FROM antismash.bgc_types WHERE term = %s)""",
                (params['region_id'], product))
    ret = cur.fetchone()
    if ret is None:
        cur.execute("""
INSERT INTO antismash.rel_regions_types (region_id, bgc_type_id)
SELECT val.region_id, f.bgc_type_id FROM ( VALUES (%s, %s) ) val (region_id, bgc_type)
LEFT JOIN antismash.bgc_types f ON val.bgc_type = f.term""",
                    (params['region_id'], product))


def handle_t2pks(data, protocluster):
    protocluster_id = data.feature_mapping[protocluster]
    t2pks_id = data.insert("INSERT INTO antismash.t2pks (protocluster_id) VALUES (%s) RETURNING t2pks_id",
                           (protocluster_id,))

    t2pks_results = data.module_results[antismash.modules.t2pks.__name__].cluster_predictions[protocluster.get_protocluster_number()]

    for starter in t2pks_results.starter_units:
        statement = """
INSERT INTO antismash.t2pks_starters (t2pks_id, name, evalue, bitscore)
VALUES (%(t2pks_id)s, %(name)s, %(evalue)s, %(bitscore)s)
RETURNING domain_id"""
        values = {
            "t2pks_id": t2pks_id,
            "name": starter.name,
            "evalue": starter.evalue,
            "bitscore": starter.score,
        }
        starter_id = data.insert(statement, values)
        for elongation_joined in t2pks_results.malonyl_elongations:
            for elongation in elongation_joined.name.split("|"):
                try:
                    weight = t2pks_results.molecular_weights["%s_%s" % (starter.name, elongation)]
                except KeyError:
                    print(t2pks_results.molecular_weights)
                    raise
                elongation_count = int(elongation)
                statement = """
    INSERT INTO antismash.t2pks_starter_elongation (domain_id, elongation, weight)
    VALUES (%s, %s, %s)
    RETURNING combo_id"""
                data.insert(statement, (starter_id, elongation_count, weight))

    for cds_name, cds_results in t2pks_results.cds_predictions.items():
        cds_id = data.feature_mapping[data.record.get_cds_by_name(cds_name)]
        assert cds_id, cds_id
        for cds_result in cds_results:
            data.cursor.execute("SELECT profile_id FROM antismash.t2pks_profiles WHERE name = %s", (cds_result.ptype,))
            profile_id = data.cursor.fetchone()[0]

            statement = """
INSERT INTO antismash.t2pks_cds_domain (t2pks_id, cds_id, profile_id, protein_type, protein_function, evalue, bitscore)
VALUES (%(t2pks_id)s, %(cds_id)s, %(profile_id)s, %(ptype)s, %(pfunc)s, %(evalue)s, %(score)s)
RETURNING domain_id"""
            values = {
                "t2pks_id": t2pks_id,
                "cds_id": cds_id,
                "profile_id": profile_id,
                "ptype": cds_result.ptype,
                "pfunc": cds_result.pfunc,
                "evalue": cds_result.evalue,
                "score": cds_result.bitscore,
            }
            data.insert(statement, values)

    for product in sorted(t2pks_results.product_classes):
        statement = """
INSERT INTO antismash.t2pks_product_classes (t2pks_id, product_class)
VALUES (%s, %s)"""
        data.insert(statement, (t2pks_id, product))


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


def test_delete():
    connection = psycopg2.connect(DB_CONNECTION)
    with connection.cursor() as cursor:
        cursor.execute("SELECT * FROM antismash.genomes")
        existing = cursor.fetchall()
        if existing:
            print("existing genomes")
            for row in existing:
                print(row)
            print("deleting genomes")
            cursor.execute("DELETE FROM antismash.genomes RETURNING *")
            rows = cursor.fetchall()
            if rows:
                print("rows deleted")
                print(rows)
            connection.commit()
        else:
            print("no existing genomes")

        print("tables with entries remaining:")
        for table in sorted([
                "filenames",
                "t2pks_starter_elongation",
                "genomes",
                "t2pks_starters",
                "t2pks",
                "t2pks_product_classes",
                "dna_sequences",
                "ripps",
                "t2pks_cds_domain",
                "monomers",
                "tta_codons",
                "rel_modules_monomers",
                "regions",
                "module_modifications",
                "rel_as_domains_substrates",
                "protoclusters",
                "rel_regions_types",
                "candidates",
                "rel_candidates_protoclusters",
                "rel_candidates_types",
                "modules",
                "rel_candidates_modules",
                "clusterblast_hits",
                "profile_hits",
                "as_domains",
                "smcog_hits",
                "cdss",
                "genes",
                "pfam_domains",
                "pfam_go_entries",
            ]):
            cursor.execute("SELECT COUNT(*) FROM antismash.%s" % table)
            count = cursor.fetchone()[0]
            if count > 0:
                print("", table, count)

        connection.commit()
    connection.close()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--db', default=DB_CONNECTION, help="DB connection string to use (default: %(default)s)")
    parser.add_argument('filenames', nargs="*")
    args = parser.parse_args()
    total_duration = 0
    total_imports = 1
    for filename in args.filenames:
        print("importing", filename)
        start_time = time.time()
        try:
            main(filename, args.db)
        except Exception as err:
            print(filename, err)
        finally:
            end_time = time.time()
            import_duration = end_time - start_time
            print("took", round(import_duration, 2), "seconds", end="\t")
            total_imports += 1
            total_duration += import_duration
            print("average:", round(total_duration/total_imports, 2))
