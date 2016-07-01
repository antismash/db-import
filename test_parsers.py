from collections import defaultdict
from import_genbank import (
    parse_smcog,
    parse_specificity,
    parse_ripp_core,
    parse_monomers,
    parse_clusterblast_line,
)


class FakeFeature(object):
    def __init__(self):
        self.qualifiers = {}


def test_parse_smcog():
    '''Test parse_smcog works as expected'''
    test_cases = [
        ("smCOG: SMCOG1094:ferredoxin (Score: 187.7; E-value: 2.562e-56);", ('SMCOG1094', '187.7', '2.562e-56')),
        ("smCOG: SMCOG1012:4'-phosphopantetheinyl_transferase (Score: 106.6; E-value: 1.9e-32);", ('SMCOG1012', '106.6', '1.9e-32')),
        ("smCOG: SMCOG1000:ABC_transporter_ATP-binding_protein (Score: 145.7; E-value: 2.9e-44);", ('SMCOG1000', '145.7', '2.9e-44')),
        ("smCOG: SMCOG1001:short-chain_dehydrogenase/reductase_SDR (Score: 220.2; E-value: 5.1e-67);", ('SMCOG1001', '220.2', '5.1e-67')),
        ("smCOG: SMCOG1171:transcriptional_regulator,_MerR_family (Score: 67.3; E-value: 3.6e-20);", ('SMCOG1171', '67.3', '3.6e-20')),
        ("smCOG: SMCOG1064:glucose-1-phosphate_adenylyl/thymidylyltransferas e (Score: 78.0; E-value: 1.1e-23);", ('SMCOG1064', '78.0', '1.1e-23')),
        ("smCOG: SMCOG1177:asparagine_synthase_(glutamine-hydrolyzing) (Score: 43.2; E-value: 1.392e-12);", ('SMCOG1177', '43.2', '1.392e-12')),
        ("smCOG: SMCOG1212:sodium:dicarboxylate_symporter (Score: 325.8; E-value: 5.8e-99);", ('SMCOG1212', '325.8', '5.8e-99')),
        ("smCOG: SMCOG1283:2`,3`-cyclic-nucleotide_2`-phosphodiesterase (Score: 162.6; E-value: 2.5e-49);", ('SMCOG1283', '162.6', '2.5e-49')),
    ]

    for test_string, expected in test_cases:
        fake = FakeFeature()
        fake.qualifiers['note'] = [test_string]
        assert parse_smcog(fake) == expected


def test_parse_specificity():
    '''Test parse_specificity'''
    params = {}
    specificities = [
        "PKS signature: mal",
        "Minowa: mal",
        "consensus: ccmal",
        "KR activity: active",
        "KR stereochemistry: ?",
        "NRPSpredictor2 SVM: orn,lys,arg",
        "Stachelhaus code: orn",
    ]


    expected = {
        'pks_signature': 'mal',
        'minowa': 'mal',
        'consensus': 'ccmal',
        'kr_activity': True,
        'kr_stereochemistry': '?',
        'nrps_predictor': 'orn,lys,arg',
        'stachelhaus': 'orn'
    }

    fake = FakeFeature()
    fake.qualifiers['specificity'] = specificities
    parse_specificity(fake, params)
    assert params == expected


def test_parse_ripp_core():
    '''Test parse_ripp_core'''
    params = defaultdict(lambda: None)
    notes = [
        'Totally unrelated: nonsense',
        'monoisotopic mass: 3333.6',
        'molecular weight: 3336.0',
        'alternative weights: 3354.0; 3372.1; 3390.1; 3408.1',
        'number of bridges: 5',
        'predicted core seq: ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK',
        'predicted class: Class-I',
        'score: 26.70',
    ]

    expected = {
        'peptide_sequence': 'ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK',
        'molecular_weight': '3336.0',
        'monoisotopic_mass': '3333.6',
        'alternative_weights': '3354.0; 3372.1; 3390.1; 3408.1',
        'bridges': '5',
        'class': 'Class-I',
        'score': '26.70',
    }

    fake = FakeFeature()
    fake.qualifiers['note'] = notes

    parse_ripp_core(fake, params)
    assert params == expected


def test_parse_monomers():
    '''Test parse_monomers'''
    tests = [
        ('Monomers prediction: (ccmal) + (pk)', (['ccmal', 'pk'], '(ccmal) + (pk)')),
        ('Monomers prediction: (orn-thr-orn)', (['orn', 'thr', 'orn'], '(orn-thr-orn)')),
        ('Monomers prediction: (ser-thr-trp-asp-asp-hpg) + (asp-gly-asn) + (nrp-trp)',
            (['ser', 'thr', 'trp', 'asp', 'asp', 'hpg', 'asp', 'gly', 'asn', 'nrp', 'trp'],
             '(ser-thr-trp-asp-asp-hpg) + (asp-gly-asn) + (nrp-trp)')),
        ('Monomers prediction: (cys) + (pk-mal)', (['cys', 'pk', 'mal'], '(cys) + (pk-mal)')),
        ('Monomers prediction: (mal-ccmal-ccmal) + (ccmal-ccmal) + (ccmal)',
            (['mal', 'ccmal', 'ccmal', 'ccmal', 'ccmal', 'ccmal'], '(mal-ccmal-ccmal) + (ccmal-ccmal) + (ccmal)')),
        ('Monomers prediction: (nrp) + (cys)', (['nrp', 'cys'], '(nrp) + (cys)')),
        ('Monomers prediction: (redmmal)', (['redmmal'], '(redmmal)')),
        ('Monomers prediction: (dhb) + (nrp-cys) + (cys)', (['dhb', 'nrp', 'cys', 'cys'], '(dhb) + (nrp-cys) + (cys)')),
        ('Total nonsense: here', ([], None)),
    ]

    for test_string, expected in tests:
        fake = FakeFeature()
        fake.qualifiers['note'] = [test_string]
        assert parse_monomers(fake) == expected


def test_parse_clusterblast_line():
    '''Test parse_clusterblast_line'''
    tests = [
        ('1. CP002365_c2       Lactococcus lactis subsp. lactis CV56, complete genome. (60% of genes show similarity)',
            {'rank': 1, 'acc': 'CP002365_c2', 'description': 'Lactococcus lactis subsp. lactis CV56, complete genome.',
             'similarity': 60}),
        ('10. AEXT01000007_c1  Streptococcus agalactiae FSL S3-026 contig07, whole genom... (26% of genes show similarity)',
            {'rank': 10, 'acc': 'AEXT01000007_c1', 'description': 'Streptococcus agalactiae FSL S3-026 contig07, whole genom...',
             'similarity': 26}),
        ('2. BGC0000536_c1     Nisin_Q_biosynthetic_gene_cluster (100% of genes show similarity)',
            {'rank': 2, 'acc': 'BGC0000536_c1', 'description': 'Nisin_Q_biosynthetic_gene_cluster', 'similarity': 100}),
        ('6. AB362350_c1    Lactococcus lactis DNA, nisin Q gene cluster (nisQ, niqB,niqT,... (26% of genes show similarity)',
            {'rank': 6, 'acc': 'AB362350_c1', 'description': 'Lactococcus lactis DNA, nisin Q gene cluster (nisQ, niqB,niqT,...',
             'similarity': 26}),
        ('1. AL939104_c2    Streptomyces coelicolor A3(2) complete genome; segment 1/29. (100% of genes show similarity)',
            {'rank': 1, 'acc': 'AL939104_c2', 'description': 'Streptomyces coelicolor A3(2) complete genome; segment 1/29.',
             'similarity': 100}),
        ('6. AJ007731_c1    Streptomyces coelicolor scbR gene, scbA gene, ORFs A,B,X & Z. (16% of genes show similarity)',
            {'rank': 6, 'acc': 'AJ007731_c1', 'description': 'Streptomyces coelicolor scbR gene, scbA gene, ORFs A,B,X & Z.',
             'similarity': 16}),
    ]

    for line, expected in tests:
        assert parse_clusterblast_line(line) == expected
