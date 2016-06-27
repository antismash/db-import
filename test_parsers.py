from import_genbank import (
    parse_smcog,
    parse_specificity,
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
