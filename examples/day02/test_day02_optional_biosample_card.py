import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from day02_optional_biosample_card import validate_record, run

BASE=json.loads(Path(__file__).with_name('day02_synthetic_biosample_card.json').read_text(encoding='utf-8'))

class Day02Test(unittest.TestCase):
    def test_default_packet(self):
        x=validate_record(deepcopy(BASE))
        self.assertEqual(x['mode'],'OFFLINE_ONLY')
        self.assertEqual(x['null_count'],1)
        self.assertEqual(x['clinical_causality'],'NOT_ESTABLISHED')
    def test_missing_remains_null(self):
        self.assertIsNone(validate_record(deepcopy(BASE))['rows'][2]['score'])
    def test_zero_not_missing(self):
        self.assertEqual(validate_record(deepcopy(BASE))['rows'][1]['interpretation'],'SYNTHETIC_ZERO')
    def test_nonzero_separate(self):
        self.assertEqual(validate_record(deepcopy(BASE))['rows'][0]['interpretation'],'SYNTHETIC_NONZERO')
    def test_rejects_other_source(self):
        x=deepcopy(BASE);x['data_kind']='PATIENT';self.assertRaises(ValueError,validate_record,x)
    def test_rejects_other_genome(self):
        x=deepcopy(BASE);x['genome_build']='GRCh38';self.assertRaises(ValueError,validate_record,x)
    def test_rejects_bad_assay(self):
        x=deepcopy(BASE);x['assay']='ATAC';self.assertRaises(ValueError,validate_record,x)
    def test_rejects_null_mislabeled_zero(self):
        x=deepcopy(BASE);x['rows'][2]['score_status']='SYNTHETIC_NUMBER';self.assertRaises(ValueError,validate_record,x)
    def test_rejects_nan(self):
        x=deepcopy(BASE);x['rows'][0]['score']=float('nan');self.assertRaises(ValueError,validate_record,x)
    def test_rejects_duplicate(self):
        x=deepcopy(BASE);x['rows'][2]['biosample']='Blood-like';self.assertRaises(ValueError,validate_record,x)
    def test_rejects_model_called(self):
        x=deepcopy(BASE);x['model_prediction_performed']=True;self.assertRaises(ValueError,validate_record,x)
    def test_rejects_fake_experiment(self):
        x=deepcopy(BASE);x['rows'][1]['experimental_observation']='PERFORMED';self.assertRaises(ValueError,validate_record,x)
    def test_written_json_csv(self):
        with TemporaryDirectory() as d:
            f=Path(d)/'out'
            o=run(Path(__file__).with_name('day02_synthetic_biosample_card.json'),f)
            self.assertEqual(o['real_api_calls_this_run'],0)
            self.assertEqual(json.loads((f/'day02_biosample_card.json').read_text())['null_count'],1)
            self.assertTrue((f/'day02_biosample_card.csv').exists())

if __name__=='__main__':
    unittest.main()
