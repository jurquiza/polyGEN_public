import unittest
import json

from engine import scarless_gg,Part,polyToJson,flattn

class TestScarlessGG(unittest.TestCase):
    def test_scarlessgg_test0_gRNAs(self):
        '''
        Test if a standard PTG ([55-65],30,[tgcc,gttt],[],ptg,bsai) from 3 gRNAs is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test0_gRNAs']]
                
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test0_gRNAs_2021-08-17/test0_gRNAs_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])
    
    
    def test_scarlessgg_test1_PE(self):
        '''
        Test if a standard PTG ([55-65],30,[tgcc,gttt],[],ptg,bsai) for a PE (pegRNA+gRNA) is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test1_PE']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test1_PE_2021-08-17/test1_PE_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])
          
    
    def test_scarlessgg_test2_smRNAs(self):
        '''
        Test if a standard PTG ([55-65],30,[tgcc,gttt],[],ptg,bsai) for smRNAs is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test2_smRNAs']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test2_smRNAs_2021-08-17/test2_smRNAs_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])


    def test_scarlessgg_test3_mixed(self):
        '''
        Test if a standard PTG ([55-65],30,[tgcc,gttt],[],ptg,bsai) for mixed gRNAs, pegRNAs and smRNAs is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test3_mixed']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test3_mixed_2021-08-17/test3_mixed_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])


    def test_scarlessgg_test4_ca(self):
        '''
        Test if a CA ([55-65],30,[tgcc,gttt],[],ca,bsai) for crRNAs is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test4_ca']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ca',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test4_ca_2021-08-17/test4_ca_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])
            

    def test_scarlessgg_test5_bsmbi(self):
        '''
        Test if a bsmbi mixed PTG ([55-65],30,[tgcc,gttt],[],ptg,bsmbi) is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test5_bsmbi']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsmbi'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test5_bsmbi_2021-08-17/test5_bsmbi_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])


    def test_scarlessgg_test6_btgzi(self):
        '''
        Test if a btgzi mixed PTG ([55-65],30,[tgcc,gttt],[],ptg,btgzi) is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test6_btgzi']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'btgzi'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test6_btgzi_2021-08-17/test6_btgzi_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])


    def test_scarlessgg_test7_bbsi(self):
        '''
        Test if a bbsi mixed PTG ([55-65],30,[tgcc,gttt],[],ptg,bbsi) is constructed correctly
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test7_bbsi']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bbsi'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test7_bbsi_2021-08-17/test7_bbsi_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])
            

    def test_scarlessgg_test8_addedLinkers(self):
        '''
        Test if a mixed PTG ([55-65],30,[tgcc,gttt],[atct],ptg,bsai) can work with additional linkers
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test8_addedLinkers']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc', 'gttt'],
                'additional_overhangs': ['atct'],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test8_addedLinkers_2021-08-17/test8_addedLinkers_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])


    def test_scarlessgg_test9_borderLinkers(self):
        '''
        Test if a mixed PTG ([55-65],30,[gaag,ttgt],[],ptg,bsai) can work with different border linkers
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test9_borderLinkers']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 30,
                'bb_overlaps': ['gaag','ttgt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test9_borderLinkers_2021-08-17/test9_borderLinkers_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])


    def test_scarlessgg_test10_lowerTemp(self):
        '''
        Test if the primers of a mixed PTG ([45-54],30,[tgcc,gttt],[],ptg,bsai) respond to a lower temperature by becoming shorter
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test10_lowerTemp']]
        
        args = {'parts_list': inpt, 
                'tm_range': [45, 54],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc','gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test10_lowerTemp_2021-08-17/test10_lowerTemp_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])
            

    def test_scarlessgg_test11_higherTemp(self):
        '''
        Test if the primers of a mixed PTG ([64-75],30,[tgcc,gttt],[],ptg,bsai) respond to a higher temperature by becoming longer
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test11_higherTemp']]
        
        args = {'parts_list': inpt, 
                'tm_range': [64, 75],
                'max_ann_len': 30,
                'bb_overlaps': ['tgcc','gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        with open('testing/pg_test11_higherTemp_2021-08-17/test11_higherTemp_raw.json', 'r') as f:
            expected = json.load(f)
            expected['polycistron']['oligos'] = []
            expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected['polycistron'])


    def test_scarlessgg_test12_annLen(self):
        '''
        Test if the primers of a mixed PTG ([55-65],25,[tgcc,gttt],[],ptg,bsai) respond to a different maximal annealing length by producing a different result
        '''
        
        inpt = []
        with open('testing/test_input_data.json', 'r') as f:
            inpt = [Part.from_json(prt) for prt in json.load(f)['test12_annLen']]
        
        args = {'parts_list': inpt, 
                'tm_range': [55, 65],
                'max_ann_len': 25,
                'bb_overlaps': ['tgcc','gttt'],
                'additional_overhangs': [],
                'poltype': 'ptg',
                'enzm': 'bsai'}
        result = {}
        result = polyToJson(scarless_gg(**args))
        
        for f in result['features']:
            f['location']['_start'] = int(f['location']['_start'])
            f['location']['_end'] = int(f['location']['_end'])
            f['qualifiers'] = {}
        
        for p in result['parts']:
            p.pop('localisation')
        result.pop('oligos')
        
        with open('testing/pg_test12_annLen_2021-03-26/test12_annLen_raw.json', 'r') as f:
            expected = json.load(f)
            #expected['polycistron']['oligos'] = []
            #expected['polycistron']['features'] = expected['polycistron']['features'][:-len(expected['polycistron']['parts']*2)]
            
            self.assertEqual(result, expected)
