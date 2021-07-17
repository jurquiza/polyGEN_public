import unittest
import polygen
import json

from engine import scarless_gg,Part,polyToJson,flattn

def equalizeResults(s):
    s = s.replace('\'', '\"')
    s = s.replace('None', 'null')
    return s
    

class TestWebApplication(unittest.TestCase):
    
    def setUp(self):
        polygen.app.testing = True
        self.app = polygen.app.test_client()
    
    def postPTG(self, poltype_input, enzm_input, PTG_name, oligo_prefix, oligo_index, sequence_spacers, bb_ovrhng, add_ovrhng, min_temp, max_temp, borderPrimers=False):
        if borderPrimers:
            return self.app.post('/ptg', data=dict(submit_button='submit', poltype_input=poltype_input, enzm_input=enzm_input, PTG_name=PTG_name, oligo_prefix=oligo_prefix, oligo_index=oligo_index, sequence_spacers=sequence_spacers, bb_ovrhng=bb_ovrhng, add_ovrhng=add_ovrhng, min_temp=min_temp, max_temp=max_temp, borderPrimers=borderPrimers), follow_redirects=True)
        else: 
            return self.app.post('/ptg', data=dict(submit_button='submit', poltype_input=poltype_input, enzm_input=enzm_input, PTG_name=PTG_name, oligo_prefix=oligo_prefix, oligo_index=oligo_index, sequence_spacers=sequence_spacers, bb_ovrhng=bb_ovrhng, add_ovrhng=add_ovrhng, min_temp=min_temp, max_temp=max_temp), follow_redirects=True)
    
    def test_scarlessgg_test0_gRNAs_flask(self):
        self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test0_gRNAs', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='gRNA;aagttttaaatcaatctaaa|gRNA;gcatcagcaccttgtcgcct|gRNA;ggatgatttctggaattcgc', 
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test0_gRNAs_2021-03-26/test0_gRNAs_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data
    
    def test_scarlessgg_test1_PE_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test1_PE', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='pegRNA;CTCGTGACCACCCTGACCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTGCACGCAGTGGGTCAGGGTGGTCA|gRNA;AGAAGTCGTGCTGCTTCATG', 
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test1_PE_2021-03-26/test1_PE_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test2_smRNAs_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test2_smRNAs', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;aacattcaacgctgtcggtgagt|smRNA;ctccttcacccgggcggtacc|smRNA;taagtgcttacctgtttgggcat', 
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test2_smRNAs_2021-03-26/test2_smRNAs_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data
            
    def test_scarlessgg_test3_mixed_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test3_mixed', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG', 
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test3_mixed_2021-03-26/test3_mixed_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test4_ca_flask(self):
        rv = self.postPTG(poltype_input='ca', 
                          enzm_input='bsai', 
                          PTG_name='test4_ca', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='crRNA;gatggtgcttcaaatgagat|crRNA;aatggttctcttcttgatga|crRNA;gaatggttctcttcttgatg', 
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test4_ca_2021-03-26/test4_ca_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test5_bsmbi_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsmbi', 
                          PTG_name='test5_bsmbi', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test5_bsmbi_2021-03-26/test5_bsmbi_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test6_btgzi_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='btgzi', 
                          PTG_name='test6_btgzi', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test6_btgzi_2021-03-26/test6_btgzi_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test7_bpii_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bpii', 
                          PTG_name='test7_bpii', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test7_bpii_2021-03-31/test7_bpii_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test8_addedLinkers_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test8_addedLinkers', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='atct', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test8_addedLinkers_2021-03-26/test8_addedLinkers_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test9_borderLinkers_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test9_borderLinkers', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='gaag;ttgt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test9_borderLinkers_2021-03-26/test9_borderLinkers_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test10_lowerTemp_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test10_lowerTemp', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='45', 
                          max_temp='54')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test10_lowerTemp_2021-03-26/test10_lowerTemp_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test11_higherTemp_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test11_higherTemp', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='64', 
                          max_temp='75')
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test11_higherTemp_2021-03-26/test11_higherTemp_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

    def test_scarlessgg_test13_staticBorderPrimers_flask(self):
        rv = self.postPTG(poltype_input='ptg', 
                          enzm_input='bsai', 
                          PTG_name='test13_staticBorderPrimers', 
                          oligo_prefix='oUU', 
                          oligo_index='0', 
                          sequence_spacers='smRNA;gctaggtgcaacaagttcaat|pegRNA;CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA|gRNA;CGATGCCCTTCAGCTCGATG',
                          bb_ovrhng='tgcc;gttt', 
                          add_ovrhng='', 
                          min_temp='55', 
                          max_temp='65',
                          borderPrimers=True)
        dt = self.app.post('/primer_list', data=dict(), follow_redirects=True)
        with open('testing/pg_test13_staticBorderPrimers_2021-05-07/test13_staticBorderPrimers_raw.json', 'r') as f:
            expected = bytes(equalizeResults(str(json.load(f))), 'utf-8')
            assert expected in dt.data

