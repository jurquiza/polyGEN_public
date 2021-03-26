import unittest

from engine import PTGbldr

class TestPTGbldr(unittest.TestCase):
    def test_PTGbldr_PTG(self):
        '''
        Test if a PTG is correctly designed for all input types
        '''
        insrts = [['gRNA', 'AAGGCCTTAAGGCCTTAAGG'], ['pegRNA', 'CACCGGGGTGGTGCCCATCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAGCTCGACGTACCAGGATGGGCACCACCCC'], ['smRNA', 'ACGTACGTACGTACGTACGTACGTACGT']]
        partsList = PTGbldr('test0', insrts, poltype='ptg')
        result = [vars(part) for part in partsList]
        expected = [{'name': 'test0_0', 'sequence': 'aacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'tRNA'}, {'name': 'test0_1', 'sequence': 'AAGGCCTTAAGGCCTTAAGGgttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcaacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'gRNA'}, {'name': 'test0_2', 'sequence': 'CACCGGGGTGGTGCCCATCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAGCTCGACGTACCAGGATGGGCACCACCCC', 'type': 'pegRNA'}, {'name': 'test0_3', 'sequence': 'aacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'tRNA'}, {'name': 'test0_4', 'sequence': 'ACGTACGTACGTACGTACGTACGTACGTaacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'smRNA'}]
        self.assertEqual(result, expected)
    
    def test_PTGbldr_CA(self):
        '''
        Test if a CA is correctly designed
        '''
        insrts = [['crRNA', 'ACGTACGTACGTACGTACGT'], ['crRNA', 'AAAACCCCGGGGTTTTAAAA'], ['crRNA', 'AACCGGTTAACCGGTTAACC']]
        partsList = PTGbldr('test1', insrts, poltype='ca')
        result = [vars(part) for part in partsList]
        expected = [{'name': 'test1_0', 'sequence': 'aatttctactgttgtagatACGTACGTACGTACGTACGT', 'type': 'crRNA'}, {'name': 'test1_1', 'sequence': 'aatttctactgttgtagatAAAACCCCGGGGTTTTAAAA', 'type': 'crRNA'}, {'name': 'test1_2', 'sequence': 'aatttctactgttgtagatAACCGGTTAACCGGTTAACC', 'type': 'crRNA'}, {'name': 'test1_3', 'sequence': 'aatttctactgttgtagat', 'type': 'DR'}]
        self.assertEqual(result, expected)
