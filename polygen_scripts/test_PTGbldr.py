import unittest

from engine import PTGbldr

class TestPTGbldr(unittest.TestCase):
    def test_PTGbldr_PTG(self):
        '''
        Test if a PTG is correctly designed for all input types
        '''
        insrts = [['gRNA0', 'gRNA', 'AAGGCCTTAAGGCCTTAAGG'], ['pegRNA0', 'pegRNA', 'CACCGGGGTGGTGCCCATCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAGCTCGACGTACCAGGATGGGCACCACCCC'], ['smRNA0', 'smRNA', 'ACGTACGTACGTACGTACGTACGTACGT']]
        partsList = PTGbldr(insrts, poltype_bldr='ptg')
        result = [vars(part) for part in partsList]
        expected = [{'name': 'tRNA', 'sequence': 'aacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'tRNA'}, {'name': 'gRNA0', 'sequence': 'AAGGCCTTAAGGCCTTAAGGgttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcaacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'gRNA'}, {'name': 'pegRNA0', 'sequence': 'CACCGGGGTGGTGCCCATCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAGCTCGACGTACCAGGATGGGCACCACCCC', 'type': 'pegRNA'}, {'name': 'tRNA', 'sequence': 'aacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'tRNA'}, {'name': 'smRNA0', 'sequence': 'ACGTACGTACGTACGTACGTACGTACGTaacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca', 'type': 'smRNA'}]
        self.assertEqual(result, expected)
    
    def test_PTGbldr_CA(self):
        '''
        Test if a CA is correctly designed
        '''
        insrts = [['crRNA0', 'crRNA', 'ACGTACGTACGTACGTACGT'], ['crRNA1', 'crRNA', 'AAAACCCCGGGGTTTTAAAA'], ['crRNA2', 'crRNA', 'AACCGGTTAACCGGTTAACC']]
        partsList = PTGbldr(insrts, poltype_bldr='ca')
        result = [vars(part) for part in partsList]
        expected = [{'name': 'crRNA0', 'sequence': 'aatttctactgttgtagatACGTACGTACGTACGTACGT', 'type': 'crRNA'}, {'name': 'crRNA1', 'sequence': 'aatttctactgttgtagatAAAACCCCGGGGTTTTAAAA', 'type': 'crRNA'}, {'name': 'crRNA2', 'sequence': 'aatttctactgttgtagatAACCGGTTAACCGGTTAACC', 'type': 'crRNA'}, {'name': 'DR', 'sequence': 'aatttctactgttgtagat', 'type': 'DR'}]
        self.assertEqual(result, expected)
