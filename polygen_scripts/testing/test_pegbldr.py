import unittest

from engine import pegbldr

seq = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacccacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa'


class TestPegbldr(unittest.TestCase):
    def test_peg_mut_single_forw(self):
        '''
        test if one pegRNA with a single mutation is designed with forward spacer
        '''
        edts = [['201', 'T', 'mut']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CTCGTGACCACCCTGACCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTGCACGCAGTGGGTCAGGGTGGTCA', 'f'], ['gRNA0', 'gRNA', 'AGAAGTCGTGCTGCTTCATG', 'f']], None)
        self.assertEqual(result, expected)
    
    def test_peg_mut_single_rev(self):
        '''
        test if one pegRNA with a single mutation is designed with reverse spacer
        '''
        edts = [['635', 'A', 'mut']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CATGTGATCGCGCTTCTCGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAAAGACCCAAACGAGAAGCGCGATCA', 'r'], ['gRNA0', 'gRNA', 'CAGAACACCCCCATCGGCGA', 'r']], None)
        self.assertEqual(result, expected)
    
    def test_peg_mut_multiple(self):
        '''
        test if one pegRNA with multiple mutations is designed reliably
        '''
        edts = [['280,283', 'C,T', 'mut']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'TGAAGAAGATGGTGCGCTCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCAGGCTACGCCCTGGAGCGCACCATCTTC', 'r'], ['gRNA0', 'gRNA', 'TTCAAGTCCGCCATGCCCGA', 'r']], None)
        self.assertEqual(result, expected)
        
    def test_peg_del_single_forw(self):
        '''
        test if one pegRNA with a single deletion is designed with forward spacer
        '''
        edts = [['342', '344', 'del']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA', 'f'], ['gRNA0', 'gRNA', 'CGATGCCCTTCAGCTCGATG', 'f']], None)
        self.assertEqual(result, expected)
        
    def test_peg_del_single_rev(self):
        '''
        test if one pegRNA with a single deletion is designed with reverse spacer
        '''
        edts = [['504', '508', 'del']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CGCTGCCGTCCTCGATGTTGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTTCAAGATCCAACATCGAGGACGGC', 'r'], ['gRNA0', 'gRNA', 'CAGCCACAACGTCTATATCA', 'r']], None)
        self.assertEqual(result, expected)
    
    def test_peg_del_multiple(self):
        '''
        test if one pegRNA with multiple deletions is designed reliably
        '''
        edts = [['425,430', '427,435', 'del']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CAACATCCTGGGGCACAAGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACGTTGTGGCTGTTGTACAGCTTGTGCCCCAGGAT', 'f'], ['gRNA0', 'gRNA', 'GATGCCGTTCTTCTGCTTGT', 'f']], None)
        self.assertEqual(result, expected)
    
    def test_peg_ins_single_forw(self):
        '''
        test if one pegRNA with a single insertion is designed with forward spacer
        '''
        edts = [['50', 'ACGT', 'ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CACCGGGGTGGTGCCCATCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAGCTCGACGTACCAGGATGGGCACCACCCC', 'f'], ['gRNA0', 'gRNA', 'CGCCGGACACGCTGAACTTG', 'f']], None)
        self.assertEqual(result, expected)
        
    def test_peg_ins_single_rev(self):
        '''
        test if one pegRNA with a single insertion is designed with reverse spacer
        '''
        edts = [['603','CCC','ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'TCAGCTTGGACTGGGTGCTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCAACCACTACCCCCTGAGCACCCAGTCCAAG', 'r'], ['gRNA0', 'gRNA', 'TACCAGCAGAACACCCCCAT', 'r']], None)
        self.assertEqual(result, expected)
    
    def test_peg_ins_multiple(self):
        '''
        test if one pegRNA with multiple insertions is designed reliably
        '''
        edts = [['120,130', 'AAAA,GGGG', 'ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'GGTGGTGCAGATGAACTTCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCACCTACAAAAGGCAAGCTGAGGGGCCCTGAAGTTCATCTGCAC', 'r'], ['gRNA0', 'gRNA', 'AAGTTCAGCGTGTCCGGCGA', 'r']], None)
        self.assertEqual(result, expected)
    
    def test_peg_multiple(self):
        '''
        test if multiple independent mutations are designed reliably on forward and reverse strand
        '''
        edts = [['201', 'T', 'mut'], ['504', '508', 'del'], ['120,130', 'AAAA,GGGG', 'ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CTCGTGACCACCCTGACCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTGCACGCAGTGGGTCAGGGTGGTCA', 'f'], ['gRNA0', 'gRNA', 'AGAAGTCGTGCTGCTTCATG', 'f'], ['pegRNA1', 'pegRNA', 'CGCTGCCGTCCTCGATGTTGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTTCAAGATCCAACATCGAGGACGGC', 'r'], ['gRNA1', 'gRNA', 'CAGCCACAACGTCTATATCA', 'r'], ['pegRNA2', 'pegRNA', 'GGTGGTGCAGATGAACTTCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCACCTACAAAAGGCAAGCTGAGGGGCCCTGAAGTTCATCTGCAC', 'r'], ['gRNA2', 'gRNA', 'AAGTTCAGCGTGTCCGGCGA', 'r']], None)
        self.assertEqual(result, expected)
