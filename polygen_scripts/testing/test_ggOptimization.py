import unittest

from engine import golden_gate_optimization, Part

tRNA = 'aacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca'
scaffld = 'gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc'
DR = 'aatttctactgttgtagat'

tRNA0 = Part('tRNA0','tRNA',tRNA)
pegRNA_BFP = Part('pegRNA_BFP', 'pegRNA', 'CTCGTGACCACCCTGACCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTGCACGCAGTGGGTCAGGGTGGTCA')
gRNA_BFP = Part('gRNA_BFP', 'gRNA', 'AGAAGTCGTGCTGCTTCATGgttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcaacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca')
smRNA0 = Part('smRNA0', 'smRNA', 'AGGCTTAGCTAGGCCTATTACGCTAaacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca')

DR0 = Part('DR0', 'DR', DR)
crRNA0 = Part('crRNA0', 'crRNA', 'aatttctactgttgtagatAAGTAGGCTTACGGTATCGT')
crRNA1 = Part('crRNA1', 'crRNA', 'aatttctactgttgtagatCGGCTAGGAAATTAGCCGAA')
crRNA2 = Part('crRNA2', 'crRNA', 'aatttctactgttgtagatGGCTTAAATCGCCTACCTTT')

oh_set0 = [['tgga', 'tcat', 'tgat', 'atgt', 'tctt', 'tgtc', 'cctt', 'cttg', 'tcgt', 'tagc', 'gccg', 'ggtc', 'tgag']]
oh_set1 = [['gggg', 'gggc', 'tctc', 'ggat', 'tacc', 'cctt', 'tgta', 'tctt', 'gatg', 'atct', 'tcgt', 'tccc', 'ctgt', 'tcct', 'gagt', 'tcat', 'tgcg', 'tgtt', 'gttg', 'cttg', 'ttca', 'cgtt', 'gctt', 'tagc', 'cctg', 'ccct', 'ggcg', 'ttgc', 'gtgt', 'tgga', 'tttg']]


class TestGoldenGateOptimization(unittest.TestCase):
    def test_GGO_PTG(self):
        '''
        Test if GGO functions for PTGs
        '''
        pList = [tRNA0, pegRNA_BFP, tRNA0, gRNA_BFP, smRNA0]
        result = golden_gate_optimization(pList, oh_set0, 'ptg')
        expected = ('ctcgt', 'actgcacgcagtgggtc', 'agaagtcgtgctgcttcat', 'aggcttagc')
        self.assertEqual(result, expected)
    
    def test_GGO_CA(self):
        '''
        Test if GGO functions for CAs
        '''
        pList = [crRNA0, crRNA1, crRNA2, DR0]
        result = golden_gate_optimization(pList, oh_set1, 'ca')
        expected = ('aagtaggctt', 'cggctaggaaattagc', 'ggcttaaatcgcctacc')
        self.assertEqual(result, expected)
