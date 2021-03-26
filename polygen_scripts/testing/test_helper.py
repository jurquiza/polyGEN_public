import unittest

from engine import reverse_complement, reverse, complement, flattn, Diff


class TestReverseComplement(unittest.TestCase):
    def test_lowercase(self):
        '''
        Test that it can rc a lowercase string
        '''
        data = 'atcg'
        result = reverse_complement(data)
        expected = 'cgat'
        self.assertEqual(result, expected)
        
    def test_uppercase(self):
        '''
        Test that it can rc an uppercase string
        '''
        data = 'GCTA'
        result = reverse_complement(data)
        expected = 'TAGC'
        self.assertEqual(result, expected)


class TestReverse(unittest.TestCase):
    def test_lowercase(self):
        '''
        Test that it can reverse a lowercase string
        '''
        data = 'atcg'
        result = reverse(data)
        expected = 'gcta'
        self.assertEqual(result, expected)
        
    def test_uppercase(self):
        '''
        Test that it can reverse an uppercase string
        '''
        data = 'GCTA'
        result = reverse(data)
        expected = 'ATCG'
        self.assertEqual(result, expected)
        

class TestComplement(unittest.TestCase):
    def test_lowercase(self):
        '''
        Test that it can complement a lowercase string
        '''
        data = 'atcg'
        result = complement(data)
        expected = 'tagc'
        self.assertEqual(result, expected)
        
    def test_uppercase(self):
        '''
        Test that it can complement an uppercase string
        '''
        data = 'GCTA'
        result = complement(data)
        expected = 'CGAT'
        self.assertEqual(result, expected)


class TestFlattn(unittest.TestCase):
    def test_oneNested(self):
        '''
        Test that it can flatten a list with one nested list inside
        '''
        data = [[1,2,3]]
        result = flattn(data)
        expected = [1,2,3]
        self.assertEqual(result, expected)
        
    def test_multipleNested(self):
        '''
        Test that it can flatten a list with multiple nested lists inside
        '''
        data = [[1,2,3],[4,5,6]]
        result = flattn(data)
        expected = [1,2,3,4,5,6]
        self.assertEqual(result, expected)
        

class TestDiff(unittest.TestCase):
    def test_firstEmpty(self):
        '''
        Test that it will simply return all values of the second list
        '''
        data1 = []
        data2 = [1,2,3]
        result = Diff(data1, data2)
        expected = [1,2,3]
        self.assertEqual(result, expected)
    
    def test_secondEmpty(self):
        '''
        Test that it will simply return all values of the second list
        '''
        data1 = [1,2,3]
        data2 = []
        result = Diff(data1, data2)
        expected = [1,2,3]
        self.assertEqual(result, expected)
    
    def test_normal(self):
        '''
        Test that it will return all unique values
        '''
        data1 = [1,2,3]
        data2 = [1,3,5]
        result = Diff(data1, data2)
        expected = [2,5]
        self.assertEqual(result, expected)
    
    def test_types(self):
        '''
        Test that it can work with integers, floats, strings
        '''
        data1 = [1,2,1.1,2.1,'a','b']
        data2 = [1,3,1.1,2.2,'a','c']
        result = Diff(data1, data2)
        expected = [2,2.1,'b',3,2.2,'c']
        self.assertEqual(set(result), set(expected))


if __name__ == '__main__':
    unittest.main()
