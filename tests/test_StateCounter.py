import unittest
from libsequence import StateCounter

class test_StateCounter(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.binaryString = '011011'
        self.DNAstring='AGAGCCTT'
        self.nonDNAstring='AZAGCT'
    def testBinaryString1(self):
        c = StateCounter()
        c(self.binaryString[0])
        self.assertEqual(c.zero,1)
        c = StateCounter()
        c(self.binaryString)
        self.assertEqual(c.zero,self.binaryString.count('0'))
        self.assertEqual(c.one,self.binaryString.count('1'))
        self.assertEqual(c.ndna,False)
    def testBinaryStringException(self):
        try:
            c = StateCounter() 
            b = self.binaryString.encode('UTF-8')
            c(b)
        except:
            self.fail("unexpected exception")
    def testDNAstring(self):
        c = StateCounter()
        c(self.DNAstring)
        self.assertEqual(c.nstates(),4)
        self.assertEqual(c.a,self.DNAstring.count('A'))
        self.assertEqual(c.ndna,False)
    def testDNAstring2(self):
        c = StateCounter()
        l=self.DNAstring.lower()
        c(l)
        self.assertEqual(c.nstates(),4)
        self.assertEqual(c.a,l.count('a'))
        self.assertEqual(c.ndna,False)
    def testNonDNAString(self):
        c = StateCounter()
        c(self.nonDNAstring)
        self.assertEqual(c.ndna,True)
