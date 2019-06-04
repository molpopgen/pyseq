import unittest
import msprime
import numpy as np
import libsequence


class testAlleleCountMatrix(unittest.TestCase):
    @classmethod
    def setUp(self):
        ts = msprime.simulate(10, mutation_rate=10, random_seed=42)
        self.vm = libsequence.VariantMatrix.from_TreeSequence(
            ts)
        self.ac = libsequence.AlleleCountMatrix(self.vm)

    def testSubset(self):
        w = self.vm.window(0.25, 0.3)
        acw = libsequence.AlleleCountMatrix(w)

        pa = np.array(self.vm.positions)
        p = np.where((pa >= 0.25) & (pa <= 0.3))[0]
        self.assertEqual(len(p), acw.nrow)

        ac_slice = self.ac[p.min():p.max()+1:1]
        self.assertTrue(np.array_equal(np.array(ac_slice),
                                       np.array(acw)))

    def testArbitrarySubset(self):
        indexes = [1, 3, 4, 11, 15]
        acw = self.ac[indexes]
        for i, j in zip(indexes, range(len(indexes))):
            row_i = self.ac.row(i)
            row_j = acw.row(j)
            self.assertTrue(all(k == l for k, l in zip(row_i, row_j)) is True)


if __name__ == "__main__":
    unittest.main()
