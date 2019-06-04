import unittest
import libsequence
import numpy as np


class testVariantMatrix(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.data = [0, 1, 1, 0, 0, 0, 0, 1]
        self.pos = [0.1, 0.2]
        self.m = libsequence.VariantMatrix(self.data, self.pos)

    def testConstruct(self):
        self.assertEqual(self.m.data, self.data)
        self.assertEqual(self.m.positions, self.pos)
        self.assertEqual(self.m.nsam, 4)
        self.assertEqual(self.m.nsites, 2)

    def testModifyCppDataViaNumpy(self):
        d = np.array(self.m, copy=False)
        d[0][2] = 4
        self.assertEqual(self.m.data[2], 4)

    def testFilterSites(self):
        from collections import Counter

        def is_singleton(x):
            c = Counter([i for i in x])
            if len(c) == 2:
                for i in c.most_common():
                    if i[1] == 1:
                        return True
            return False
        self.assertEqual(self.m.nsites, 2)
        libsequence.filter_sites(
            self.m, is_singleton)
        self.assertEqual(self.m.nsites, 1)

    def testPickle(self):
        import pickle
        p = pickle.dumps(self.m, -1)
        up = pickle.loads(p)
        self.assertEqual(up.data, self.m.data)
        self.assertEqual(up.positions, self.m.positions)
        self.assertEqual(up.nsam, self.m.nsam)
        self.assertEqual(up.nsites, self.m.nsites)


class testColumnViews(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.data = [0, 1, 1, 0, 0, 0, 0, 1]
        self.pos = [0.1, 0.2]
        self.m = libsequence.VariantMatrix(self.data, self.pos)

    def testIterateColumns(self):
        for i in range(self.m.nsites):
            c = self.m.site(i)
            self.assertEqual(len(c), self.m.nsam)
        s = [j for j in self.m.site(0)]
        self.assertEqual(s, self.data[:self.m.nsam])
        s = [j for j in self.m.site(1)]
        self.assertEqual(s, self.data[self.m.nsam:])

    def testIterateRows(self):
        d = np.array(self.data, dtype=np.int8)
        d = d.reshape((self.m.nsites, self.m.nsam))
        for i in range(self.m.nsam):
            c = self.m.sample(i)
            s = np.array([j for j in c], dtype=np.int8)
            self.assertTrue(np.array_equal(s, d[:, i]))


class testCreationFromNumpy(unittest.TestCase):
    def testConstruct(self):
        d = np.array([0, 1, 1, 0, 0, 0, 0, 1], dtype=np.int8).reshape((2, 4))
        pos = np.array([0.1, 0.2])
        m = libsequence.VariantMatrix(d, pos)
        ma = np.array(m)
        self.assertTrue(np.array_equal(np.sum(ma, axis=0), np.sum(d, axis=0)))
        self.assertTrue(np.array_equal(np.sum(ma, axis=1), np.sum(d, axis=1)))


class testDataFromMsprime(unittest.TestCase):
    def testDirectConversionOfData(self):
        """
        Test is only run is msprime is
        available for import
        """
        try:
            import msprime
            try:
                ts = msprime.simulate(10, mutation_rate=10, random_seed=666)
                gm = ts.genotype_matrix()
                pos = np.array([i.position for i in ts.sites()])
                m = libsequence.VariantMatrix(gm, pos)
                # If the data conversion is correct, the row and
                # column sums must match
                ma = np.array(m, copy=False)
                self.assertTrue(np.array_equal(
                    np.sum(ma, axis=0), np.sum(gm, axis=0)))
                self.assertTrue(np.array_equal(
                    np.sum(ma, axis=1), np.sum(gm, axis=1)))
            except:
                raise
        except:
            pass


if __name__ == "__main__":
    unittest.main()
