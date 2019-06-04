import msprime
import libsequence
import seaborn as sns

D = []

for ts in msprime.simulate(100, num_replicates=100,
                           mutation_rate=25.0,
                           random_seed=42):
    m = libsequence.VariantMatrix.from_TreeSequence(ts)
    ac = m.count_alleles()
    D.append(libsequence.tajd(ac))

sns.distplot(D)
