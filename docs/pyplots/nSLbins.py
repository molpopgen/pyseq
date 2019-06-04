from collections import namedtuple
import msprime
import libsequence
import numpy as np
import pandas as pd


Datum = namedtuple("Datum", ['dafbin', 'mean', 'sd'])


def get_bin_stats(vm, binsize, minfreq):
    # Get unstandardized nSL for all markers:
    raw = libsequence.nsl(vm, 0)
    # Filter out any NaN values and values
    # not passing our DAF filter
    raw = [(i.nsl, i.core_count)
           for i in raw if i.core_count >= minfreq and np.isfinite(i.nsl)]
    # Convert to numpy array
    # for binning
    raw_array = np.array(raw)
    # Assign data to
    # bins based on DAF count
    bins = np.digitize(raw_array[:, 1], np.arange(0, vm.nsam, binsize))

    # Get mean and SD per bin
    # and return data
    rv = []
    for b in np.unique(bins):
        w = np.where(bins == b)[0]
        stats = raw_array[:, 0][w]
        m = stats.mean()
        sd = stats.std()
        rv.append(Datum(b, m, sd))
    return rv


binsize = 10
minfreq = 3
nreps = 100
nsam = 100
seed = 42*666
results = []
for ts in msprime.simulate(nsam, mutation_rate=25.0,
                           recombination_rate=25.0,
                           num_replicates=nreps,
                           random_seed=seed):
    vm = libsequence.VariantMatrix.from_TreeSequence(ts)
    bins = get_bin_stats(vm, binsize, minfreq)
    results.extend(bins)
df = pd.DataFrame(results, columns=Datum._fields)
g = df.groupby(['dafbin']).mean().reset_index()

# plot DAF bin vs mean nSL statistic in bin
g.plot('dafbin', 'mean', kind='line')
