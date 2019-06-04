.. _variantmatrix:

The VariantMatrix
===============================

Variation data are stored in instances of :class:`libsequence.VariantMatrix`.  This class stores state
data as signed 8-bit integers and position data as floats (C/C++ doubles).  Missing data are encoded as negative
numbers.  You can use any negative number for missing data except the minimum possible value of an 8-bit integer, which
is reserved for internal use as a "mask" value when filtering data out of variant matrices.

.. ipython:: python

    import libsequence

    # You can construct variant matrices from lists
    states = [0, 1, 1, 0, 0, 0, 0, 1]
    pos = [0.1, 0.2]

    m = libsequence.VariantMatrix(states, pos)
    print(m.nsites)
    print(m.nsam)

A variant matrix supports Python's buffer protocol, meaning that it can be converted directly to a numpy array:

.. ipython:: python

    import numpy as np

    ma = np.array(m, copy=False)
    print(ma)

Our `ma` object is a thin wrapper to the underlying memory allocated by C++, so we can change the data:

.. ipython:: python

    x = ma[0,2]
    ma[0,2] = -1
    print(ma)
    ma[0,2] = x
    print(ma)

You can construct from numpy arrays:

.. ipython:: python

    m = libsequence.VariantMatrix(np.array(states, dtype=np.int8).reshape((2,4)), np.array(pos))
    print(m.nsites)
    print(m.nsam)

You can access the site and sample data via loops:

.. ipython:: python

    for i in range(m.nsites):
        s = m.site(i)
        print([i for i in s])

    for i in range(m.nsam):
        s = m.sample(i)
        print([i for i in s])

The types involved in the last example are:

.. ipython:: python

    print(type(m.site(0)))
    print(type(m.sample(0)))

Counting the states at a site
-------------------------------------

The most straightforward way to get the allele counts at all sites is via
:class:`libsequence.AlleleCountMatrix`:

.. ipython:: python

    import msprime

    ts = msprime.simulate(10, mutation_rate=25, random_seed=666)
    m = libsequence.VariantMatrix.from_TreeSequence(ts)
    ac = m.count_alleles()
    print(np.array(ac)[:5])

    # Confirm that the counts are the same as 
    # what msprime thinks:
    vi = ts.variants()
    for i in range(5):
        v = next(vi)
        ones = np.count_nonzero(v.genotypes)
        print(len(v.genotypes)-ones, ones)

    # These count objects are sliceable...
    print(np.array(ac[1:ac.nrow:25]))

    # ...and indexable via lists
    print(np.array(ac[[0,1,2,3,4]]))

The allele count data are stored in order of allele label, starting with zero.  The sum
of allele counts at a site is the sample size at that site.

:class:`libsequence.StateCounts` provides a means to generate allele counts
on-demand for a site:

.. ipython:: python

    c = libsequence.StateCounts()
    # These objects are callable classes:
    c(m.site(0))
    print(c.counts[:3])
    # The sample size at this site
    print(c.n)
    # c is iterable...
    for i in c:
        if i > 0:
            print(i)
    #...and indexable...
    for i in range(len(c)):
        if c[i] > 0:
            print(i,c[i])
    #...and supports the buffer protocol
    ca = np.array(c)
    nonzero_states = np.where(ca > 0)
    print(nonzero_states[0])
    ca[nonzero_states[0]]
            
By convention, missing data affects the sample size at a site:

.. ipython:: python

    ma = np.array(m, copy=False)
    ma[0,2] = -1
    c(m.site(0))
    print(c.counts[:3])
    print(c.n)

    # restore our object
    ma[0,2] = x

    print(c.refstate)

You may specify a reference state when counting.  Depending on the analysis, that may mean a literal reference genome
state, an ancestral state, a minor allele state, etc.

.. ipython:: python

    # Above, no reference state was specified, 
    # so it is considered missing:

    print(c.refstate)

    # Let's let 0 be the reference state:
    c = libsequence.StateCounts(refstate = 0)
    c(m.site(0))
    print(c.counts[:3])
    print(c.refstate)

You may get all of the counts at all sites in three different ways:

.. ipython:: python

    # Without respect to reference state
    lc = libsequence.process_variable_sites(m)
    for i in lc[:5]:
        print(i.counts[:2], i.refstate)
    
    # With a single reference state for all sites
    lc = libsequence.process_variable_sites(m, 0)
    for i in lc[:5]:
        print(i.counts[:2], i.refstate)

    # With a reference specified state for each site
    rstats = [0 for i in range(m.nsites)]
    rstats[0:len(rstats):2] = [1 for i in range(0,len(rstats),2)] 
    lc = libsequence.process_variable_sites(m, rstats)
    for i in lc[:5]:
        print(i.counts[:2], i.refstate)


Encoding missing data
-------------------------------------

.. ipython:: python
    :okexcept:

    # This is the value of the reserved state:
    print(libsequence.VariantMatrix.mask)

    # Attempting to construct an object with this
    # value is allowed, but is an error.
    # Downstream analyses will see this and raise exceptions.

    x = libsequence.VariantMatrix([0, 1, libsequence.VariantMatrix.mask, 2], [0.2, 0.5])
    print(x.data)

    # For example:
    c(x.site(1))

Filtering VariantMatrix data
-------------------------------------

You may remove sites and/or samples via the application of functions written in Python.  To filter sites, a function
must take the return value of :func:`libsequence.variant_matrix.VariantMatrix.site` as an argument:

.. ipython:: python

    class RemoveNonRefSingletons(object):
        def __init__(self):
            # Treat 0 as the reference state
            self.__c = libsequence.StateCounts(0)
        def __call__(self, x):
            self.__c(x)
            n=np.array(self.__c, copy=False)
            singletons = np.where(n == 1)
            if len(singletons[0])>0:
                return True
            return False

    # Copy our data
    m2 = libsequence.VariantMatrix(m.data, m.positions)

    rv = libsequence.filter_sites(m2, RemoveNonRefSingletons())
    print(np.array(m))
    print(np.array(m2))

    # This is the number of sites removed:
    print(rv)

Performance tip: I wrote the callable as a class so that a StateCounts could be stored as member data.  The reason is
that :attr:`libsequence.StateCounts.counts` is a buffer whose memory is re-used for each call.  Thus,
storing an instance saves repeated memory allocation/deallocation events for each site.

Similarly, we can remove samples:

.. ipython:: python

    # Treat 0 as the reference state
    def remove_all_ref_samples(x):
        if all([i==0 for i in x]):
            return True
        return False

    m2 = libsequence.VariantMatrix(m.data, m.positions)

    rv = libsequence.filter_haplotypes(m2, remove_all_ref_samples)

    print(rv)
    print(np.array(m))
    print(np.array(m2))
