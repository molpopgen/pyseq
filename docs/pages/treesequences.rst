.. _treeseqeuences:

Tree Sequences
================================

pylibseq is able to generate instances of :class:`libsequence.VariantMatrix` from a Python object
representing a "tree sequence" :cite:`Kelleher2018-fu`. Currently, the primary method of interacting with these data
structures is via msprime_.

There are two methods to get data into a VariantMatrix.  The first is via numpy arrays as described in
:ref:`variantmatrix`:

.. ipython:: python

    import msprime
    import libsequence

    # Get a TreeSequence from msprime:
    ts = msprime.simulate(100, mutation_rate = 100)

    # Use numpy arrays to make VariantMatrix
    m = libsequence.VariantMatrix(ts.genotype_matrix(), ts.tables.sites.position)

The second method is to call a class method of VariantMatrix:

.. ipython:: python

    m2 = libsequence.VariantMatrix.from_TreeSequence(ts)
    assert(m.data == m2.data)
    assert(m.positions == m2.positions)

The second method is a touch slower than the first, but it uses about half as much memory.  The reason is that getting
the genotype matrix from the TreeSequence requires allocating the entire matrix.  The second method only asks msprime to
generate a 1d numpy array of length `ts.num_samples`, but does so once for each of `ts.num_variants`.

.. _msprime: http://msprime.readthedocs.io

