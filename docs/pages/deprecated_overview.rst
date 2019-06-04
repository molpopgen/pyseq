
Example workflow
================

We'll work with pylibseq's wrapper to libsequence's SimData, which is
used to process bi-allele data encoded as 0/1 = ancestral/derived,
respectively

.. ipython:: python

    from __future__ import print_function
    import libsequence

Assigning to an object
----------------------

You may assign from a list of tuples.

Each tuple is a site: (pos:genotypes)

Here, there are 2 sites and a sample size of :math:`n=4`

.. ipython:: python

    rawData1 = [(0.1,'0101'),(0.2,'1010')]
    #We can construct objects straight from these tuples
    sd = libsequence.SimData(rawData1)
    sd.size()
    sd.pos()
    sd.data()
, you can assign from separate list of positions and haplotypes

.. ipython:: python

    rawDataPos = [0.1,0.2]
    rawDataGenos = ['01','10','01','10']
    sd.assign(rawDataPos,rawDataGenos)
    sd.numsites()
    sd.size()
    sd.pos()
    sd.data()

Summary statistics
------------------

Let's calculate some basic summary statistics

See :class:`libsequence.summstats.PolySIM` for more documentation

.. ipython:: python

    #ms 10 1 -s 10 -I 2 5 5 0.05
    rawDataPos=[0.0997, 0.2551, 0.3600, 0.4831, 0.5205, 0.5668, 0.5824, 0.6213, 0.7499, 0.9669]
    rawDataGenos=['0000001010',
                  '1000000011',
                  '1000001010',
                  '1000000010',
                  '1000001010',
                  '1111010100',
                  '1111010100',
                  '1111110100',
                  '1111010100',
                  '1111010100']
    sd.assign(rawDataPos,rawDataGenos)
    ps = libsequence.PolySIM(sd)
    ps.thetapi()
    ps.thetaw()
    ps.tajimasd()

Filtering sites
------------------

You may also filter sites from your variation tables using :func:`libsequence.removeColumns`, which operates on
objects of type :class:`libsequence.StateCounter`.

Our data look like this right now:

.. ipython:: python

    print(sd)

Let's remove derived singletons:

.. ipython:: python

    sd2 = libsequence.removeColumns(sd,lambda x : x.one != 1)
    print(sd2)

Let's remove all singletons:

.. ipython:: python

    sd3 = libsequence.removeColumns(sd,lambda x: x.one !=1 and x.zero != 1)
    print(sd3)

.. note::
    Yes, it is odd that the column removal function removes sites for which the lambda
    returns False.  I'll fix that in a future release, which requires an upstream change
    to libsequence.
    
Sliding windows
---------------

.. ipython:: python

    w = libsequence.Windows(sd,window_size=0.1,step_len=0.05,starting_pos=0.,ending_pos=1.0)
    len(w)
    for i in range(len(w)):
        #Each window is a simData
        wi = w[i]
        pswi = libsequence.PolySIM(wi)
        print(pswi.thetaw())


Linkage disequilibrium
----------------------

The function ``libsequence.ld`` returns pairwise LD stats as a
``list`` of ``dict``\ s. The return value is easily coerced into a
``pandas.DataFrame``.

:math:`F_{ST}`
--------------

Let's pretend that our data are from two demes of sizes n/2 each.

Note that most flavors of :math:`F_{ST}` are very similar to one
another. See Charlesworth, B. (1998) Mol. Biol. Evol. 15(5): 538-543 for
a great overview.

.. ipython:: python

    import libsequence
    sd.size()
    f = libsequence.Fst(sd,[5,5])

    #Hudson, Slatkin, and Maddison's FST:
    f.hsm()

    #Slatkin's
    f.slatkin()

    #Hudson, Boos, and Kaplan, which is also Nei's Gst:
    f.hbk()

    #Positions of snps shared b/w demes 0 and 1
    f.shared(0,1)

    #Positions of private mutations in deme 0 and 1:
    f.priv(0,1)

    #Positions of fixed differences between demes 0 and 1:
    f.fixed(0,1)
