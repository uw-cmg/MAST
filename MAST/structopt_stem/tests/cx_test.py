
def cx_test(ind1, ind2, A):
    '''Crossover function unit tests
    Key to test is:
        - imports and executes without errors
        - preserves number of atoms
        - preserves concentration of atoms
        - maintains relative positions of atoms
        - produces two distinct individuals
    '''
    print 'Beginning unit testing of crossovers'
    try:
        from MAST.structopt_stem.crossover.clustbx import clustbx
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = clustbx(i1,i2,A)
        print 'clustbx crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: clustbx crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.cxTP import cxTP
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = cxTP(i1,i2,A)
        print 'cxTP crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: cxTP crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.cxTPA import cxTPA
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = cxTPA(i1,i2,A)
        print 'cxTPA crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: cxTPA crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.cxTPC import cxTPC
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = cxTPC(i1,i2,A)
        print 'cxTPC crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: cxTPC crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.NewClus import NewClus
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = NewClus(i1,i2,A)
        print 'NewClus crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: NewClus crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.randalloybox import randalloybox
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = randalloybox(i1,i2,A)
        print 'randalloybox crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: randalloybox crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.rotct_rand_clus import rotct_rand_clus
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = rotct_rand_clus(i1,i2,A)
        print 'rotct_rand_clus crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: rotct_rand_clus crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.rotct_rand import rotct_rand
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = rotct_rand(i1,i2,A)
        print 'rotct_rand crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: rotct_rand crossover test FAILED'
        print e
        pass
    try:
        from MAST.structopt_stem.crossover.rotct import rotct
        i1 = ind1.duplicate()
        i2 = ind2.duplicate()
        nc1, nc2 = rotct(i1,i2,A)
        print 'rotct crossover test SUCCESSFUL'
    except Exception, e:
        print 'ERROR: rotct crossover test FAILED'
        print e
        pass
    print 'Crossover Function Testing Complete'
    return
