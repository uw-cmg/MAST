from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
import pdb

def crossover_switch(child1, child2, Optimizer):
    """Functions for selecting and pairing individuals 
    for crossovers"""
    #logger = initialize_logger(Optimizer.loggername)
    logger = logging.getLogger(Optimizer.loggername)
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    if debug:
        s1 = child1[0].copy()
        s2 = child2[0].copy()
        if Optimizer.structure=='Defect':
            s1.extend(child1.bulki.copy())
            s2.extend(child2.bulki.copy())
        write_xyz(Optimizer.debugfile,s1,'First Cx Individual - Pre ')
        write_xyz(Optimizer.debugfile,s2,'Second Cx Individual - Pre')
    passflag = True
    scheme = Optimizer.cx_scheme
    try:
        exec "from MAST.structopt.crossover.{0} import {0}".format(scheme)
        nchild1, nchild2 = eval('{0}(child1, child2, Optimizer)'.format(scheme))
    except NameError, e:
        logger.warning('Specified Crossover not one of the available options. Please check documentation and spelling! Crossover : {0}. {1}'.format(Optimizer.cx_scheme,e), exc_info=True)
        print 'Name Error:', e
        passflag = False
    except Exception, e:
        logger.error('Error in Crossover fuction : {0}. {1}'.format(Optimizer.cx_scheme,e), exc_info=True)
        print 'Exception: ', e
        passflag = False
    if not passflag:
        try:
            from MAST.structopt.crossover.cxtp import cxtp
            nchild1, nchild2 = cxtp(child1, child2, Optimizer)
            passflag = True
        except NameError, e:
            logger.error('Error attempting to load back up cxTP Crossover.  Something very wrong!! {1}'.format(e), exc_info=True)
            print 'CXTP Name Error: ', e
        except Exception, e:
            logger.error('Error in cxtp Crossover {0}'.format(e), exc_info=True)
            print 'Exception: ',e
    if passflag:
        child1=nchild1
        child2=nchild2
        child1.energy=0
        child2.energy=0
        child1.fitness=0
        child2.fitness=0 
        ch1hi=child1.history_index
        child1.history_index = '(' + repr(child1.index) + '+' + repr(child2.index) + ')'
        child2.history_index = '(' + repr(child2.index) + '+' + repr(child1.index) + ')'
    if debug:
        s1 = child1[0].copy()
        s2 = child2[0].copy()
        if Optimizer.structure=='Defect':
            s1.extend(child1.bulki.copy())
            s2.extend(child2.bulki.copy())
        write_xyz(Optimizer.debugfile,s1,'First Cx Individual - Post')
        write_xyz(Optimizer.debugfile,s2,'Second Cx Individual - Post')
    
    return child1, child2
	