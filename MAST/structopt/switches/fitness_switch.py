from MAST.structopt.tools import rattle
import logging
try:
    from mpi4py import MPI
except ImportError:
    pass
import pdb

def fitness_switch(input):
    """Performs fitness evaluation on individual"""
    if input[0]==None:
        indiv=0
        rank = MPI.COMM_WORLD.Get_rank()
        stro='Evaluated none individual on {0}\n'.format(rank)
    else:
        Optimizer,indiv=input
        #logger = initialize_logger(Optimizer.loggername)
        logger = logging.getLogger(Optimizer.loggername)
        stro='Evaluating individual {0}\n'.format(indiv.index)
        if Optimizer.rattle_atoms:
            indiv[0]=rattle(indiv)
        scheme = Optimizer.fitness_scheme
        strn = ''
        try:
            exec "from MAST.structopt.fitness.{0} import {0}".format(scheme)
            indiv, strn = eval('{0}(indiv, Optimizer)'.format(scheme))
        except NameError, e:
            print e
            print 'WARNING: Specified fitness not one of the available options. Please check documentation and spelling! Fitness : '+repr(scheme)
            logger.warn('NameError in fitness switch for scheme = {0}. {1}'.format(scheme,e), exc_info=True)
        except Exception, e:
            print 'Error in Fitness fuction : '+repr(scheme)
            print e
            logger.warn('Exception in fitness switch for {0} scheme: {1}'.format(scheme,e),exc_info=True)
        stro+=strn
    return indiv,stro

# if __name__ == "__main__":
#     import sys
#     indiv = sys.argv[1]
#     Optimizer = sys.argv[2]
#     output = fitness_switch([Optimizer,indiv])
#     output[0].write()
#     Optimizer.output.write(output[1])
