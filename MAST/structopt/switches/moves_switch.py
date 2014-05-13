import random
import logging

def moves_switch(indiv, Optimizer):
    """Select mutation to preform on individual from list of options"""
    #logger = initialize_logger(Optimizer.loggername)
    logger = logging.getLogger(Optimizer.loggername)
    scheme = random.choice(Optimizer.mutation_options)
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    try:
       exec "from MAST.structopt.moves.{0} import {0}".format(scheme)
       mutant = eval('{0}(indiv, Optimizer)'.format(scheme))
       mutant.energy = 0
       mutant.fitness = 0
    except NameError, e:
        logger.warning('Specified mutation not one of the available options. Please check documentation and spelling! SKIPPING. Mutation : {0}. {1}'.format(scheme,e), exc_info=True)
        print 'Mutation Name error: ',e
        print scheme
        mutant = indiv
    except Exception, e:
        logger.error('Problem with mutation! SKIPPING. Mutation = {0}. {1}'.format(scheme,e), exc_info=True)
        print 'Mutation Exception: ', e
        mutant = indiv
    return mutant, scheme
