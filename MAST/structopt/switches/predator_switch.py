from MAST.structopt.tools import get_best
from MAST.structopt.switches import selection_switch
from MAST.structopt.switches import lambdacommamu
from MAST.structopt.tools import remove_duplicates
import logging
import random
import pdb

def predator_switch(pop,Optimizer):
    """Function for removing individuals from the population"""
    #logger = initialize_logger(Optimizer.loggername)
    logger = logging.getLogger(Optimizer.loggername)
    scheme = Optimizer.predator
    logger.info('Applying predator to population with initial size = {0}'.format(len(pop)))
    STR = 'PREDATOR\n'
    try:
       exec "from MAST.structopt.predator.{0} import {0}".format(scheme)
       pop, STR = eval('{0}(pop, Optimizer)'.format(scheme))
       passflag = True
    except NameError, e:
        logger.warning('Specified predator not one of the available options. Please check documentation and spelling! Predator : {0}. {1}'.format(scheme,e), exc_info=True)
        passflag = False
        STR+='Specified predator not one of the available options. Please check documentation and spelling! Predator : '+repr(scheme)
        STR+=repr(e)+'\n'
    except Exception, e:
        logger.error('ERROR: Issue in Predator Scheme. Predator = {0}. {1}'.format(scheme,e), exc_info=True)
        print 'ERROR: Issue in Predator Scheme. Predator = '+repr(scheme)
        print e
        passflag = False
        STR+=''
    if not passflag:
        logger.warning('Issue in predator. Attempting basic Fitpred')
        fitlist = [one.fitness for one in pop]
        nfitlist, nindices = remove_duplicates(fitlist, Optimizer.demin)
        STR += 'Issue in predator. Attempting basic Fitpred\n'
        newpop = []
        if len(nfitlist) != len(fitlist):
            STR+='Predator: Removed total of '+repr(len(fitlist)-len(nfitlist))+' from population\n'
        otherlist = []
        for i in range(len(pop)):
            if i not in nindices:
                STR+='Predator: Removed '+repr(pop[i].history_index)+'\n'
                otherlist.append(pop[i])
            else:
                newpop.append(pop[i])
        while len(newpop) < Optimizer.nindiv:
            STR+='Predator: Adding duplicates back\n'
            choice = random.choice(otherlist)
            if choice.index not in nindices:
                newpop.append(choice)
                nindices.append(choice.index)
        nindices.sort()
        if Optimizer.natural_selection_scheme=='fussf':
            for ind in newpop:
                if ind.fingerprint == 0:
                    ind.fingerprint = get_fingerprint(Optimizer,ind,Optimizer.fpbin,Optimizer.fpcutoff)
        if 'lambda,mu' in Optimizer.algorithm_type:
            try:
                mark = [ index for index,n in enumerate(nindices) if n > Optimizer.nindiv-1][0]
            except:
                mark = Optimizer.nindiv
            Optimizer.mark = mark
            pop,str = lambdacommamu.lambdacommamu(newpop, Optimizer)
            STR+=str
        else:
            pop = selection_switch(newpop, Optimizer.nindiv, Optimizer.natural_selection_scheme, Optimizer)
        pop = get_best(pop,len(pop))
    Optimizer.output.write(STR)
    return pop

