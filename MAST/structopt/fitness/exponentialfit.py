from MAST.structopt.tools import eval_energy
from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
import math
import os

def exponentialfit(indiv, Optimizer):
    """Fitness function to evaluate the exponential fitness of an individual structure based on its energy.
    Inputs:
        indiv = Individual class object to be evaluated
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Evaluated Individual class object
    """
    #logger = initialize_logger(Optimizer.loggername)
    logger = logging.getLogger(Optimizer.loggername)
    starting = indiv.duplicate()
    cwd = os.getcwd()
    try:
        outs = eval_energy(Optimizer,indiv)
    except Exception, e:
        logger.warn('Error in energy evaluation: {0}'.format(e), exc_info=True)
        stro = 'ERROR: Problem in Energy Evaluation'
        print stro
        print e
        stro += '\n' + repr(e)
        os.chdir(cwd)
        f=open('problem-structures.xyz','a')
        totalsol = indiv[0].copy()
        totalsol.extend(indiv.bulki)
        write_xyz(f,totalsol,data='Starting structure hindex={0}'.format(indiv.history_index))
        indiv.energy = 10000
        f.close()
        print '    Writing structure to problemstructures.xyz file. Structure (hindex) : '+indiv.history_index
        print '    Setting individual energy to 50000.'
        outs = [10000, starting.bulki, starting, stro]
    indiv.energy = outs[0]
    stro=outs[3]
    if Optimizer.structure == 'Defect' or Optimizer.structure=='Surface':
        indiv.bulki = outs[1]
    try:
        f=float(indiv.energy-indiv.tenergymin)/float(indiv.tenergymx-indiv.tenergymin)
    except ZeroDivisionError:
        logger.warn('Population is converging received ZeroDivisionError', exc_info=True)
        stro+='WARNING: Population is converging\n'
        f=float('inf')
    f=math.exp(-f)
    indiv.fitness = f
    return indiv, stro
