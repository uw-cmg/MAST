from MAST.structopt.tools import eval_energy
from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
import math
import os

def stem_cost(indiv, Optimizer):
    '''Function to calculate STEM_Cost fitness of individual.
    Input:
        indiv = structopt Individual class object to be evaluated
        Optimizer = structopt Optimizer class object
            Must have STEM Calc function attached
    Output:
        indiv = structopt Individual class object with new fitness.
    '''
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
    fit = indiv.energy
    if abs(fit) > Optimizer.energy_cutoff_factor*(len(indiv[0])+len(indiv.bulki)):
        message = 'Warning: Found oddly large energy from Lammps in structure HI={0}'.format(indiv.history_index)
        logger.warn(message)
        print message
        print '    Setting fitness to 10000'
    if math.isnan(fit):
        logger.warn('Found NAN energy structure HI={0}'.format(indiv.history_index))
        indiv.energy = 10000
    
    #Calculate the chisq
    chisq = Optimizer.stemcalc.run(indiv[0])
    #Calculate the coefficient to ensure same order of magnitude
    if Optimizer.stem_coeff == None:
        logger.warn('Not able to identify a stem_coeff')
        print 'Looking for stem_coeff...'
        aflag=True
        alpha = 1.0
        while True:
            value=alpha*chisq
            div=abs(indiv.energy)/value
            if div <1:
                alpha*=0.1
            elif div >10:
                alpha*=10
            else:
                break
        Optimizer.stem_coeff = alpha
    else:
        aflag=False
    #Calculate the fitness of the individual
    indiv.fitness=indiv.energy+Optimizer.stem_coeff*chisq
    #Output values of fitness, energy, rms, and alpha for comparison
    stro+='Individual '+repr(indiv.history_index)+': Fitness '+repr(indiv.fitness)+', Energy '+repr(indiv.energy)+', Chi^2 '+repr(chisq)+', Alpha '+repr(Optimizer.stem_coeff)+'\n'
    return indiv, stro