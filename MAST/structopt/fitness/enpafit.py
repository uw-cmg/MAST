from MAST.structopt.tools import eval_energy
from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
import math
import os

def enpafit(indiv,Optimizer):
    """Fitness function to evaluate the energy per atom of an individual structure.
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
    fit = indiv.energy
    passflag = True
    if abs(fit) > Optimizer.energy_cutoff_factor*(len(indiv[0])+len(indiv.bulki)):
        fit=10000
        message = 'Warning: Found oddly large energy from Lammps in structure HI={0}'.format(indiv.history_index)
        logger.warn(message)
        print message
        print '    Setting fitness to 10000'
        passflag = False
    if math.isnan(fit):
        logger.warn('Found NAN energy structure HI={0}'.format(indiv.history_index))
        fit=10000
        passflag = False
        indiv.energy = 10000
    if passflag:
        if Optimizer.structure == 'Defect' or Optimizer.structure=='Surface':
            indiv.fitness = indiv.energy/(indiv[0].get_number_of_atoms()+indiv.bulki.get_number_of_atoms())
        else:
            indiv.fitness = indiv.energy/indiv[0].get_number_of_atoms()
    else:
        indiv.fitness = fit
    return indiv, stro
