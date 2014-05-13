from MAST.structopt.tools.eval_energy import eval_energy
from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
import math
import os
from MAST.structopt.tools.find_defects import find_defects
from ase import Atom, Atoms

def sibias2(indiv, Optimizer):
    """Function to calculate total energy fitness of individual and bias away 
    from Si-Interstitials.
    Input:
        indiv = structopt Individual class object to be evaluated
        Optimizer = structopt Optimizer class object
    Output:
        indiv = structopt Individual class object with new fitness.
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
    indiv = outs[2]
    indiv.energy = outs[0]
    stro=outs[3]
    if Optimizer.structure == 'Defect' or Optimizer.structure=='Surface':
        indiv.bulki = outs[1]
    indiv.fitness = indiv.energy
    if abs(indiv.fitness) > Optimizer.energy_cutoff_factor*(len(indiv[0])+len(indiv.bulki)):
        indiv.fitness=10000
        message = 'Warning: Found oddly large energy from Lammps in structure HI={0}'.format(indiv.history_index)
        logger.warn(message)
        print message
        print '    Setting fitness to 10000'
    if math.isnan(indiv.fitness):
        logger.warn('Found NAN energy structure HI={0}'.format(indiv.history_index))
        indiv.fitness=10000
        indiv.energy = 10000
    indc,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
    syms = indc.get_chemical_symbols()
    for sym in syms:
        if sym=='Si':
            message = 'Identified Si atom in defect added energy bias +50eV'
            logger.info(message)
            stro+='    '+message+'\n'
            indiv.fitness+=50
            indiv.energy+=50
    nindiv = Atoms(cell=indiv[0].get_cell(),pbc=True)
    nbulk = indiv.bulki.copy()
    for at in indiv[0]:
        if at.symbol=='Si':
            nbulk.append(at)
        else:
            nindiv.append(at)
    indiv[0] = nindiv
    indiv.bulki = nbulk
    return indiv, stro