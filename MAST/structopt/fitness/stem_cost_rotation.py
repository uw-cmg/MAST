from MAST.structopt.tools import eval_energy
from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
import math
import os

def stem_cost_rotation(indiv, Optimizer):
    '''Function to calculate STEM_Cost fitness of individual.
    Input:
        indiv = structopt Individual class object to be evaluated
        Optimizer = structopt Optimizer class object
            Must have STEM Calc function attached
    Output:
        indiv = structopt Individual class object with new fitness.
    '''
    if 'FIT' in Optimizer.debug:
        debug=True
    else:
        debug=False
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
    chisq0 = Optimizer.stemcalc.run(indiv[0])
    #Identify the 12 <110> directions
    indatoms = indiv[0].copy()
    zdir_list = get_z_directions(indatoms)
    chisqlist = []
    for direction in zdir_list:
        #align atoms so z direction is is in direction
        indatoms.rotate_euler(direction[0], direction[1], direction[2])
        chisqn, phi = find_optimial_phi_degree(atoms, debug=debug)
        chisqlist.append([chisqn,phi,direction])
    minchi = min([chi for chi, phi,direction in chisqlist])
    minphi = [one[1] for one in chisqlist if one[0]==minchi][0]
    minz = [one[2] for one in chisqlist if one[0]==minchi][0]
    atomsnew = indiv[0].copy().rotate_euler(minz[0],minz[1],minz[2])
    atomsnew.rotate(minphi)
    indiv[0]=atomsnew.copy()
    #Calculate the fitness of the individual
    indiv.fitness=indiv.energy+Optimizer.stem_coeff*minchi
    #Output values of fitness, energy, rms, and alpha for comparison
    stro+='Individual '+repr(indiv.history_index)+': Fitness '+repr(indiv.fitness)+', Energy '+repr(indiv.energy)+', Chi^2 '+repr(chisq)+', Alpha '+repr(Optimizer.stem_coeff)+'\n'
    return indiv, stro

def get_z_directions(atoms):
    """Funciton to calculate the z directions from an atoms structure
    Input:
        atoms = atoms structure
    Output:
        z_dirs = list of angles for euler rotation representing 12 z directions"""
    #Do function
    return z_dirs

def find_optimal_psi_degree(atoms, debug):
    chisq = Optimizer.stemcalc.run(indiv[0])
    if debug:
        print angle
    return chisq, atomsn
