from MAST.structopt.tools import eval_energy
from MAST.structopt.fitness import totalenfit
from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
import math
import os

def formationenergy(indiv, Optimizer):
    """Function to calculate formation energy fitness of individual.
    Input:
        indiv = structopt Individual class object to be evaluated
        Optimizer = structopt Optimizer class object
    Output:
        indiv = structopt Individual class object with new fitness.
    """
    if Optimizer.structure != 'Defect':
        if Optimizer.generation <1:
            print 'WARNING: Formation Energy fitness calculator only available for Defect structure simulation.'
            print '     Using total energy fitness calculator instead'
        indiv = totalenfit(indiv,Optimizer)
        stro = ''
    else:
        logger = logging.getLogger(Optimizer.loggername)
        #logger = initialize_logger(Optimizer.loggername)
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
        indiv.bulki = outs[1]
        #fit = (indiv.purebulkenpa*(indiv[0].get_number_of_atoms()+indiv.bulki.get_number_of_atoms()) - indiv.energy)/(indiv[0].get_number_of_atoms()+indiv.bulki.get_number_of_atoms()-indiv.natomsbulk)
        solid=Atoms()
        solid.extend(indiv[0])
        solid.extend(indiv.bulki)
        fit=indiv.energy
        passflag = True
        if abs(fit) > Optimizer.energy_cutoff_factor*(len(indiv[0])+len(indiv.bulki)):
            fit=10000
            message = 'Warning: Found oddly large energy from Lammps in structure HI={0}'.format(indiv.history_index)
            logger.warn(message)
            print message
            print '    Setting fitness to 10000'
            passflag=False
        if math.isnan(fit):
            logger.warn('Found NAN energy structure HI={0}'.format(indiv.history_index))
            fit=10000
            passflag=False
            indiv.energy = 10000
        if passflag:
            for sym,c,m,u in Optimizer.atomlist:
                nc=len([atm for atm in solid if atm.symbol==sym])
                fit-= float(nc)*float(u)
        indiv.fitness=fit
    return indiv, stro
