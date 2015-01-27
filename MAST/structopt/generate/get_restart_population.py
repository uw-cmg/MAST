from MAST.structopt.inp_out import read_xyz
from MAST.structopt.generate.surface import get_surface_restart_indiv
from MAST.structopt.generate.defect import get_defect_restart_indiv
from MAST.structopt.generate.crystal import get_crystal_restart_indiv
from MAST.structopt.generate.Individual import Individual
import logging

def get_restart_population(Optimizer):
    """
    Function to generate a population from a folder containing existing structures.
    Inputs:
        Optimizer = structopt Optimizer class object
    Outputs:
        pop = List of structopt Individual class objects containing existing structures.
    """
    logger = logging.getLogger(Optimizer.loggername)
    index1 = 0
    Optimizer.output.write('Loading structures from old run\n')
    pop = []
    for i in range(Optimizer.nindiv):
        successflag = False
        try:
            indiv = read_xyz(Optimizer.files[i].name)
            successflag = True
        except IOError,e:
            logger.error('Not enough files in restart to generate population. Resetting nindiv to {0}'.format(i-1),exc_info=True)
            Optimizer.output.write('WARNING: Not enough files in restart to generate population\n')
            Optimizer.nindiv=i-1
            Optimizer.output.write('Resetting nindiv = {0}\n'.format(Optimizer.nindiv))
            Optimizer.output.flush()
            break
        except Exception,e:
            Optimizer.output.write('WARNING: Trouble reading file: {0}'.format(Optimizer.files[i].name),exc_info=True)
            Optimizer.output.write('Error: {0}'.format(e))
            Optimizer.output.flush()
            Optimizer.nindiv-=1
        if successflag:
            Optimizer.output.write('Found good individual = {0}\n'.format(indiv))
            if Optimizer.structure == 'Defect':
                individ = get_defect_restart_indiv(Optimizer,indiv)
            elif Optimizer.structure == 'Surface':
                individ = get_surface_restart_indiv(Optimizer, indiv)
            else:
                indiv.set_cell([Optimizer.size,Optimizer.size,Optimizer.size])
                individ = Individual(indiv)
            individ.index = index1
            if Optimizer.genealogy: 
            	individ.history_index = repr(index1)
            pop.append(individ)
            index1 = index1+1
    if len(pop) == 0:
        raise RuntimeError('Unable to load any structures for Restart')
    # Generate new atomlist concentrations based on cluster+box
    if Optimizer.structure == 'Defect':
        if Optimizer.alloy:
            concents=[]
            for ind in pop:
                cs=[]
                for sym,c,u,m in Optimizer.atomlist:
                    sylen=[atm for atm in ind[0] if atm.symbol==sym]
                    cs.append(len(sylen))
                concents.append(cs)
            natmlist=[0]*len(Optimizer.atomlist)
            for i in range(len(concents[0])):
                alls=[cs[i] for cs in concents]
                avgall=int(sum(alls)/len(alls))
                if avgall>=Optimizer.atomlist[i][1]:
                    natmlist[i]=(Optimizer.atomlist[i][0], avgall,\
                    Optimizer.atomlist[i][2],Optimizer.atomlist[i][3])
                else:
                    natmlist[i]=(Optimizer.atomlist[i][0], Optimizer.atomlist[i][1],\
                    Optimizer.atomlist[i][2],Optimizer.atomlist[i][3])
        else:
            natmlist=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                atms1=[inds for inds in pop[0][0] if inds.symbol==Optimizer.atomlist[i][0]]
                natmlist[i]=(Optimizer.atomlist[i][0], len(atms1),Optimizer.atomlist[i][2],\
                Optimizer.atomlist[i][3])
        Optimizer.atomlist=natmlist
        Optimizer.output.write('\n\nNew atomlist concentrations based on cluster+box = {0}\n'.format(
            Optimizer.atomlist))
    return pop

