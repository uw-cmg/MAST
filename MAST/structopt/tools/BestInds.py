import os
from operator import attrgetter
from ase import Atom, Atoms
from MAST.structopt.inp_out.write_xyz import write_xyz
import logging
try:
    from mpi4py import MPI
except:
    pass

def BestInds(pop, bests, Optimizer, writefile=False, fname=None):
    """Function to keep track of the best individuals throughout an optimization
    Inputs:
        pop = list of Individual class structures to be compared
        bests = list of previously obtained best structures
        Optimizer = Optimizer class structure that contains key parameters
        writefile = True/False boolean to tell function whether or not to output structures
        fname = filename to output structures
    Outputs:
        bests = list of Individual class structures updated with new options
    """
    #logger = initialize_logger(Optimizer.loggername)
    logger = logging.getLogger(Optimizer.loggername)
    logger.info('best_inds_list recieved population with length = {0} and best list with length = {1} for generation {2}'.format(
        len(pop),len(bests),Optimizer.generation))
    pop = sorted(pop, key=attrgetter('fitness'))
    tol=Optimizer.demin
    if len(bests)!=0:
        bfits=[ind.fitness for ind in bests]
        #Identify the fitnesses in the population
        popfits=[ind.fitness for ind in pop]
        #Initialize a list for individuals to not be added to best list
        indices=[]
        for i in range(len(popfits)):
            for j in range(len(bfits)):
                if popfits[i]==bfits[j]:
                    indices.append(i)
                if popfits[i] < bfits[j]:
                    diff = abs(bfits[j] - popfits[i])
                    diff2 = abs(bfits[j-1] - popfits[i])
                    if diff > tol:
                        if diff2 > tol:
                            if i not in indices:
                                bests.append(pop[i])
                                logger.info('Best list found new structure for List HI = {0}'.format(pop[i].history_index))
                                indices.append(i)
                                b=bfits[0:j]
                                b.append(pop[i].fitness)
                                b.extend(bfits[j::])
                                bfits=b
                        else:
                            indices.append(i)
        bests = sorted(bests, key=attrgetter('fitness'))
        if len(bests) < Optimizer.number_of_bests:
            bmax = bests[len(bests)-1].fitness
            for i in range(len(popfits)):
                if popfits[i] > bmax:
                    diff = abs(popfits[i]-bmax)
                    if diff > tol:
                        bests.append(pop[i])
                        logger.info('Appending best list with lower energy structure')
                        bmax = pop[i].fitness
        else:
            logger.info('Removing extra {0} structures from best list'.format(len(bests)-Optimizer.number_of_bests))
            bests = bests[0:Optimizer.number_of_bests]
        bfits=[ind.fitness for ind in bests]
    else:
        bests = []
        bfits = []
        for one in pop:
            fit = one.fitness
            if len(bfits) == 0:
                bests.append(one)
                logger.info('Adding new structure to best list HI = {0}'.format(one.history_index))
                bfits.append(fit)
            else:
                diff = abs(bfits[-1]-fit)
                if diff > tol:
                    bests.append(one)
                    bfits.append(fit)
        if len(bests) > Optimizer.number_of_bests:
            logger.info('Removing extra {0} structures from best list'.format(len(bests)-Optimizer.number_of_bests))
            bests = bests[0:Optimizer.number_of_bests]
    if writefile==True:
        if fname==None:
            try:
                rank = MPI.COMM_WORLD.Get_rank()
            except:
                rank = 0
            path = os.path.join(os.getcwd(),'{0}-rank{1}'.format(Optimizer.filename,rank))
            logger.info('Writing Best lists to {0}'.format(path))
            fxyzname=os.path.join(path,'Bests-{0}.xyz'.format(Optimizer.filename))
            ffitname=os.path.join(path,'Bests-fitnesses-{0}.txt'.format(Optimizer.filename))
        else:
            fxyzname='{0}.xyz'.format(fname)
            ffitname='{0}-fitnesses.txt'.format(fname)
        bestfile=open(fxyzname,'w')
        bestfitfile=open(ffitname,'w')
        for one in bests:
            if Optimizer.structure=='Defect':
                cbulk=Atoms()
                cbulk.extend(one[0])
                cbulk.extend(one.bulki)
                write_xyz(bestfile,cbulk,one.energy)
            else:
                write_xyz(bestfile,one[0],one.energy)
            bestfitfile.write(repr(one.fitness)+'\n')
        bestfile.close()
        bestfitfile.close()
    return bests
    