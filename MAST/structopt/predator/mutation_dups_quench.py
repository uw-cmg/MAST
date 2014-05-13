from MAST.structopt.fingerprinting import get_fingerprint
from MAST.structopt.tools import get_best
from MAST.structopt.tools import remove_duplicates
from MAST.structopt.switches import selection_switch, moves_switch, lambdacommamu
import random

def mutation_dups_quench(pop, Optimizer):
    """Predator function that removes individuals based on fitness and mutates replacements
    Also quenches top individuals
    """
    fitlist = [one.fitness for one in pop]
    nfitlist, nindices = remove_duplicates(fitlist, Optimizer.demin)
    STR = ''
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
        indiv = random.choice(otherlist).duplicate()
        indiv, scheme = moves_switch(indiv,Optimizer)
        indiv.energy = 1000
        indiv.fitness = 1000
        newpop.append(indiv)
        STR+='Predator: Adding mutated duplicates to new pop history='+indiv.history_index+'\n'
        nindices.append(indiv.index)
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
        pop, str1 = lambdacommamu.lambdacommamu(newpop, Optimizer)
        STR+=str1
    else:
        pop = selection_switch(newpop, Optimizer.nindiv, Optimizer.natural_selection_scheme, Optimizer)
    pop = get_best(pop,len(pop))
    
    if Optimizer.genrep >10:
        from MAST.structopt.moves.quench import quench
        import os
        olammpsvar = os.environ['LAMMPS_COMMAND']
        try:
            from mpi4py import MPI
            if '-n' in olammpsvar:
                lcommand = olammpsvar.split('-n')
                lcommand[1]=lcommand[1].split()
                nproc = MPI.COMM_WORLD.Get_size()
                os.environ['LAMMPS_COMMAND'] = '{0}-n {1} {2}'.format(lcommand[0],nproc,lcommand[1][1])
        except:
            pass
        oqns2 = Optimizer.quench_n_steps_2
        Optimizer.quench_n_steps_2 = 100000
        opar = Optimizer.parallel
        Optimizer.parallel = False
        for i in range(3):
            pop[i] = quench(pop[i],Optimizer)
        Optimizer.quench_n_steps_2 = oqns2
        os.environ['LAMMPS_COMMAND'] = olammpsvar
        Optimizer.parallel = opar
    return pop, STR
