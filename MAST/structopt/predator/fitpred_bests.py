from MAST.structopt.tools import get_best
from MAST.structopt.switches import moves_switch, lambdacommamu
from MAST.structopt.generate import gen_pop_box
from MAST.structopt.generate.Individual import Individual
from MAST.structopt.tools import remove_duplicates

def fitpred_bests(pop,Optimizer):
    """Predator function to identify similar structures based on energy and replace one 
    with structure from BESTS List.
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
    count = 0
    while len(newpop) < Optimizer.nindiv:
        try:
            Optimizer.BESTS
        except:
            Optimizer.BESTS=[]
        if len(Optimizer.BESTS) > 0:
            idx = random.choice(range(len(Optimizer.BESTS)))
            newpop.append(Optimizer.BESTS[idx])
            STR+='Predator: Adding in structure from Best List from position = {0} with fitness = {1}\n'.format(idx,Optimizer.BESTS[idx].fitness)
            newindices.append(len(pop)+count)
            count+=1
        else:
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
    return pop, STR
