from MAST.structopt.tools import get_best
from MAST.structopt.switches import moves_switch, lambdacommamu
from MAST.structopt.generate import gen_pop_box
from MAST.structopt.generate.Individual import Individual
from MAST.structopt.tools import remove_duplicates

def fitpred_new(pop,Optimizer):
    """Predator function to identify similar structures based on energy and replace one with new structure.
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
        if Optimizer.structure=='Defect' or Optimizer.structure=='Cluster':
            ind=gen_pop_box(Optimizer.atomlist,Optimizer.size)
        elif Optimizer.structure=='Crystal':
            outts=gen_pop_box(Optimizer.atomlist,Optimizer.size,Optimizer.cell_shape_options)
            ind=outts[0]
        elif Optimizer.structure=='Surface':
            mutopto=Optimizer.mutation_options
            Optimizer.mutation_options=['Lattice_Alteration_rdrd']
            topind=random.choice(pop)[0].copy()
            ind, scheme = moves_switch(topind,Optimizer)
            Optimizer.mutation_options=mutopto
        individ=Individual(ind)
        #CHECK THIS LATER!! MAY NEED TO ADD MORE PROPERTIES!!
        individ.energy=1000
        individ.fitness=1000
        newpop.append(individ)
        STR+='Predator: Adding mutated duplicates to new pop history='+individ.history_index+'\n'
        nindices.append(individ.index)
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
