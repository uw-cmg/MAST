from MAST.structopt.fingerprinting import get_fingerprint
from MAST.structopt.switches import selection_switch
from MAST.structopt.tools import remove_duplicates
import math

def adapting(pop, Optimizer):
    """Function to provide an adapting fitness function for GA evaluation
    Input:
        pop = population consisting of list of Individual Class objects to be evaluated
        Optimizer = Optimizer class object with fitness parameters
    Output:
        pop = new population updated based on fitness evaluation
    *** needs work ***
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
        Optimizer.output.write('Predator: Adding duplicates back')
        indiv = random.choice(otherlist)
        newpop.append(indiv)
        STR+='Predator: Adding mutated duplicates to new pop history='+indiv.history_index+'\n'
    if Optimizer.natural_selection_scheme=='fussf':
        for ind in newpop:
            if ind.fingerprint == 0:
                ind.fingerprint = get_fingerprint(Optimizer,ind,Optimizer.fpbin,Optimizer.fpcutoff)
    if genrep >= Optimizer.reqrep*Optimizer.adaptbegin:
        ofusslim = Optimizer.fusslimit
        nfusslim = ofusslim*math.exp(-Optimizer.adaptmultiplier*float(Optimizer.genrep)/float(Optimizer.reqrep))
        Optimizer.fusslimit = nfusslim
    else:
        ofusslim = Optimizer.fusslimit
    pop = selection_switch(newpop, Optimizer.nindiv,Optimizer.natural_selection_scheme,Optimizer)
    pop = get_best(pop,len(pop))
    Optimizer.fusslimit=ofusslim
    return pop, STR
