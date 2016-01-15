import random

def swap_species(indiv, Optimizer):
    """Move function to perform atom structure swap based on atomlist
    used in strucutre with more than one atom specy.
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Swap Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')

    syms = [sym for sym,c,m,u in Optimizer.atomlist]
    numatom = [c for sym,c,m,u in Optimizer.atomlist]

    if min(numatom) > 10:
        natomsswap=random.randint(1,min(numatom)/10)
    else:
        natomsswap=1
    Optimizer.output.write('Number of swaps = '+repr(natomsswap)+'\n')

    if len(syms)==1:
        Optimizer.output.write('WARNING: Swap Mutation attempted on single atom structure system\n')
    else:
        sym1 = random.choice(syms)
        nsyms=[sym for sym in syms if sym != sym1]
        sym2 = random.choice(nsyms)

        atomlist1 = []
        atomlist2 = []
        for i in range(len(indiv[0])):
            if indiv[0][i].symbol == sym1:
               atomlist1.append(i)
            elif indiv[0][i].symbol == sym2:
               atomlist2.append(i)
        Optimizer.output.write('Swap Symbol' + repr(sym1) +  'from totally' +  repr(len(atomlist1)) +'atoms \n')
        Optimizer.output.write('with Symbol' + repr(sym2) +  'from totally' +  repr(len(atomlist2)) +'atoms \n')

        for i in range(natomsswap):
            atomindex = random.choice(atomlist1)
            a1=indiv[0][atomindex]
            a1.symbol = sym2
            atomlist1.remove(atomindex)
        for i in range(natomsswap):
            atomindex = random.choice(atomlist2)
            a2=indiv[0][atomindex]
            a2.symbol = sym1
            atomlist2.remove(atomindex)

    #Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='S'+repr(natomsswap)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
