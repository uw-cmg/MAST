import random

def swap(indiv, Optimizer):
    """Move function to perform atom structure swap based on atomlist
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
    if len(indiv[0]) > 4:
        natomsswap=random.randint(1,len(indiv[0])/5)
    else:
        natomsswap=1
    Optimizer.output.write('Number of swaps = '+repr(natomsswap)+'\n')
    syms=list(set(indiv[0].get_chemical_symbols()))
    if len(syms)<len(Optimizer.atomlist):
        syms=[sym for sym,c,m,u in Optimizer.atomlist]
    if len(syms)==1:
        Optimizer.output.write('WARNING: Swap Mutation attempted on single atom structure system\n')
    else:
        for i in range(natomsswap):
            if len(indiv[0])>1:
                a1=indiv[0][random.randint(0,indiv[0].get_number_of_atoms()-1)]
            else:
                a1=indiv[0][0]
            osym=a1.symbol
            nsymlist=[sym for sym in syms if sym != osym]
            a1.symbol=random.choice(nsymlist)
            nsym=a1.symbol
            Optimizer.output.write('Swapped '+osym+' atom with '+nsym+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='S'+repr(natomsswap)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv