import random

def swap_int(indiv, Optimizer):
    """Move function to perform an atom structure swap based on swaplist for cluster
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
    Optimizer.output.write('Cluster Swap Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    swapmax=sum([c for sym,c in indiv.swaplist])
    if swapmax >1:
        natomsswap=random.randint(1,swapmax)
    elif swapmax==1:
        natomsswap=1
    else:
        Optimizer.output.write('Maximum number of swaps is '+repr(swapmax)+'. Unable to perform mutation')
        natomsswap=0
    if len(indiv.swaplist)==1:
        Optimizer.output.write('WARNING: Swap Mutation attempted on single atom structure system\n')
        natomsswap=0
    Optimizer.output.write('Number of swaps = '+repr(natomsswap)+'\n')
    syms=[sym for sym,c in indiv.swaplist]
    if debug: print 'Starting swaplist = '+repr(indiv.swaplist)
    
    #Extend to bulk
    if Optimizer.structure=='Defect':
        indatms=indiv[0].copy()
        brmark=len(indatms)
        indatms.extend(indiv.bulki.copy())
    else:
        indatms=indiv[0].copy()
        brmark=int(len(indatms)/2.0)
    
    #Sanity check
    sanch = [[sym,0] for sym in syms]
    for i in range(len(sanch)):
        nc = len([atm for atm in indatms if atm.symbol==sanch[i][0]])
        nc += [c for sym,c in indiv.swaplist if sym==sanch[i][0]][0]
        sanch[i][1]=nc
    
    # Swap Atoms
    for i in range(natomsswap):
        while True:
            sind=random.choice(indiv.swaplist)
            nsymlist=[sym for sym in syms if sym != sind[0]]
            natlist=[atm for atm in indatms if atm.symbol in nsymlist]
            if sind[1]>0:
                if len(natlist)>0:
                    sym1=sind[0]
                    break
        #Preferentially choose the atoms nearest the defect to swap
        pb = random.random()
        n=0
        while True:
            a1 = random.choice(natlist)
            n += 1
            if pb < 0.5:
                if n < 20:
                    if a1.index < brmark:
                        break
            else:
                if n < 20:
                    if a1.index > brmark:
                        break
            if n > 20:
                break
        nsym=indatms[a1.index].symbol
        indatms[a1.index].symbol=sym1
        Optimizer.output.write('Swapped '+nsym+' atom with '+sym1+'\n')
        for i in range(len(indiv.swaplist)):
            if indiv.swaplist[i][0]==nsym:
                indiv.swaplist[i][1]+=1
            elif indiv.swaplist[i][0]==sym1:
                indiv.swaplist[i][1]-=1
    Optimizer.output.write('New swaplist for individual: '+repr(indiv.swaplist)+'\n')
    if debug: print 'New swaplist = '+repr(indiv.swaplist)
    #Sanity check
    sanchn = [[sym,0] for sym in syms]
    for i in range(len(sanchn)):
        nc = len([atm for atm in indatms if atm.symbol==sanchn[i][0]])
        nc += [c for sym,c in indiv.swaplist if sym==sanchn[i][0]][0]
        sanchn[i][1]=nc
    
    if debug:
        for i in range(len(sanch)):
            print 'INTSWAP: Starting Atom structure '+repr(sanch[i][0])+' = '+repr(sanch[i][1])
            print 'INTSWAP: After Mutation Atom structure '+repr(sanchn[i][0])+' = '+repr(sanchn[i][1])
            
    if Optimizer.structure=='Defect':
        indiv[0]=indatms[0:brmark].copy()
        indiv.bulki=indatms[brmark::].copy()
    else:
        indiv[0]=indatms.copy()
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='SC'+repr(natomsswap)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv