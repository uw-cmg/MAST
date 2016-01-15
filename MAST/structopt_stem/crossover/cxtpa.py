import random
from ase import Atom, Atoms

def cxtpa(ind1, ind2, Optimizer):
    """ Simple Two point Crossover of atoms
    Maintains concentration
    Maintains total number of atoms
    Random clusters not necessarily geometry based
    """
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    # Get individual atoms
    Optimizer.output.write('Two Pt Atom Cx between individual '+repr(ind1.index)+' and individual '+repr(ind2.index)+'\n')
    indi1=ind1[0]
    indi2=ind2[0]
    #Identify range of atoms in individual 1 to exchange
    cxpoint1 = random.randint(0, len(indi1)-1)
    cxpoint2 = random.randint(1, len(indi1)-1)
    if cxpoint2 == cxpoint1:
        cxpoint2 += 1
    elif cxpoint2 < cxpoint1:
        # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1
    group1 = Atoms()
    symlist = []
    indices1 = []
    for one in range(cxpoint1, cxpoint2):
        group1.append(indi1[one])
        indices1.append(one)
        if indi1[one].symbol not in symlist:
            symlist.append(indi1[one].symbol)
    #Identify atoms in individual 2 of similar structure as in individual 1
    group2 = Atoms()
    seplist = [[atm for atm in indi2 if atm.symbol == sym] for sym in symlist]
    indices2 = []
    dellist = []
    for one in group1:
        sym1=one.symbol
        listpos=[i for i,value in enumerate(symlist) if value==sym1][0]
        if len(seplist[listpos]) > 0:
            pos = random.choice(range(len(seplist[listpos])))
            group2.append(seplist[listpos][pos])
            indices2.append(seplist[listpos][pos].index)
            del seplist[listpos][pos]
        else:
            dellist.append(one.index)
    if len(dellist) != 0:
        dellist.sort(reverse=True)
        for one in dellist:
            del group1[one]
            del indices1[one]
    #Identify atoms not to be exchanged
    other1=Atoms()
    for one in indi1:
        if one.index not in indices1:
            other1.append(one)
    other2 = Atoms()
    for one in indi2:
        if one.index not in indices2:
            other2.append(one)
    #Perform exchange
    indi1n = other1.copy()
    indi1n.extend(group2)
    indi2n = other2.copy()
    indi2n.extend(group1)
    #Perserve cell and boundary conditions
    indi1n.set_cell(indi1.get_cell())
    indi1n.set_pbc(indi1.get_pbc())
    indi2n.set_cell(indi2.get_cell())
    indi2n.set_pbc(indi2.get_pbc())
    ind1[0]=indi1n
    ind2[0]=indi2n
    #Debug options
    if debug:
        print 'DEBUG CX: Number of atoms exchanged from Individual 1 = ', len(group1)
        print 'DEBUG CX: Number of atoms exchanged from Individual 2 = ', len(group2)
        for s,c,m,u in Optimizer.atomlist:
            print 'DEBUG CX: Individual 1 has ', len([atm for atm in indi1n if atm.symbol==s]), ' atoms of structure ', s
            print 'DEBUG CX: Individual 2 has ', len([atm for atm in indi2n if atm.symbol==s]), ' atoms of structure ', s
        print 'DEBUG CX: Number of atoms in Individual 1 = ', len(indi1n)
        print 'DEBUG CX: Number of atoms in Individual 2 = ', len(indi2n)
    return ind1, ind2
