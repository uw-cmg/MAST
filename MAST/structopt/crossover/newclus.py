import random
import numpy
from ase import Atom, Atoms
from ase.calculators.neighborlist import NeighborList
from MAST.structopt.inp_out import write_xyz

def newclus(ind1, ind2, Optimizer):
    """Select a box in the cluster configuration"""
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Box Cluster Cx between individual '+repr(ind1.index)+' and individual '+repr(ind2.index)+'\n')
    
    #Perserve starting conditions of individual
    solid1 = ind1[0].copy()
    solid2 = ind2[0].copy()
    cello1 = ind1[0].get_cell()
    cello2 = ind2[0].get_cell()
    cell1 = numpy.maximum.reduce(solid1.get_positions())
    cell1m = numpy.minimum.reduce(solid1.get_positions())
    cell2 = numpy.maximum.reduce(solid2.get_positions())
    cell2m = numpy.minimum.reduce(solid2.get_positions())
    cell = numpy.minimum(cell1,cell2)
    pbc1 = solid1.get_pbc()
    pbc2 = solid2.get_pbc()
    #Get starting concentrations and number of atoms
    nat1 = len(solid1)
    nat2 = len(solid2)

    # Pick a origin point for box in the cell
    pt1 = random.choice(solid1)
    pt1f = [(pt1.position[i]-cell1m[i])/cell1[i] for i in range(3)]
    pt2 = [pt1f[i]*cell2[i]+cell2m[i] for i in range(3)]
    solid2.append(Atom(position=pt2))
    pt2=solid2[len(solid2)-1]
    #Find max neighborsize of circle cut
    r = random.uniform(0,min(nat1,nat2)/5.0)
    if debug:
        print 'DEBUG CX: Point one =', pt1.position
        print 'DEBUG CX: Point two =', pt2.position
    #Find atoms within sphere of neighborsize r for both individuals
    #Make sure that crossover is only selection of atoms not all
    while True:
        ctoff = [r for on in solid1]
        nl = NeighborList(ctoff, bothways=True, self_interaction=False)
        nl.update(solid1)
        indices1, offsets = nl.get_neighbors(pt1.index)
        if len(indices1)==0:
            r = r*1.2
        elif len(indices1) < nat1*.75:
            break
        else:
            r = r*0.8
    if debug:
        print 'Neighborsize of box = '+repr(r)+'\nPosition in solid1 = '+repr(pt1.position)+'\nPosition in solid2 = '+repr(pt2.position)
    group1 = Atoms(cell=solid1.get_cell(),pbc=solid1.get_pbc())
    group1.append(pt1)
    indices1a=[pt1.index]
    for index, d in zip(indices1,offsets):
        if index not in indices1a:
            index = int(index)
            pos = solid1[index].position + numpy.dot(d,solid1.get_cell())
            group1.append(Atom(symbol=solid1[index].symbol,position=pos))
            indices1a.append(index)
    indices1=indices1a
    ctoff = [r for on in solid2]
    nl = NeighborList(ctoff, bothways=True, self_interaction=False)
    nl.update(solid2)
    indices2, offsets = nl.get_neighbors(pt2.index)
    group2 = Atoms(cell=solid2.get_cell(),pbc=solid2.get_pbc())
    indices2a = []
    for index, d in zip(indices2,offsets):
        if index not in indices2a:
            index = int(index)
            pos = solid2[index].position + numpy.dot(d,solid2.get_cell())
            group2.append(Atom(symbol=solid2[index].symbol,position=pos))
            indices2a.append(index)
    indices2=indices2a
    if len(indices2)==0:
        for one in group1:
            while True:
                sel = random.choice(solid2)
                if sel.symbol==one.symbol:
                    if sel.index not in indices2:
                        group2.append(sel)
                        indices2.append(sel.index)
                        break

    if Optimizer.forcing=='Concentration':
        symlist = list(set(group1.get_chemical_symbols()))
        seplist = [[atm for atm in group2 if atm.symbol == sym] for sym in symlist]
        group2n=Atoms(cell=group2.get_cell(),pbc=group2.get_pbc())
        indices2n = []
        dellist = []
        for one in group1:
            sym1=one.symbol
            listpos=[i for i,s in enumerate(symlist) if s==sym1][0]
            if len(seplist[listpos]) > 0:
                pos = random.choice(range(len(seplist[listpos])))
                group2n.append(seplist[listpos][pos])
                indices2n.append(indices2[seplist[listpos][pos].index])
                del seplist[listpos][pos]
            else:
                dellist.append(one.index)
        if len(dellist) != 0:
            dellist.sort(reverse=True)
            for one in dellist:
                del group1[one]
                del indices1[one]
        indices2n.append(pt2.index)
        indices2=indices2n
        group2=group2n.copy()
    else:
        dellist = []
        while len(group2) < len(group1)-len(dellist):
            #Too many atoms in group 1
            dellist.append(random.choice(group1).index)
        if len(dellist) != 0:
            dellist.sort(reverse=True)
            for one in dellist:
                del group1[one]
                del indices1[one]
        dellist = []
        while len(group1) < len(group2)-len(dellist):
            #Too many atoms in group 2
            dellist.append(random.choice(group2).index)
        if len(dellist) != 0:
            dellist.sort(reverse=True)
            for one in dellist:
                del group2[one]
                del indices2[one]
            
    other2 = Atoms(cell=solid2.get_cell(),pbc=solid2.get_pbc())
    for one in solid2:
        if one.index not in indices2:
            other2.append(one)
    other1 = Atoms(cell=solid1.get_cell(),pbc=solid1.get_pbc())
    for one in solid1:
        if one.index not in indices1:
            other1.append(one)

    #Exchange atoms in sphere and build new solids
    nsolid1 = other1.copy()
    nsolid1.extend(group2.copy())
    nsolid2 = other2.copy()
    nsolid2.extend(group1.copy())

    #DEBUG: Write crossover to file
    if debug: 
        write_xyz(Optimizer.debugfile, nsolid1,'CX(randalloybx):nsolid1')
        write_xyz(Optimizer.debugfile, nsolid2,'CX(randalloybx):nsolid2')


    #DEBUG: Check structure of atoms exchanged
    for sym,c,m,u in Optimizer.atomlist:
        if Optimizer.structure=='Defect':
            nc=len([atm for atm in nsolid1 if atm.symbol==sym])
            nc+=len([atm for atm in ind1.bulki if atm.symbol==sym])
            oc=len([atm for atm in solid1 if atm.symbol==sym])
            oc+=len([atm for atm in ind1.bulki if atm.symbol==sym])
        else:
            nc=len([atm for atm in nsolid1 if atm.symbol==sym])
            oc=len([atm for atm in solid1 if atm.symbol==sym])
        Optimizer.output.write('CX(clustbx):New solid1 contains '+repr(nc)+' '+repr(sym)+' atoms\n')
        if debug: print 'DEBUG CX: New solid1 contains '+repr(nc)+' '+repr(sym)+' atoms'
        if oc != nc:
            #pdb.set_trace()
            print 'CX: Issue in maintaining atom concentration\n Dropping new individual'
            Optimizer.output.write('CX: Issue in maintaining atom concentration\n Dropping new individual 1\n')
            nsolid1 = solid1
    for sym,c,m,u in Optimizer.atomlist:
        if Optimizer.structure=='Defect':
            nc=len([atm for atm in nsolid2 if atm.symbol==sym])
            nc+=len([atm for atm in ind2.bulki if atm.symbol==sym])
            oc=len([atm for atm in solid2 if atm.symbol==sym])
            oc+=len([atm for atm in ind2.bulki if atm.symbol==sym])
        else:
            nc=len([atm for atm in nsolid2 if atm.symbol==sym])
            oc=len([atm for atm in solid2 if atm.symbol==sym])
        Optimizer.output.write('CX(clustbx):New solid2 contains '+repr(nc)+' '+repr(sym)+' atoms\n')
        if debug: print 'DEBUG CX: New solid2 contains '+repr(nc)+' '+repr(sym)+' atoms'
        if oc != nc:
            #pdb.set_trace()
            print 'CX: Issue in maintaining atom concentration\n Dropping new individual'
            Optimizer.output.write('CX: Issue in maintaining atom concentration\n Dropping new individual 2\n')
            solid2.pop()
            nsolid2 = solid2
    if Optimizer.forcing !='Concentration':
        for i in range(len(Optimizer.atomlist)):
            atms1=[inds for inds in nsolid1 if inds.symbol==Optimizer.atomlist[i][0]]
            atms2=[inds for inds in nsolid2 if inds.symbol==Optimizer.atomlist[i][0]]
            if len(atms1)==0:
                if len(atms2)==0:
                    nsolid1[random.randint(0,len(indi1)-1)].symbol==Optimizer.atomlist[i][0]
                    nsolid2[random.randint(0,len(indi2)-1)].symbol==Optimizer.atomlist[i][0]
                else:
                    nsolid1.append(atms2[random.randint(0,len(atms2)-1)])
                    nsolid1.pop(random.randint(0,len(nsolid1)-2))
            else:
                if len(atms2)==0:
                    nsolid2.append(atms1[random.randint(0,len(atms1)-1)])
                    nsolid2.pop(random.randint(0,len(nsolid2)-2))	

    nsolid1.set_cell(cello1)
    nsolid2.set_cell(cello2)
    nsolid1.set_pbc(pbc1)
    nsolid2.set_pbc(pbc2)
    
    ind1[0]=nsolid1.copy()
    ind2[0]=nsolid2.copy()

    return ind1, ind2
