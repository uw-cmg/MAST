import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.inp_out import write_xyz

def rotct(ind1, ind2, Optimizer):
    """Rotate atoms cut and splice
    Rotates atoms randomly around center of mass and cuts with xy plane
    Maintains number of atoms
    Maintains concentration of atoms
    Returns individuals to standard positions at end (un-rotates)
    """
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Rotate Cut/Splice Cx between individual '+repr(ind1.index)+' and individual '+repr(ind2.index)+'\n')
    
    #Perserve starting conditions of individual
    indi1 = ind1[0].copy()
    indi2 = ind2[0].copy()
    
    #Translate individuals so COM is at (0,0,0)
    com1 = indi1.get_center_of_mass()
    indi1.translate(-1*com1)
    com2 = indi2.get_center_of_mass()
    indi2.translate(-1*com2)
    #Select random axis and random angle and rotate individuals
    n=0
    while n<10:
        rax = random.choice(['x', '-x','y','-y','z','-z'])
        rang = random.random()*90
        indi1.rotate(rax,a=rang,center='COM',rotate_cell=False)
        #Search for atoms in individual 1 that are above the xy plane
        group1 = Atoms(cell=ind1[0].get_cell(),pbc=ind1[0].get_pbc())
        indices1=[]
        for one in indi1:
            if one.position[2] >= 0:
                group1.append(one)
                indices1.append(one.index)
        if len(group1) > 2 and len(group1) < len(indi1):
            break
        else:
            indi1.rotate(rax,a=-1*rang,center='COM', rotate_cell=False)
        n+=1
    indi2.rotate(rax,a=rang,center='COM', rotate_cell=False)
    if debug: 
        print 'Group1 size = ', len(group1)
        print 'Position = ', [0,0,0]
        print 'Angle = ', rang
        print 'Axis = ', rax
        print 'Number of tries = ', n+1
    #Apply concentration forcing if needed
    group2 = Atoms(cell=ind2[0].get_cell(),pbc=ind2[0].get_pbc())
    indices2 = []
    dellist = []
    if Optimizer.forcing=='Concentration':
        symlist = list(set(group1.get_chemical_symbols()))
        seplist = [[atm for atm in indi2 if atm.symbol == sym] for sym, count in symlist]
        for one in group1:
            sym1=one.symbol
            listpos=[i for i,s in enumerate(symlist) if s==sym1][0]
            if len(seplist[listpos]) > 0:
                dist=[z for x,y,z in [one.position for one in seplist[listpos]]]
                pos = [i for i,value in enumerate(dist) if value==max(dist)][0]
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
        #indices2.append(pt2.index)
    else:
        for one in indi2:
            if one.position[2] >= 0:
                group2.append(one)
                indices2.append(one.index)
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
    
    other1 = Atoms()
    other2 = Atoms()
    other2 = Atoms(cell=ind2[0].get_cell(),pbc=ind2[0].get_pbc())
    for one in indi2:
        if one.index not in indices2:
            other2.append(one)
    other1 = Atoms(cell=ind1[0].get_cell(),pbc=ind1[0].get_pbc())
    for one in indi1:
        if one.index not in indices1:
            other1.append(one)
    

    indi1 = group2.copy()
    indi1.extend(other1)
    indi2 = group1.copy()
    indi2.extend(other2)
    
    #DEBUG: Write crossover to file
    if debug: 
        write_xyz(Optimizer.debugfile, group1,'group1')
        write_xyz(Optimizer.debugfile, other1,'other1')
        write_xyz(Optimizer.debugfile, group2,'group2')
        write_xyz(Optimizer.debugfile, other2,'other2')
        print 'Length of group1 = ',len(group1),'Length of group2',len(group2)
    
    
    #DEBUG: Check structure of atoms exchanged
    for sym,c,m,u in Optimizer.atomlist:
        nc=len([atm for atm in indi1 if atm.symbol==sym])
        Optimizer.output.write('CX ROTCT: Individual 1 contains '+repr(nc)+' '+repr(sym)+' atoms\n')
        nc=len([atm for atm in indi2 if atm.symbol==sym])
        Optimizer.output.write('CX ROTCT: Individual 2 contains '+repr(nc)+' '+repr(sym)+' atoms\n')

    #Need to have at least one atom of each structure in atomlist to prevent Lammps for erroring
    if Optimizer.forcing !='Concentration':
        for i in range(len(Optimizer.atomlist)):
            atms1=[inds for inds in indi1 if inds.symbol==Optimizer.atomlist[i][0]]
            atms2=[inds for inds in indi2 if inds.symbol==Optimizer.atomlist[i][0]]
            if len(atms1)==0:
                if len(atms2)==0:
                    indi1[random.randint(0,len(indi1)-1)].symbol==Optimizer.atomlist[i][0]
                    indi2[random.randint(0,len(indi2)-1)].symbol==Optimizer.atomlist[i][0]
                else:
                    indi1.append(atms2[random.randint(0,len(atms2)-1)])
                    indi1.pop(random.randint(0,len(indi1)-2))
            else:
                if len(atms2)==0:
                    indi2.append(atms1[random.randint(0,len(atms1)-1)])
                    indi2.pop(random.randint(0,len(indi2)-2))	
    indi1.rotate(rax,a=-1*rang,center='COM', rotate_cell=False)
    indi2.rotate(rax,a=-1*rang,center='COM', rotate_cell=False)
    indi1.translate(com1)
    indi2.translate(com2)

    if Optimizer.structure != 'Defect':
        cm = indi1.get_center_of_mass()
        cell = numpy.maximum.reduce(indi1.get_cell())
        cop = [cell[0]/float(2), cell[1]/float(2), cell[2]/float(2)]
        indi1.translate(-1*cm)
        indi1.translate(cop)
        cm = indi2.get_center_of_mass()
        cell = numpy.maximum.reduce(indi2.get_cell())
        cop = [cell[0]/float(2), cell[1]/float(2), cell[2]/float(2)]
        indi2.translate(-1*cm)
        indi2.translate(cop)

    ind1[0]=indi1
    ind2[0]=indi2
    
    #Check structure and number of atoms in crystal
    if Optimizer.structure=='Defect':
        solid1=Atoms()
        solid1.extend(ind1[0])
        solid1.extend(ind1.bulki)
        solid2=Atoms()
        solid2.extend(ind1[0])
        solid2.extend(ind2.bulki)
        for sym,c,m,u in Optimizer.atomlist:
                nc=len([atm for atm in solid1 if atm.symbol==sym])
                Optimizer.output.write('CX ROTCT: CIBS 1 configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n')
                nc=len([atm for atm in solid2 if atm.symbol==sym])
                Optimizer.output.write('CX ROTCT: CIBS 2 configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n')
    return ind1, ind2
