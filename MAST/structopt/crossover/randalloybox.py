import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.tools.find_defects import find_defects
from MAST.structopt.inp_out import write_xyz

def randalloybox(ind1, ind2, Optimizer):
    """Select a box in the alloy configuration
    *** Note: CX may be obsolete. In development ***
    """
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    if Optimizer.structure != 'Defect':
        Optimizer.output.write('WARNING: Box Random alloy Cx attempted on nondefect structure. SKIPPING CX.\n')
    else:
        Optimizer.output.write('Box Random alloy Cx between individual '+repr(ind1.index)+' and individual '+repr(ind2.index)+'\n')
    
        #Perserve starting conditions of individual
        indi1 = ind1[0].copy()
        indi2 = ind2[0].copy()
        cell1 = numpy.maximum.reduce(indi1.get_cell())
        cell2 = numpy.maximum.reduce(indi2.get_cell())
        cell = numpy.minimum(cell1,cell2)
        pbc1 = indi1.get_pbc()
        pbc2 = indi2.get_pbc()

        #Build solids
        solid1 = indi1.copy()
        solid1.extend(ind1.bulki.copy())
        solid2 = indi2.copy()
        solid2.extend(ind2.bulki.copy())

        #Get starting concentrations and number of atoms
        nat1 = len(solid1)
        nat2 = len(solid2)
        symlist = list(set(solid1.get_chemical_symbols()))
        #Assumes same types of atoms in both solid1 and solid2
        concent1 = []
        concent2 = []
        for one in symlist:
            atmss = [atm for atm in solid1 if atm.symbol==one]
            concent1.append(len(atmss))
            atmss = [atm for atm in solid2 if atm.symbol==one]
            concent2.append(len(atmss))

        # Pick a origin point for box in the cell
        #pt=[random.uniform(0,cell[0]),random.uniform(0,cell[1]),random.uniform(0,cell[2])]
        pt1 = solid1[0].position
        pt2 = solid2[0].position

        #Find max radius of circle cut
        r1 = min([min(cell1[0]-pt1[0],pt1[0]),min(cell1[1]-pt1[1],pt1[1]),min(cell1[2]-pt1[2],pt1[2])])
        r2 = min([min(cell2[0]-pt2[0],pt2[0]),min(cell2[1]-pt2[1],pt2[1]),min(cell2[2]-pt2[2],pt2[2])])
        r = min(r1,r2)
        mcell = min(cell)
        if r > mcell*0.125:
            r = mcell*0.125
        if debug:
            Optimizer.output.write('Radius of box = '+repr(r)+'\nPosition in solid1 = '+repr(pt1)+'\nPosition in solid2 = '+repr(pt2)+'\n')

        #Find atoms within sphere of radius r
        solid1.append(Atom(position=pt1))
        dist1 = []
        blist1 = []
        for i in range(len(solid1)-1):
            d = solid1.get_distance(i,len(solid1)-1)
            if d < r:
                dist1.append((d,solid1[i]))
            else:
                blist1.append((d,solid1[i]))
        solid1.pop()

        solid2.append(Atom(position=pt2))
        dist2 = []
        blist2 = []
        for i in range(len(solid2)-1):
            d = solid2.get_distance(i,len(solid2)-1)
            if d < r:
                dist2.append((d,solid2[i]))
            else:
                blist2.append((d,solid2[i]))
        solid2.pop()

        #Translate spheres to opposite location
        dats1 = Atoms()
        for d,atm in dist1:
            dats1.append(atm)
        dats1.translate(-pt1)
        dats1.translate(pt2)
        dats2 = Atoms()
        for d,atm in dist2:
            dats2.append(atm)
        dats2.translate(-pt2)
        dats2.translate(pt1)

        #Exchange atoms in sphere and build new solids
        nsolid1 = dats2.copy()
        for d,atm in blist1:
            nsolid1.append(atm)
        nsolid2 = dats1.copy()
        for d,atm in blist2:
            nsolid2.append(atm)

        #Identify new number of atoms
        nnat1 = len(nsolid1)
        nnat2 = len(nsolid2)

        if nnat1 > nat1:
            ds = []
            #More atoms in new solid1 means solid2 has less atoms
            #Need to transfer atoms from solid1 to solid2
            #Find atoms that are too close
            for i in range(len(dist2)):
                for j in range(len(dist2),len(nsolid1)):
                    ds.append((nsolid1.get_distance(i,j),i))
            #Sort distances by longest to shortest
            ds = sorted(ds, key=lambda one:one[0])
            diff = nnat1-nat1
            exchangeindices = []
            while len(exchangeindices) < diff:
                #Grab shortest distance in ds list
                d, index = ds.pop()
                if index not in exchangeindices:
                    exchangeindices.append(index)
            nnsolid1 = Atoms()
            for i in range(len(nsolid1)):
                if i not in exchangeindices:
                    nnsolid1.append(nsolid1[i])
            nnsolid2 = nsolid2.copy()
            for i in exchangeindices:
                #for j in range(3):
                #	position[j]=nsolid2[i].position[j]-pt1[j]+pt2[j]
                position = [random.uniform(pt2[0]-r/(2**0.5),pt2[0]+r/(2**0.5)),random.uniform(pt2[1]-r/(2**0.5),pt2[1]+r/(2**0.5)),random.uniform(pt2[2]-r/(2**0.5),pt2[2]+r/(2**0.5))]
                nnsolid2.append(Atom(symbol=nsolid1[i].symbol, position=position))
            nsolid1 = nnsolid1.copy()
            nsolid2 = nnsolid2.copy()

        elif nnat1 < nat1:
            ds = []
            #More atoms in new solid2 means solid1 has less atoms
            #Need to transfer atoms from solid2 to solid1
            #Find atoms that are too close
            for i in range(len(dist1)):
                for j in range(len(dist1),len(nsolid2)):
                    ds.append((nsolid2.get_distance(i,j),i))
            #Sort distances by longest to shortest
            ds = sorted(ds, key=lambda one:one[0])
            diff = nnat2-nat2
            exchangeindices = []
            while len(exchangeindices) < diff:
                #Grab shortest distance in ds list
                d, index = ds.pop()
                if index not in exchangeindices:
                    exchangeindices.append(index)
            nnsolid2 = Atoms()
            for i in range(len(nsolid2)):
                if i not in exchangeindices:
                    nnsolid2.append(nsolid2[i])
            nnsolid1 = nsolid1.copy()
            for i in exchangeindices:
                #for j in range(3):
                #	position[j]=nsolid2[i].position[j]-pt1[j]+pt2[j]
                position = [random.uniform(pt1[0]-r/(2**0.5),pt1[0]+r/(2**0.5)),random.uniform(pt1[1]-r/(2**0.5),pt1[1]+r/(2**0.5)),random.uniform(pt1[2]-r/(2**0.5),pt1[2]+r/(2**0.5))]
                nnsolid1.append(Atom(symbol=nsolid2[i].symbol, position=position))
            nsolid1 = nnsolid1.copy()
            nsolid2 = nnsolid2.copy()

        #Number of atoms in each individual should now be preserved
        if Optimizer.forcing=='Concentration':
            #Need to check concentrations
            nconcent1 = []
            for one in symlist:
                atmss = [atm for atm in nsolid1 if atm.symbol==one]
                nconcent1.append(len(atmss))
            nconcent2 = []
            for one in symlist:
                atmss = [atm for atm in nsolid2 if atm.symbol==one]
                nconcent2.append(len(atmss))
            #Let's assume a random perturbation to concentration in order to correct this issue
            posd = []
            negd = []
            for i in range(len(nconcent1)):
                diff = nconcent1[i]-concent1[i]
                if diff >0:
                    posd.append((diff,symlist[i]))
                elif diff<0:
                    negd.append((diff,symlist[i]))
            for c,sym in negd:
                while c !=0:
                    symr = posd[0][1]
                    rlist = [atm.index for atm in nsolid1 if atm.symbol==symr]
                    index1 = random.choice(rlist)
                    nsolid1[index1].symbol=sym
                    c += 1
                    posd[0] = (posd[0][0]-1,posd[0][1])
                    if posd[0][0]==0:
                        posd = posd[1::]
            #Do the same for solid2
            posd = []
            negd = []
            for i in range(len(nconcent2)):
                diff = nconcent2[i]-concent2[i]
                if diff >0:
                    posd.append((diff,symlist[i]))
                elif diff<0:
                    negd.append((diff,symlist[i]))
            for c,sym in negd:
                while c !=0:
                    symr = posd[0][1]
                    rlist = [atm.index for atm in nsolid2 if atm.symbol==symr]
                    index1 = random.choice(rlist)
                    nsolid2[index1].symbol=sym
                    c += 1
                    posd[0] = (posd[0][0]-1,posd[0][1])
                    if posd[0][0]==0:
                        posd = posd[1::]

        #DEBUG: Write crossover to file
        if debug: 
            write_xyz(Optimizer.debugfile, nsolid1,'CX(randalloybx):nsolid1')
            write_xyz(Optimizer.debugfile, nsolid2,'CX(randalloybx):nsolid2')

        #DEBUG: Check structure of atoms exchanged
        for sym,c,m,u in Optimizer.atomlist:
                nc = len([atm for atm in nsolid1 if atm.symbol==sym])
                Optimizer.output.write('CX(randalloybx):New solid1 contains '+repr(nc)+' '+repr(sym)+' atoms\n')
                nc = len([atm for atm in nsolid2 if atm.symbol==sym])
                Optimizer.output.write('CX(randalloybx):New solid2 contains '+repr(nc)+' '+repr(sym)+' atoms\n')
        if Optimizer.forcing !='Concentration':
            for i in range(len(Optimizer.atomlist)):
                atms1 = [inds for inds in nsolid1 if inds.symbol==Optimizer.atomlist[i][0]]
                atms2 = [inds for inds in nsolid2 if inds.symbol==Optimizer.atomlist[i][0]]
                if len(atms1)==0:
                    if len(atms2)==0:
                        nsolid1[random.randint(0,len(indi1)-1)].symbol = Optimizer.atomlist[i][0]
                        nsolid2[random.randint(0,len(indi2)-1)].symbol = Optimizer.atomlist[i][0]
                    else:
                        nsolid1.append(atms2[random.randint(0,len(atms2)-1)])
                        nsolid1.pop(random.randint(0,len(nsolid1)-2))
                else:
                    if len(atms2)==0:
                        nsolid2.append(atms1[random.randint(0,len(atms1)-1)])
                        nsolid2.pop(random.randint(0,len(nsolid2)-2))	

        nsolid1.set_cell(cell1)
        nsolid2.set_cell(cell2)
        nsolid1.set_pbc(pbc1)
        nsolid2.set_pbc(pbc2)

        outs = find_defects(nsolid1,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=False)
        ind1[0] = outs[0].copy()
        ind1.bulki = outs[1].copy()
        ind1.vacancies = outs[2].copy()
        ind1.swaps = outs[3].copy()
        outs = find_defects(nsolid2,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=False)
        ind2[0] = outs[0].copy()
        ind2.bulki = outs[1].copy()
        ind2.vacancies = outs[2].copy()
        ind2.swaps = outs[3].copy()
    
    return ind1, ind2
