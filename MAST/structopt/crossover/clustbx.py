import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.inp_out import write_xyz

def clustbx(ind1, ind2, Optimizer):
    """Select a box in the cluster configuration
    """
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Box Cluster Cx between individual '+repr(ind1.index)+' and individual '+repr(ind2.index)+'\n')
    
    #Perserve starting conditions of individual
    solid1=ind1[0].copy()
    solid2=ind2[0].copy()
    cello1=ind1[0].get_cell()
    cello2=ind2[0].get_cell()
    cell1=numpy.maximum.reduce(solid1.get_positions())#-numpy.minimum.reduce(solid1.get_positions())
    cell2=numpy.maximum.reduce(solid2.get_positions())#-numpy.minimum.reduce(solid2.get_positions())
    cell=numpy.minimum(cell1,cell2)
    pbc1=solid1.get_pbc()
    pbc2=solid2.get_pbc()

    #Get starting concentrations and number of atoms
    nat1=len(solid1)
    nat2=len(solid2)
    
    if Optimizer.alloy==False:
        #Atomlist-based
        concent1=[c for sym,c,m,u in Optimizer.atomlist]
        concent2=[c for sym,c,m,u in Optimizer.atomlist]
        symlist=[sym for sym,c,m,u in Optimizer.atomlist]
    else:
        symlist=list(set(solid1.get_chemical_symbols()))
        #Assumes same types of atoms in both solid1 and solid2
        concent1=[]
        for one in symlist:
            atmss=[atm for atm in solid1 if atm.symbol==one]
            concent1.append(len(atmss))
        concent2=[]
        for one in symlist:
            atmss=[atm for atm in solid2 if atm.symbol==one]
            concent2.append(len(atmss))

    # Pick a origin point for box in the cell
    #pt=[random.uniform(0,cell[0]),random.uniform(0,cell[1]),random.uniform(0,cell[2])]
    pt1=random.choice(solid1).position
    pt2=random.choice(solid2).position

    #Find max radius of circle cut
    r1=min([min(cell1[0]-pt1[0],pt1[0]),min(cell1[1]-pt1[1],pt1[1]),min(cell1[2]-pt1[2],pt1[2])])
    r2=min([min(cell2[0]-pt2[0],pt2[0]),min(cell2[1]-pt2[1],pt2[1]),min(cell2[2]-pt2[2],pt2[2])])
    r=min(r1,r2)
    mcell=min(cell)
    if r > mcell*0.125:
        r=mcell*0.125
    elif r==0.0:
        r=1.0
    if debug:
        print 'Radius of box = '+repr(r)+'\nPosition in solid1 = '+repr(pt1)+'\nPosition in solid2 = '+repr(pt2)
    
    #Find atoms within sphere of radius r
    solid1.append(Atom(position=pt1))
    dist1=[]
    blist1=[]
    for i in range(len(solid1)-1):
        d=solid1.get_distance(i,len(solid1)-1)
        if d < r:
            dist1.append((d,solid1[i]))
        else:
            blist1.append((d,solid1[i]))
    solid1.pop()

    solid2.append(Atom(position=pt2))
    dist2=[]
    blist2=[]
    for i in range(len(solid2)-1):
        d=solid2.get_distance(i,len(solid2)-1)
        if d < r:
            dist2.append((d,solid2[i]))
        else:
            blist2.append((d,solid2[i]))
    solid2.pop()
    
    #Translate spheres to opposite location
    dats1=Atoms()
    for d,atm in dist1:
        dats1.append(atm)
    dats1.translate(-pt1)
    dats1.translate(pt2)
    dats2=Atoms()
    for d,atm in dist2:
        dats2.append(atm)
    dats2.translate(-pt2)
    dats2.translate(pt1)

    #Exchange atoms in sphere and build new solids
    nsolid1=dats2.copy()
    for d,atm in blist1:
        nsolid1.append(atm)
    nsolid2=dats1.copy()
    for d,atm in blist2:
        nsolid2.append(atm)

    #Identify new number of atoms
    nnat1=len(nsolid1)
    nnat2=len(nsolid2)

    if nnat1 > nat1:
        ds=[]
        #More atoms in new solid1 means solid2 has less atoms
        #Need to transfer atoms from solid1 to solid2
        #Find atoms that are too close
        for i in range(len(dist2)):
            for j in range(len(dist2),len(nsolid1)):
                ds.append((nsolid1.get_distance(i,j),i))
        #Sort distances by longest to shortest
        ds=sorted(ds, key=lambda one:one[0])
        diff=nnat1-nat1
        exchangeindices=[]
        while len(exchangeindices) < diff:
            #Grab shortest distance in ds list
            d, index=ds.pop()
            if index not in exchangeindices:
                exchangeindices.append(index)
        nnsolid1=Atoms()
        for i in range(len(nsolid1)):
            if i not in exchangeindices:
                nnsolid1.append(nsolid1[i])
        nnsolid2=nsolid2.copy()
        for i in exchangeindices:
            #for j in range(3):
            #	position[j]=nsolid2[i].position[j]-pt1[j]+pt2[j]
            position=[random.uniform(pt2[0]-r/(2**0.5),pt2[0]+r/(2**0.5)),random.uniform(pt2[1]-r/(2**0.5),pt2[1]+r/(2**0.5)),random.uniform(pt2[2]-r/(2**0.5),pt2[2]+r/(2**0.5))]
            nnsolid2.append(Atom(symbol=nsolid1[i].symbol, position=position))
        nsolid1=nnsolid1.copy()
        nsolid2=nnsolid2.copy()

    elif nnat1 < nat1:
        ds=[]
        #More atoms in new solid2 means solid1 has less atoms
        #Need to transfer atoms from solid2 to solid1
        #Find atoms that are too close
        for i in range(len(dist1)):
            for j in range(len(dist1),len(nsolid2)):
                ds.append((nsolid2.get_distance(i,j),i))
        #Sort distances by longest to shortest
        ds=sorted(ds, key=lambda one:one[0])
        diff=nnat2-nat2
        exchangeindices=[]
        while len(exchangeindices) < diff:
            #Grab shortest distance in ds list
            d, index=ds.pop()
            if index not in exchangeindices:
                exchangeindices.append(index)
        nnsolid2=Atoms()
        for i in range(len(nsolid2)):
            if i not in exchangeindices:
                nnsolid2.append(nsolid2[i])
        nnsolid1=nsolid1.copy()
        for i in exchangeindices:
            #for j in range(3):
            #	position[j]=nsolid2[i].position[j]-pt1[j]+pt2[j]
            position=[random.uniform(pt1[0]-r/(2**0.5),pt1[0]+r/(2**0.5)),random.uniform(pt1[1]-r/(2**0.5),pt1[1]+r/(2**0.5)),random.uniform(pt1[2]-r/(2**0.5),pt1[2]+r/(2**0.5))]
            nnsolid1.append(Atom(symbol=nsolid2[i].symbol, position=position))
        nsolid1=nnsolid1.copy()
        nsolid2=nnsolid2.copy()

    #Number of atoms in each individual should now be preserved
    try:
        if Optimizer.forcing=='Concentration':
            #Need to check concentrations
            nconcent1=[]
            for one in symlist:
                atmss=[atm for atm in nsolid1 if atm.symbol==one]
                nconcent1.append(len(atmss))
            nconcent2=[]
            for one in symlist:
                atmss=[atm for atm in nsolid2 if atm.symbol==one]
                nconcent2.append(len(atmss))
            #Let's assume a random perturbation to concentration in order to correct this issue
            posd=[]
            negd=[]
            for i in range(len(nconcent1)):
                diff=nconcent1[i]-concent1[i]
                if diff >0:
                    posd.append((diff,symlist[i]))
                elif diff<0:
                    negd.append((diff,symlist[i]))
            for c,sym in negd:
                while c !=0:
                    symr=posd[0][1]
                    rlist=[atm.index for atm in nsolid1 if atm.symbol==symr]
                    index1=random.choice(rlist)
                    nsolid1[index1].symbol=sym
                    c+=1
                    posd[0]=(posd[0][0]-1,posd[0][1])
                    if posd[0][0]==0:
                        posd=posd[1::]
            #Do the same for solid2
            posd=[]
            negd=[]
            for i in range(len(nconcent2)):
                diff=nconcent2[i]-concent2[i]
                if diff >0:
                    posd.append((diff,symlist[i]))
                elif diff<0:
                    negd.append((diff,symlist[i]))
            for c,sym in negd:
                while c !=0:
                    symr=posd[0][1]
                    rlist=[atm.index for atm in nsolid2 if atm.symbol==symr]
                    index1=random.choice(rlist)
                    nsolid2[index1].symbol=sym
                    c+=1
                    posd[0]=(posd[0][0]-1,posd[0][1])
                    if posd[0][0]==0:
                        posd=posd[1::]
    except:
        f=open('problem-structures.xyz','a')
        write_xyz(f,nsolid1,data='Failed - CX(randalloybx):nsolid1 '+ind1.history_index)
        write_xyz(f,nsolid2,data='Failed - CX(randalloybx):nsolid2 '+ind2.history_index)
        nsolid1=ind1[0].copy()
        nsolid2=ind2[0].copy()
        f.close()
    #DEBUG: Write crossover to file
    if debug: 
        write_xyz(Optimizer.debugfile, nsolid1,'CX(randalloybx):nsolid1')
        write_xyz(Optimizer.debugfile, nsolid2,'CX(randalloybx):nsolid2')


    #DEBUG: Check structure of atoms exchanged
    for sym,c,m,u in Optimizer.atomlist:
            nc=len([atm for atm in nsolid1 if atm.symbol==sym])
            Optimizer.output.write('CX(clustbx):New solid1 contains '+repr(nc)+' '+repr(sym)+' atoms\n')
            nc=len([atm for atm in nsolid2 if atm.symbol==sym])
            Optimizer.output.write('CX(clustbx):New solid2 contains '+repr(nc)+' '+repr(sym)+' atoms\n')
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