import os
import math
import numpy
from MAST.structopt.inp_out import read_xyz
from MAST.structopt.tools import calc_dist

def get_lattice_concentration(bulkfile,indivfile):
    """Function to identify the lattice concentration of atoms in a bulk structure compared to a 
    structure with a defect.
    Inputs:
        bulkfile = filename for starting structure with original lattice atoms
        indivfile = filename for structure to compare
    Outputs:
        File: LatticeConcentration.txt in working directory.  Includes summary of concentration
        of each atom type and vacancies
    ** Note: Currently limited to cubic structures **
    """
    # Load Bulk Solid File
    solid=read_xyz(bulkfile)

    # Get lattice sites for bulk

    #Identify nearest neighbor distance
    solid.set_pbc(True)
    distmin=[]
    for i in range(20):
        dist=[]
        for j in range(len(solid)):
            if i !=j:
                d = calc_dist(solid[i],solid[j])
                dist.append(d)
        distmin.append(min(dist))
    nndist=sum([one for one,x,y,z in distmin])/len(distmin)
    nnxd=sum([x for one,x,y,z in distmin])/len(distmin)
    nnyd=sum([y for one,x,y,z in distmin])/len(distmin)
    nnzd=sum([z for one,x,y,z in distmin])/len(distmin)

    solid.translate([nnxd/2.0,nnyd/2.0,nnzd/2.0])

    # Get size of cell for pbc
    cell = numpy.maximum.reduce(solid.get_positions())
    cell += [nnxd/2.0,nnyd/2.0,nnzd/2.0]

    # Initialize boxes
    nx=int(math.ceil(float(cell[0])/float(nnxd)))
    ny=int(math.ceil(float(cell[1])/float(nnyd)))
    nz=int(math.ceil(float(cell[2])/float(nnzd)))
    bxarray0=[[0,[]] for i in range(nx*ny*nz)]

    # Identify which atoms are in which box
    positions=solid.get_positions()
    for i in range(len(solid)):
        box=[math.floor(positions[i][0]/nnxd),math.floor(positions[i][1]/nnyd),math.floor(positions[i][2]/nnzd)]
        bxarray0[int((nx*ny)*(box[2])+nx*(box[1])+box[0])][0]+=1
        bxarray0[int((nx*ny)*(box[2])+nx*(box[1])+box[0])][1]+=[solid[i].symbol]

    #Get types of atoms for bulk and bulk lattice concentration
    nlatsites=len(solid)
    syms =list(set([atm.symbol for atm in solid]))
    nsyms = []
    for one in syms:
        numberofsym=len([atm for atm in solid if atm.symbol==one])
        nsyms.append(float(numberofsym)/float(nlatsites))
    concentbulk=zip(syms,nsyms)

    #Get Lattice sites for individual
    onlatcon=[[-1,concentbulk]]
    n=0
    while True:
        try:
            indiv=read_xyz(indivfile,n)
        except:
            break
        indiv.translate([nnxd/2.0,nnyd/2.0,nnzd/2.0])
        bxarray=[[0,[],[]] for i in range(nx*ny*nz)]
        positions=indiv.get_positions()
    
        # Wrap positions in individual to cell size
        for i in range(len(positions)):
            while positions[i][0] > cell[0]:
                positions[i][0]=positions[i][0]-cell[0]
            while positions[i][1] > cell[1]:
                positions[i][1]=positions[i][1]-cell[1]
            while positions[i][2] > cell[2]:
                positions[i][2]=positions[i][2]-cell[2]
            while positions[i][0] < 0:
                positions[i][0]=positions[i][0]+cell[0]
            while positions[i][1] < 0:
                positions[i][1]=positions[i][1]+cell[1]
            while positions[i][2] < 0:
                positions[i][2]=positions[i][2]+cell[2]
    
        for i in range(len(indiv)):
            box=[math.floor(positions[i][0]/nnxd),math.floor(positions[i][1]/nnyd),math.floor(positions[i][2]/nnzd)]
            bxarray[int((nx*ny)*(box[2])+nx*(box[1])+box[0])][0]+=1
            bxarray[int((nx*ny)*(box[2])+nx*(box[1])+box[0])][1]+=[indiv[i].symbol]
            bxarray[int((nx*ny)*(box[2])+nx*(box[1])+box[0])][2]+=[i]

        # Get on-lattice concentration
        latsyms=[]
        for i in range(len(bxarray0)):
            if bxarray0[i][0]!=0:
                if len(bxarray[i][1])==1:
                    latsyms.extend(bxarray[i][1])
                elif len(bxarray[i][1]) >1:
                    if bxarray[i][1][0] in bxarray[i][1]:
                        latsyms.extend(bxarray0[i][1])
                else:
                    latsyms.append('Vacancy')
        reducedsyms=list(set(latsyms))
        concenti=[]
        for one	in reducedsyms:
            numberofsym=len([atm for atm in latsyms if atm==one])
            concenti.append(float(numberofsym)/float(nlatsites))
        concents=zip(reducedsyms,concenti)
        onlatcon.append([n,concents])
        n+=1

    maxlen=max([len(con) for n,con in onlatcon])
    maxsyms=[con for n,con in onlatcon if len(con)==maxlen]
    symlist=[sym for sym,n in maxsyms[0]]
    output=open('LatticeConcentration.txt','a')
    output.write('Generation ')
    for one in symlist:
        output.write(repr(one)+' ')
    output.write('TotalSites \n')
    for n,con in onlatcon:
        output.write(repr(n)+' ')
        for sym in symlist:
            num=[count for atm,count in con if atm==sym]
            num=sum(num)
            output.write(repr(num)+' ')
        output.write(repr(nlatsites)+'\n')
    output.close()

