import random
import numpy
from ase import Atom, Atoms

def read_Ltraj_file(filename,atomlist=False,timestep=-1,atomlist=False,ratomlist=False,writefile=False):
    """Function to convert LAMMPS Trajectory file to ASE atoms object
    Inputs:
        filename = string of LAMMPS Trajectory file to read
        timestep = Which timestep to read from the trajectory file.
            Default is last time step available
            If timestep is set to all then outputs become lists
        atomlist = list of element strings corresponding to number type
        ratomlist = Boolean for whether or not to return the atom list
        write_file = Boolean for whether or not to write ASE Atoms object to file
    Ouput:
        ASE atoms object containing structure from timestep specified
        if ratomlist: list containing element strings applied to structure corresponding to LAMMPS number
        vlist = list of velocities for each atom in structure at timestep
        flist = list of forcer for each atom in structure at timestep
    """
    f = open(filename,'r')
    if atomlist==False:
        atnum = []
        attype = []
    alist = []
    vlist = []
    flist = []
    timeflag = False
    numberflag = False
    boxflag = False
    atomsflag = False
    for line in f.readlines():
        sp = line.split()
        if timeflag == True:
            timeflag = False
        elif numberflag == True:
            noa = float(sp[0])
            numberflag = False
        elif boxflag == True:
            if boxcount <= 2:
                box.append(float(sp[1]))
                boxcount += 1
            if boxcount > 2:
                boxflag = False
        elif atomsflag == True:
            if na < noa:
                if atomlist:
                    sym = atomlist[[i for i,sym in enumerate(atomlist) if str(i+1)==sp[1]][0]]
                    at = Atom(symbol=sym, position=[float(sp[2]),float(sp[3]),float(sp[4])])
                    a.append(at)
                    v = [float(sp[i]) for i in range(5,8)]
                    vs.append(v)
                    f1 = [float(sp[i]) for i in range(8,11)]
                    fs.append(f1)
                    na += 1
                else:
                    if sp[1] in atnum:
                        sym = attype[[i for i,num in enumerate(atnum) if num==sp[1]][0]]
                        at=Atom(symbol=sym,position=[float(sp[2]),float(sp[3]),float(sp[4])])
                        a.append(at)
                    else:
                        symn = random.choice(range(1,100))
                        at=Atom(symbol=symn,position=[float(sp[2]),float(sp[3]),float(sp[4])])
                        a.append(at)
                        atnum.append(sp[1])
                        attype.append(at.symbol)
                    v = [float(sp[i]) for i in range(5,8)]
                    vs.append(v)
                    f1 = [float(sp[i]) for i in range(8,11)]
                    fs.append(f1)
                    na += 1
            if na >= noa:
                atomsflag = False
                if atomlist==False:
                    atomlist=list(numpy.zeros(len(attype)))
                    for i in range(len(atnum)):
                        atomlist[int(atnum[i])-1]=attype[i]
                alist.append(a)
                vlist.append(vs)
                flist.append(fs)
        elif sp[1] == 'TIMESTEP':
            timeflag = True
        elif sp[1] == 'NUMBER':
            numberflag = True
        elif sp[1] == 'BOX':
            boxflag = True
            box = []
            boxcount = 0
            pbcc = [False,False,False]
            if sp[3] == 'pp':
                pbcc[0] = True
            if sp[4] == 'pp':
                pbcc[1] = True
            if sp[5] == 'pp':
                pbcc[2] = True
        elif sp[1] == 'ATOMS':
            atomsflag = True
            na = 0
            a = Atoms(cell=box,pbc=pbcc)
            vs = []
            fs = []
    f.close()
    if timestep == 'all':
        writefile=False
    if writefile == True:
        write_xyz(filename+'.xyz',alist[timestep],'xyz')
    if timestep == 'all':
        if ratomlist:
            return alist, atomlist, vlist, flist
        else:
            return alist, vlist, flist
    else:
        if ratomlist:
            return alist[timestep], atomlist, vlist[timestep], flist[timestep]
        else:
            return alist[timestep], vlist[timestep], flist[timestep]

