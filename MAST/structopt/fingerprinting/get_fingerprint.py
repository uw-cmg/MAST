import numpy
import math
from ase import Atom, Atoms
from MAST.structopt.fingerprinting import dirac

def get_fingerprint(Optimizer,indiv,binsize,cutoffdist):
    """Function to calculate the fingerprint of a structure
    """
    rs=numpy.linspace(0.0,cutoffdist,cutoffdist/binsize)
    indi=indiv[0]
    Vuc=indi.get_volume()
    if Optimizer.structure=='Defect':
        solid=Atoms()
        solid.extend(indi)
        solid.extend(indiv.bulki)
    elif Optimizer.structure=='Crystal':
        solid=indi.repeat([3,3,3])
    else:
        solid=indi.copy()
    syms = sorted(list(set([atm.symbol for atm in solid])))
    fingerprints=[]
    for i in range(len(syms)):
        for j in range(i,len(syms)):
            indl=[atm for atm in indi if atm.symbol==syms[i]]
            ind=Atoms()
            for one in indl:
                ind.append(one)
            soll=[atm for atm in solid if atm.symbol==syms[j]]
            sol=Atoms()
            for one in soll:
                sol.append(one)
            soli=[atm for atm in solid if atm.symbol==syms[i]]
            value=[]
            for R in rs:
                value2=[]
                for k in range(len(ind)):
                    value1=[]
                    for m in range(len(sol)):
                        if k !=m:
                            rij=sol.get_distance(k,m,mic=True)
                            if rij==0: 
                                pass
                                #pdb.set_trace()
                            value1.append(dirac(R,a=rij,sig=0.02)*1/(4*math.pi*rij**2*binsize*len(soli)*len(sol)/Vuc))
                    value2.append(sum(value1))
                value.append(sum(value2))
            fingerprints.append(value)
    fpt=[]
    for one in fingerprints:
        fpt.extend(one)
    return fpt

