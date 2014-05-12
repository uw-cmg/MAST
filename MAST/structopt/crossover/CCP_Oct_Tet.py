import numpy
from ase import Atoms, Atom

def CCP_Oct_Tet(ind1, ind2, Optimizer):
    """CX Function underdevelopment for preferentially exchanging octahedral and tetrahedral defect sites
    Non-Functional
    """
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    uc = numpy.maximum.reduce(Optimizer.solidbulk.get_cell()/Optimizer.supercell)[0]
    interstitials = [[0.0,0.5,0.0],[0.5,0.0,0.0],[0.0,0.0,0.5],[0.5,0.5,0.5],
                    [0.25,0.25,0.25],[0.25,0.75,0.25],[0.75,0.25,0.25],[0.75,0.75,0.25],
                    [0.25,0.25,0.75],[0.25,0.75,0.75],[0.75,0.25,0.75],[0.75,0.75,0.75]]
    sites=[]
    for one in interstitials:
        pos = [uc*p for p in one]
        npos = [0,0,0]
        for i in range(Optimizer.supercell[0]):
            npos[0]=pos[0]+i*uc
            for j in range(Optimizer.supercell[1]):
                npos[1]=pos[1]+j*uc
                for k in range(Optimizer.supercell[2]):
                    npos[2]=pos[2]+k*uc
                    sites.append(copy.copy(npos))
    solidsites = Atoms()
    for one in sites:
        solidsites.append(Atom(position=one))
    return ind1, ind2
