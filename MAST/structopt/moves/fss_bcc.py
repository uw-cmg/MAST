from ase.calculators.neighborlist import NeighborList
from scipy.stats import mode
from ase import Atoms, Atom
import copy
import random
import numpy
import math
from MAST.structopt.tools.calc_dist import calc_dist
from MAST.structopt.tools.eval_energy import run_energy_eval as REE

def fss_bcc(indiv, Optimizer):
    defected = indiv[0].copy()
    defected.extend(indiv.bulki)
    #Identify nearest neighbor cutoff distance
    nndists = []
    for i in range(5):
        one = random.choice(defected)
        distances = [defected.get_distance(one.index, j) for j in range(len(defected)) if j != one.index]
        distances.sort()
        nndists.extend(distances[0:3])
    nndist = sum(nndists)/len(nndists)
    cutoff = nndist*0.6
    #Create nearest neighbor list from cutoff distance
    ctflist = [cutoff for one in defected]
    nl = NeighborList(ctflist, bothways=True, self_interaction=False)
    nl.update(defected)
    #Identify most common number of nearest neighbors for each atom
    nneigh = []
    for one in defected:
        indices, offsets = nl.get_neighbors(one.index)
        nneigh.append(len(indices))
    avn = mode(nneigh)
    #Identify those atoms that have a different number of nearest neighbors
    defs = [i for i in range(len(nneigh)) if nneigh[i] != avn[0][0]]
    #Create new structure from translated defected atoms
    defsat = Atoms(pbc = defected.get_pbc(), cell=defected.get_cell())
    for i in defs:
        defsat.append(defected[i])
    #Identify center of mass of defected group and translate to center of cell
    cop = position_average(defsat)
    ndefsat = shift_atoms(defsat, cop)
    ndefected = shift_atoms(defected, cop)
    #Identify bounds of defected portion of structure
    maxpos = max(numpy.maximum.reduce(ndefsat.positions))
    minpos = min(numpy.minimum.reduce(ndefsat.positions))
    #Identify size of structure that will encompass defected structure
    osc = copy.deepcopy(Optimizer.supercell)
    latcont = numpy.maximum.reduce(ndefsat.get_cell())[0]/osc[0]
    deltapos = abs(maxpos-minpos)
    newsupercell = round(deltapos/latcont)+1
    bxlen = newsupercell*latcont
    #Identify those atoms within a box of encompassing length centered at the center of the cell
    tol = .1
    cell = numpy.maximum.reduce(ndefsat.get_cell())
    boxlow = [cell[i]/2.0-bxlen/2.0-tol for i in range(3)]
    boxhigh = [cell[i]/2.0+bxlen/2.0+tol for i in range(3)]
    atlist = []
    otherlist = []
    for one in ndefected:
        pos = one.position
        if boxlow[0] < pos[0] < boxhigh[0]:
            if boxlow[1] < pos[1] < boxhigh[1]:
                if boxlow[2] < pos[2] < boxhigh[2]:
                    atlist.append(one)
                else:
                    otherlist.append(one)
            else:
                otherlist.append(one)
        else:
            otherlist.append(one)
    ncell = [bxlen,bxlen,bxlen]
    #Create a new atoms object from the atoms within the box
    ssats = Atoms(pbc = True, cell=ncell)
    for one in atlist:
        ssats.append(one)
    #Attach a calculator for the atoms
    ssats.set_calculator(Optimizer.calc)
    #Calculate the energy
    out = REE(ssats)

totalsol, energy, pressure, volume, STR

    return defsat  

def position_average(atoms, mode='median', trimp = 0.05):
    """Function to calculate the trimmed mean or median
    of the positions of atoms
    Input:
        atoms = ASE Atoms object
        mode = method to use either median or trimmedmean
            default value is median
        trimp = percentage of high and lows to leave out
            default value is 5% == 0.05
    Output:
        [x,y,z] = position of atoms center
    """
    positions = atoms.get_positions()
    xs = []
    ys = []
    zs = []
    for x,y,z in positions:
        xs.append(x)
        ys.append(y)
        zs.append(z)
    xs.sort()
    ys.sort()
    zs.sort()
    if mode=='median':
        xc = xs[len(xs)/2]
        yc = ys[len(ys)/2]
        zc = zs[len(zs)/2]
    else:
        ndel = int(round(trimp*len(xs)))
        xc = sum(xs[ndel:len(xs)-ndel])/(len(xs)-2*ndel)
        yc = sum(ys[ndel:len(ys)-ndel])/(len(ys)-2*ndel)
        zc = sum(zs[ndel:len(zs)-ndel])/(len(zs)-2*ndel)
    return [xc,yc,zc]

def shift_atoms(atoms, position=None):
    """Function to shift the positions of the atoms so that
    the given location is in the center of the cell.
    Input:
        atoms = ASE atoms object to be shifted
        position = [x,y,z] position to shift to center of cell
            default is center of mass of atoms
    Output:
        atoms = ASE atoms object that has been shifted
    """
    atoms = atoms.copy()
    if not position:
        position = atoms.get_center_of_mass
    cell = numpy.maximum.reduce(atoms.get_cell())
    trans = [cell[i]/2.0-position[i] for i in range(3)]
    atoms.translate(trans)
    positions = atoms.get_positions()
    for i in range(len(positions)):
        if atoms.pbc[0]:
            if trans[0] > 0:
                if positions[i][0] > cell[0]:
                    positions[i][0] -= cell[0]
            else:
                if positions[i][0] < 0:
                    positions[i][0] += cell[0]
        if atoms.pbc[1]:
            if trans[1] > 0:
                if positions[i][1] > cell[1]:
                    positions[i][1] -= cell[1]
            else:
                if positions[i][1] < 0:
                    positions[i][1] += cell[1]
        if atoms.pbc[2]:
            if trans[2] >0:
                if positions[i][2] > cell[2]:
                    positions[i][2] -= cell[2]
            else:
                if positions[i][2] < 0:
                    positions[i][2] += cell[2]
    atoms.set_positions(positions)
    return atoms

