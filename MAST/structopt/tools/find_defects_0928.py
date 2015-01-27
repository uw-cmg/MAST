from ase import Atom, Atoms
import numpy
import copy
from ase.calculators.neighborlist import NeighborList
from MAST.structopt.inp_out.write_xyz import write_xyz
from MAST.structopt.tools.calc_dist import calc_dist

def find_defects(solid, bulko, rcutoff, atomlistcheck=False, trackvacs=False, trackswaps=False, debug=False, dcheck = 0.6):
    """Function to find interstitials, vacancies, and substitutional atoms (swaps) in a defected structure.
    Identifies species by comparison to perfect structure.
    Inputs:
        solid = ASE atoms class for defected structure
        bulko = ASE atoms class for perfect structure
        rcutoff = float value of distance to surrounding atoms to include
        atomlistcheck = False/list of atom types and concentrations according to atomlist format
        trackvacs = True/False whether or not to identify vacancies in defect
        trackswaps = True/False whether or not to identify substitutional defects
        debug = False/file object to write debug structures"""
    # Combine perfect and defect structures together
    b = bulko.copy()
    b.extend(solid)
    b.set_pbc(True)
    #Debug: Write solid and bulko to file
    if debug:
        print len(bulko)
        write_xyz(debug,b,'Find Ints: Solid and Bulko')
    # Identify nearest neighbor atoms for each atom in perfect structure
    ntot = len(bulko)
    ctoff1 = [1.2 for one in b]
    nl = NeighborList(ctoff1, bothways=True, self_interaction=False)
    nl.update(b)
    slist = []
    blist = []
    wlist = []
    #Loop over each atom in perfect structure
    for one in range(ntot):
        indices, offsets = nl.get_neighbors(one)
        for index, d in zip(indices,offsets):
            index = int(index)
            if index >= ntot:
                pos = b[index].position + numpy.dot(d,bulko.get_cell())
                natm1 = Atom(position=pos)
                dist, dx, dy, dz = calc_dist(b[one],natm1)
                if dist <= dcheck:
                    #Assume atoms closer than 0.6 Angstroms to be equivalent
                    slist.append(index-ntot)
                    blist.append(one)
                    if b[one].symbol == b[index].symbol:
                        wlist.append(index-ntot)
    #Identify those atoms corresponding to interstitials, vacancies, and substitutions
    oslist = [atm.index for atm in solid if atm.index not in slist]
    vlist = [atm.index for atm in bulko if atm.index not in blist]
    swlist = [atm.index for atm in solid if atm.index not in wlist and atm.index not in oslist]
    # Create Atoms objects for each identified defect
    ntot = len(solid)
    cluster = Atoms()
    for one in oslist:
        cluster.append(solid[one])
    vacant = Atoms()
    if trackvacs==True:
        for one in vlist:
            vacant.append(Atom(symbol=bulko[one].symbol, position=bulko[one].position))
            solid.append(Atom(symbol='X', position=bulko[one].position))
            oslist.append(len(solid)-1)
        stro = 'Cluster Identified with length = {0}\nIdentified {1} vacancies\n'.format(len(cluster),len(vlist))
    swaps = Atoms()
    if trackswaps==True:
        for one in swlist:
            swaps.append(solid[one])
            oslist.append(one)
        stro = 'Cluster Identified with length = {0}\nIdentified {1} swaps\n'.format(len(cluster),len(swlist))
    else:
        stro = 'Cluster Identified with length = {0}\n'.format(len(cluster))
    #Debug: write cluster to file
    if debug:
        b = cluster.copy()
        write_xyz(debug,b,'Find Ints: Cluster')
        debug.flush()
        print 'Found cluster size = ',len(b)
    # Identify atoms surrounding the identified defects in the defected structure
    box=Atoms()
    bulki=Atoms()
    if rcutoff != 0:
        if rcutoff > 2.0:
            cutoffs = [rcutoff for one in solid]
        else:
            cutoffs = [2.0 for one in solid]
        solid.set_pbc(True)
        nl = NeighborList(cutoffs,bothways=True,self_interaction=False)
        nl.update(solid)
        nbatmsd = []
        repinds = []
        for one in oslist:
            if one not in repinds:
                if one < ntot:
                    nbatmsd.append((0,one))
                    repinds.append(one)
            indices, offsets = nl.get_neighbors(one)
            for index,d in zip(indices,offsets):
                index = int(index)
                if index not in repinds and index < ntot:
                    opos = copy.copy(solid[index].position)
                    solid[index].position = solid[index].position + numpy.dot(d,solid.get_cell())
                    dist = solid.get_distance(one,index)
                    solid[index].position = opos
                    if dist <= rcutoff:
                        nbatmsd.append((dist,index))
                        repinds.append(index)
    else:
        nbatmsd = []
        repinds = []
        for one in oslist:
            if one not in repinds:
                if one < ntot:
                    nbatmsd.append((0,one))
                    repinds.append(one)
    nbatmsd=sorted(nbatmsd, key=lambda one:one[0], reverse=True)
    indices=[]
    i=0
    # Select only atoms closest to defects that satisfy concentrations specified by atomlist given in atomlistcheck
    if atomlistcheck:
        natomsbox=sum([c for sym,c,m,u in atomlistcheck])
        if len(nbatmsd) > natomsbox:
            while i<natomsbox:
                a=nbatmsd.pop()
                box.append(solid[a[1]])
                indices.append(a[1])
                i+=1
        elif len(nbatmsd) <= natomsbox:
            for a in nbatmsd:
                box.append(solid[a[1]])
                indices.append(a[1])
                i+=1
            if len(box) < natomsbox:
                try:
                    while True:
                        for n in range(len(nbatmsd)-1,-1,-1):
                            inds, offsets=nl.get_neighbors(nbatmsd[n][1])
                            for one,d in zip(inds,offsets):
                                if len(box) < natomsbox:
                                    if one not in indices and one < ntot:
                                        opos = copy.copy(solid[one].position)
                                        solid[one].position = solid[one].position + numpy.dot(d,solid.get_cell())
                                        dist = solid.get_distance(nbatmsd[n][1],one)
                                        solid[one].position = opos
                                        if dist <= rcutoff*5.0:
                                            box.append(solid[one])
                                            indices.append(one)
                                else:
                                    raise StopIteration()
                            for one,d in zip(inds,offsets):
                                if len(box) < natomsbox:
                                    if one not in indices and one < ntot:
                                        opos = copy.copy(solid[one].position)
                                        solid[one].position = solid[one].position + numpy.dot(d,solid.get_cell())
                                        dist = solid.get_distance(nbatmsd[n][1],one)
                                        solid[one].position = opos
                                        box.append(solid[one])
                                        indices.append(one)
                                else:
                                    raise StopIteration()
                except StopIteration:
                    pass
    # If atomlistcheck is False then use all the atoms in the given cutoff distance
    else:
        for a in nbatmsd:
            box.append(solid[a[1]])
            indices.append(a[1])
    # Add remaining atoms in defect to defected bulk atoms object
    for one in range(len(solid)):
        if one not in indices and one < ntot:
            bulki.append(solid[one])
    #Check for accidental vacancy admissions
    #checklist=[atm for atm in box if atm.symbol=='X']
    #checklist.extend([atm for atm in bulki if atm.symbol=='X'])
    #Set up new individual
    indiv=box.copy()
    bulki.set_cell(bulko.get_cell())
    indiv.set_cell(bulko.get_cell())
    bulki.set_pbc(True)
    indiv.set_pbc(True)
    stro+='New individual ({0} atoms) : {1}\n'.format(len(indiv), indiv)
    stro+='New bulki ({0} atoms) : {1}\n'.format(len(bulki), bulki)
    #Debug: write new indiv to file
    if debug:
        b=indiv.copy()
        write_xyz(debug,b,'Find Ints: New Individual')
        #Debug: write new bulki to file
        b=bulki.copy()
        write_xyz(debug,b,'Find Ints: New Bulki')
        debug.flush()
        print len(bulko)
    
    return indiv,bulki,vacant,swaps,stro