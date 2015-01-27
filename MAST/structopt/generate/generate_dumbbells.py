from ase import Atom, Atoms
import random

def generate_dumbbells(ndumbbells,dumbbellsym,nindiv,solid,symbol=None,size=None):
    """Function to generate initial dumbbells 
    Inputs:
        ndumbbells = number of dumbbells to generate
        dumbbellsym = symbol of atom for dumbell
        nidiv = number of individuals to generate
        solid = ASE Atoms class object with initial bulk structure for dumbbells to reside
        symbol = dumbbell symbol matching
        size = restrict location of dumbbells
    Outputs:
        indivs = list of ASE atoms class objects with dumbbells
    *** Note: function is still in development.  Has not been implemented ***
    """
    if size==None:
        size = ndumbbells**0.333*3.0
    bulk = solid.copy()
    bulkcom = bulk.get_center_of_mass()
    bulk.translate(-bulkcom)
    bulk.append(Atom(position=[0,0,0]))
    nr2 = Atoms(pbc=True, cell=bulk.get_cell())
    for i in range(len(bulk)-1):
        dist = bulk.get_distance(-1,i)
        if dist <= size:
            nr2.append(bulk[i])
    indivs = []
    for one in range(nindiv):
        indices = []
        if symbol == None:
            opts = [atm for atm in nr2]
        else:
            opts = [atm for atm in nr2 if atm.symbol==symbol]
        indiv = Atoms()
        if len(opts) <= ndumbbells:
            ndum = len(opts)
            flag = True
        else:
            ndum = ndumbbells
            flag = False
        for one in range(ndum):
            while True:
                atm = random.choice(opts)
                if atm.index not in indices:
                    break
            datm = Atom(symbol=dumbbellsym, position=atm.position)
            indiv.append(datm)
        if flag==True:
            lack = ndumbbells-ndum
            for one in range(lack):
                indiv.append(Atom(symbol=dumbbellsym,position=random.choice(opts).position))
        indiv.translate(bulkcom)
        indiv.rattle(stdev=0.5)
        indivs.append(indiv)
    return indivs
