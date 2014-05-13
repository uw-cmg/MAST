from MAST.structopt.tools.find_defects import find_defects
from ase import Atom, Atoms

def rattle(indiv):
    """Function to slightly alter atoms in structure. Intended for use in defect function.
    """
    atms,indb,vacant,swap,stro = find_defects(indiv[0],self.solidbulk,0.0)
    atmsl,indbl,vacantl,swapl,strol = find_defects(indiv[0],self.solidbulk,2.0)
    atmslist = []
    for atm1 in atmsl:
        for atm2 in atms:
            if atm1.symbol==atm2.symbol:
                if atm1.position[0]==atm2.position[0] and atm1.position[1]==atm2.position[1] and atm1.position[2]==atm2.position[2]:
                    atmslist.append(atm1.index)
    atmolist=[atm for atm in atmsl if atm.index not in atmslist]
    rat=Atoms(cell=atms.get_cell(), pbc=atms.get_pbc)
    for one in atmolist:
        rat.append(one)
    rat.rattle(stdev=0.2)
    ind=atms.copy()
    ind.extend(rat)
    ind.extend(indbl)
    return ind
	
	