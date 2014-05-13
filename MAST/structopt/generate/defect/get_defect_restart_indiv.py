import os
from ase import Atom, Atoms
from MAST.structopt.generate.defect import gen_solid
from MAST.structopt.tools.find_defects import find_defects
from MAST.structopt.generate.Individual import Individual
try:
	from mpi4py import MPI
except:
	pass

def get_defect_restart_indiv(Optimizer, indiv):
    """
    Function to generate an structopt Individual class object containing 
    	a defect structure from a previously existing structure
    Inputs:
    	Optimizer = structopt Optimizer class
    	indiv = ASE Atoms object containing the previously existing structure
    Outputs:
    	individ = structopt Individual class object containing defect structure data
    """
    if not Optimizer.solidbulk:
		#Initialize Bulk - Generate or load positions of bulk solid
		try:
		    rank = MPI.COMM_WORLD.Get_rank()
		except:
		    rank = 0
		outfilename = os.path.join(os.path.join(os.getcwd(),Optimizer.filename+'-rank'+repr(rank)),'Bulkfile.xyz')
		if Optimizer.evalsolid:
			bulk1, PureBulkEnpa, stro = gen_solid(Optimizer.solidfile,
				Optimizer.solidcell,outfilename,Optimizer.calc,Optimizer.calc_method)
			Optimizer.output.write(stro)
		else:
			bulk1 = gen_solid(Optimizer.solidfile,Optimizer.solidcell,outfilename)
			PureBulkEnpa = 0
		natomsbulk = len(bulk1)
		Optimizer.solidbulk = bulk1.copy()
		Optimizer.summary.write('CIBS Run Pure Bulk Energy per Atom:'+
			repr(PureBulkEnpa)+'\n')
		Optimizer.purebulkenpa = PureBulkEnpa
		Optimizer.natomsbulk = natomsbulk
    indiv.set_cell(Optimizer.solidcell)
    indiv.set_pbc(True)
    if Optimizer.restart_ints == 0:
        outt = find_defects(indiv,Optimizer.solidbulk,Optimizer.sf)
    else:
        indicop = [atm for atm in indiv if atm.symbol != 'X']
        indiv = Atoms(cell=Optimizer.solidcell, pbc=True)
        for atm in indicop:
            indiv.append(atm)
        outt=[indiv[0:Optimizer.restart_ints],indiv[Optimizer.restart_ints::], Atoms(),
        	Atoms(),'Assuming first '+repr(Optimizer.restart_ints)+' are interstitials\n']
    indi = outt[0].copy()
    bulki = outt[1].copy()
    individ = Individual(indi)
    individ.bulko = bulki.copy()
    individ.bulki = bulki.copy()
    individ.purebulkenpa = Optimizer.purebulkenpa
    individ.natomsbulk = Optimizer.natomsbulk
    individ.vacancies = outt[2].copy()
    individ.swaps = outt[3].copy()
    Optimizer.output.write(outt[4])
    return individ
