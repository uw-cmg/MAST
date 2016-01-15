import random
import numpy
import os
try:
    from ase import Atom, Atoms
except ImportError:
    print "NOTE: ASE is not installed. To use Structopt get_defect_indiv.py, ASE must be installed."
from MAST.structopt_stem.generate.defect import gen_solid
from MAST.structopt_stem.generate import gen_pop_box, gen_pop_sphere, generate_dumbbells, gen_pop_plate
try:
    from MAST.structopt_stem.generate.Individual import Individual
except NameError:
    print "NOTE: ASE is not installed. ASE must be installed for Structopt Individual.py to work correctly."
try:
    from mpi4py import MPI
except ImportError:
    pass

def get_defect_indiv(Optimizer):
    """
    Function to generate a structopt Individual class structure with a defect structure.
    Inputs:
    	Optimizer = structopt Optimizer class object
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
			if Optimizer.parallel:
			    from MAST.structopt_stem.tools.setup_calculator import setup_calculator
			    Optimizer.calc = setup_calculator(Optimizer)
			bulk1, PureBulkEnpa, stro = gen_solid(Optimizer.solidfile,
				Optimizer.solidcell,outfilename,Optimizer.calc,Optimizer.calc_method)
			Optimizer.output.write(stro)
		else:
			bulk1 = gen_solid(Optimizer.solidfile,Optimizer.solidcell,outfilename)
			PureBulkEnpa = 0
		natomsbulk = len(bulk1)
		Optimizer.solidbulk = bulk1.copy()
		Optimizer.purebulkenpa = PureBulkEnpa
		Optimizer.natomsbulk = natomsbulk
    # Identify nearby atoms for region 2 inclusion
    bulk = Optimizer.solidbulk.copy()
    if not Optimizer.random_loc_start:
        bulkcom = bulk.get_center_of_mass()
    else:
        cellsize = numpy.maximum.reduce(bulk.get_cell())
        dsize = max([Optimizer.size,Optimizer.sf])
        bulkcom = numpy.array([random.uniform(dsize,cellsize[one]-dsize) for one in range(3)])
    bulk.translate(-bulkcom)
    if Optimizer.sf != 0:
        bulk.append(Atom(position=[0,0,0]))
        nbulk = Atoms(pbc=True, cell=bulk.get_cell())
        nr2 = Atoms(pbc=True, cell=bulk.get_cell())
        for i in range(len(bulk)-1):
            dist = bulk.get_distance(-1,i)
            if dist <= Optimizer.sf:
                nr2.append(bulk[i])
            else:
                nbulk.append(bulk[i])
    else:
        nbulk = bulk.copy()
        nr2 = Atoms(pbc=True, cell=bulk.get_cell())
    	# Generate random individual
    if 'sphere' in Optimizer.generate_flag:
        ind = gen_pop_sphere(Optimizer.atomlist,Optimizer.size)
    elif 'dumbbell' in Optimizer.generate_flag:
        ind = Atoms(cell=[Optimizer.size for i in range(3)], pbc=True)
        for sym,c,m,u in Optimizer.atomlist:
            if c > 0:
                dums = generate_dumbbells(c, dumbbellsym=sym, nindiv=1, solid = Optimizer.solidbulk, size=Optimizer.size)[0]
                ind.extend(dums)
	## HKK: adding plate
    elif 'plate' in Optimizer.generate_flag:
        #print 'HKK:: " plate" as Generate_flag'
        #print type(gen_pop_plate)
        #print gen_pop_plate.__file__
        ind = gen_pop_plate(Optimizer.atomlist,Optimizer.size)
	## HKK: adding plate		
    else:
        #print 'HKK:: "box" as Generte flag'
        ind = gen_pop_box(Optimizer.atomlist,Optimizer.size)
    if 'random' in Optimizer.generate_flag:
        #Update atom list with atoms in region 2
        natlist = []
        for sym,c,m,u in Optimizer.atomlist:
            atsym = [atm for atm in nr2 if atm.symbol==sym]
            natlist.append((sym,len(atsym),m,u))
        nnr2 = gen_pop_sphere(natlist, Optimizer.sf*2.0)
        nnr2.translate([-Optimizer.sf,-Optimizer.sf,-Optimizer.sf])
        nnr2.set_pbc(True)
        nnr2.set_cell(bulk.get_cell())
        nr2 = nnr2.copy()
    # Initialize class individual with known values
    individ = Individual(ind)
    individ.purebulkenpa = Optimizer.purebulkenpa
    individ.natomsbulk = Optimizer.natomsbulk
    
    # Combine individual with R2
    icom = ind.get_center_of_mass()
    ind.translate(-icom)
    ind.extend(nr2)
    ind.set_pbc(True)
    ind.set_cell(bulk.get_cell())
    # Recenter structure
    nbulk.translate(bulkcom)
    ind.translate(bulkcom)
    individ[0] = ind.copy()
    individ.bulki = nbulk.copy()
    individ.bulko = nbulk.copy()
    bulk = nbulk.copy()
    bul = bulk.copy()
    for atm in individ[0]:
        bul.append(atm)
    indices = []
    for sym,c,m,u in Optimizer.atomlist:
        if c < 0:
            if Optimizer.random_vac_start:
                alist = [one for one in bul if one.symbol==sym]
                count = abs(c)
                while count > 0:
                    indices.append(random.choice(alist).index)
                    count -= 1
            else:
                pos = individ[0][0:Optimizer.natoms].get_center_of_mass()
                count = abs(c)
                bul.append(Atom(position=pos))
                alist = [one for one in bul if one.symbol==sym]
                alistd = [(bul.get_distance(len(bul)-1,one.index),one.index)
                			for one in alist]
                alistd.sort(reverse=True)
                bul.pop()
                while count > 0:
                    idx = alistd.pop()[1]
                    indices.append(idx)
                    count-=1
    if len(indices) !=0:
        nbulklist = [at for at in bul if at.index not in indices and at.index<len(bulk)]
        nalist = [at for at in bul if at.index not in indices and at.index>=len(bulk)]
        bulkn = Atoms(cell=bulk.get_cell(),pbc=True)
        for atm in nbulklist:
            bulkn.append(atm)
        individ.bulki = bulkn.copy()
        individ.bulko = bulkn.copy()
        newind = Atoms()
        for atm in nalist:
            newind.append(atm)
        newind.set_cell(individ[0].get_cell())
        newind.set_pbc(True)									
        individ[0] = newind
    return individ
