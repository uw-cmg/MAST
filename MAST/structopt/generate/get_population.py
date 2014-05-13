from MAST.structopt.generate.surface import get_surface_indiv
from MAST.structopt.generate.defect import get_defect_indiv
from MAST.structopt.generate.crystal import get_crystal_indiv
from MAST.structopt.generate import gen_pop_box, gen_pop_sphere
from MAST.structopt.generate.Individual import Individual

def get_population(Optimizer):
    """
	Function to generate a population of structures.
	Inputs:
		Optimizer = structopt Optimizer class object
	Outputs:
		pop = List of structopt Individual class objects containing new structures.
	"""
    index1 = 0
    pop = []
    for i in range(Optimizer.nindiv):
        if Optimizer.structure == 'Defect':
            individ = get_defect_indiv(Optimizer)
        elif Optimizer.structure == 'Surface':
            individ = get_surface_indiv(Optimizer)
        elif Optimizer.structure == 'Crystal':
            individ = get_crystal_indiv(Optimizer)
        else:
            if 'sphere' in Optimizer.generate_flag:
                ind = gen_pop_sphere(Optimizer.atomlist,Optimizer.size)
            else:
                ind = gen_pop_box(Optimizer.atomlist,Optimizer.size)
            individ = Individual(ind)
        individ.index = index1
        if Optimizer.genealogy: 
            individ.history_index = repr(index1)
        Optimizer.output.write('Generated cluster individual with natoms = '+
            repr(individ[0].get_number_of_atoms())+'\n')
        pop.append(individ)
        index1=index1+1
    # Generate new atomlist concentrations based on cluster+box
    if Optimizer.structure == 'Defect':
        if Optimizer.alloy:
            concents=[]
            for ind in pop:
                cs=[]
                for sym,c,u,m in Optimizer.atomlist:
                    sylen=[atm for atm in ind[0] if atm.symbol==sym]
                    cs.append(len(sylen))
                concents.append(cs)
            natmlist=[0]*len(Optimizer.atomlist)
            for i in range(len(concents[0])):
                alls=[cs[i] for cs in concents]
                avgall=int(sum(alls)/len(alls))
                if avgall>=Optimizer.atomlist[i][1]:
                    natmlist[i]=(Optimizer.atomlist[i][0], avgall,
                        Optimizer.atomlist[i][2],Optimizer.atomlist[i][3])
                else:
                    natmlist[i]=(Optimizer.atomlist[i][0], Optimizer.atomlist[i][1],
                        Optimizer.atomlist[i][2],Optimizer.atomlist[i][3])
        else:
            natmlist=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                atms1=[inds for inds in pop[0][0] if inds.symbol==Optimizer.atomlist[i][0]]
                natmlist[i]=(Optimizer.atomlist[i][0], len(atms1),Optimizer.atomlist[i][2],
                    Optimizer.atomlist[i][3])
        Optimizer.atomlist=natmlist
        Optimizer.output.write('\n\nNew atomlist concentrations based on cluster+box = '+
            repr(Optimizer.atomlist)+'\n')
    return pop
