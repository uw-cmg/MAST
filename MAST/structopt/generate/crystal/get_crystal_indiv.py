from MAST.structopt.generate import gen_pop_box, gen_pop_sphere
from MAST.structopt.generate.Individual import Individual

def get_crystal_indiv(Optimizer):
    """
    Function to generate an structopt Individual class object containing a crystal structure.
    Inputs:
    	Optimizer = structopt Optimizer class
    Outputs:
    	individ = structopt Individual class object containing crystal structure data
    """
    if 'sphere' in Optimizer.generate_flag:
        outts = gen_pop_sphere(Optimizer.atomlist,Optimizer.size,Optimizer.cell_shape_options)
    else:
        outts = gen_pop_box(Optimizer.atomlist,Optimizer.size,Optimizer.cell_shape_options)
    ind = outts[0]
    Optimizer.output.write(outts[1])
    individ = Individual(ind)
    return individ
