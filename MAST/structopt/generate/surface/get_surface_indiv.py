from MAST.structopt.inp_out import read_xyz
from MAST.structopt.tools import find_top_layer
from MAST.structopt.switches import moves_switch
try:
    from MAST.structopt.generate.Individual import Individual
except NameError:
    print "NOTE: ASE is not installed. ASE must be installed for Structopt Individual.py to work correctly."
def get_surface_indiv(Optimizer):
    """
    Function to generate an structopt Individual class object containing as surface structure.
    Inputs:
        Optimizer = structopt Optimizer class
    Outputs:
        individ = structopt Individual class object containing surface structure data
    """
    #Load surface structure
    surfs = read_xyz(Optimizer.surfacefile)
    cells = Optimizer.surfacecell
    surfs.set_cell(cells)
    surf.set_pbc([True,True,False])
    #Find top layer
    top,bulks=find_top_layer(surfs,Optimizer.surftopthick)
    mutopto = Optimizer.mutation_options
    Optimizer.mutation_options = ['lattice_alteration_rdrd']
    topind = top.copy()
    ind = moves_switch(topind,Optimizer)
    Optimizer.mutation_options = mutopto
    individ = Individual(ind)
    individ.bulki = bulks.copy()
    individ.bulko = bulks.copy()
    return individ
