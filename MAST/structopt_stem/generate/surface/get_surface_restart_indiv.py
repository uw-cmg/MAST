from MAST.structopt_stem.inp_out import read_xyz
try:
    from MAST.structopt_stem.generate.Individual import Individual
except NameError:
    print "ASE is not installed. ASE must be installed for Structopt Individual.py to work correctly."
from MAST.structopt_stem.tools import find_top_layer

def get_surface_restart_indiv(Optimizer, indiv):
    """
    Function to generate an structopt Individual class object containing 
        a surface structure from a previously existing structure
    Inputs:
        Optimizer = structopt Optimizer class
        indiv = ASE Atoms object containing the previously existing structure
    Outputs:
        individ = structopt Individual class object containing surface structure data
    """
    #Load surface structure
    surfs = read_xyz(Optimizer.surfacefile)
    cells = Optimizer.surfacecell
    surfs.set_cell(cells)
    surf.set_pbc([True,True,False])
    top,bulks=find_top_layer(indiv,Optimizer.surftopthick)
    individ=Individual(top)
    individ.bulki=bulks.copy()
    individ.bulko=bulks.copy()
    return individ
