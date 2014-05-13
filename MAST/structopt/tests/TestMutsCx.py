
def test_mutations(indiv):
    #Test BH_LA
    from MAST.structopt.moves.Basin_Hop_LA import Basin_Hop_LA
    A.FindDefects=False
    A.BHsteps=10
    A.BHtemp=1000*8.617385692256675e-05
    ind = indiv.duplicate()
    nind = Basin_Hop_LA(ind,A,False)
    
    #Test BH_Permute
    from MAST.structopt.moves.Basin_Hop_Permute import Basin_Hop_Permute
    A.FindDefects=False
    A.BHsteps=10
    A.BHtemp=1000*8.617385692256675e-05
    ind = indiv.duplicate()
    nind = Basin_Hop_Permute(ind,A,False)
    
    #Test BH_R_Aatoms
    from MAST.structopt.moves.Basin_Hop_R_Aatoms import Basin_Hop_R_Aatoms
    A.FindDefects=False
    A.BHsteps=10
    A.BHtemp=1000*8.617385692256675e-05
    ind = indiv.duplicate()
    nind = Basin_Hop_R_Aatoms(ind,A,False)
    
    #Test BH_rotate
    from MAST.structopt.moves.Basin_Hop_Rotate import Basin_Hop_Rotate
    A.FindDefects=False
    A.BHsteps=10
    A.BHtemp=1000*8.617385692256675e-05
    ind = indiv.duplicate()
    nind = Basin_Hop_Rotate(ind,A,False)
    
    #Test LAC
    from MAST.structopt.moves.Lattice_Alteration_Crystal import Lattice_Alteration_Crystal
    ind = indiv.duplicate()
    nind = Lattice_Alteration_Crystal(ind,A,False)
    
    #Test LAG
    from MAST.structopt.moves.Lattice_Alteration_Group import Lattice_Alteration_Group
    ind = indiv.duplicate()
    nind = Lattice_Alteration_Group(ind,A,False)
    
    #Test LAnn
    from MAST.structopt.moves.Lattice_Alteration_nn import Lattice_Alteration_nn
    ind = indiv.duplicate()
    nind = Lattice_Alteration_nn(ind,A,False)
    
    #Test LArdrd
    from MAST.structopt.moves.Lattice_Alteration_rdrd import Lattice_Alteration_rdrd
    ind = indiv.duplicate()
    nind = Lattice_Alteration_rdrd(ind,A,False)
    
    #Test LAsmall
    from MAST.structopt.moves.Lattice_Alteration_small import Lattice_Alteration_small
    ind = indiv.duplicate()
    nind = Lattice_Alteration_small(ind,A,False)
    
    #Test LA
    from MAST.structopt.moves.Lattice_Alteration import Lattice_Alteration
    ind = indiv.duplicate()
    nind = Lattice_Alteration(ind,A,False)
    
    #Test Move_LA
    from MAST.structopt.moves.Move_LA import Move_LA
    ind = indiv.duplicate()
    nind = Move_LA(ind,A,False)
    
    #Test Permutation
    from MAST.structopt.moves.Permutation import Permutation
    ind = indiv.duplicate()
    nind = Permutation(ind,A,False)
    
    #Test Random_Replacement
    from MAST.structopt.moves.Random_Replacement import Random_Replacement
    ind = indiv.duplicate()
    nind = Random_Replacement(ind,A,False)
    
    #Test Rotation_geo
    from MAST.structopt.moves.Rotation_geo import Rotation_geo
    ind = indiv.duplicate()
    nind = Rotation_geo(ind,A,False)
    
    #Test Rotation
    from MAST.structopt.moves.Rotation import Rotation
    ind = indiv.duplicate()
    nind = Rotation(ind,A,False)
    
    #Test Scale_size
    from MAST.structopt.moves.Scale_size import Scale_size
    ind = indiv.duplicate()
    nind = Scale_size(ind,A,False)
    
    #Test Swap_Int_Local
    from MAST.structopt.moves.Swap_Int_Local import Swap_Int_Local
    ind = indiv.duplicate()
    nind = Swap_Int_Local(ind,A,False)
    
    #Test Swap_Int
    from MAST.structopt.moves.Swap_Int import Swap_Int
    ind = indiv.duplicate()
    nind = Swap_Int(ind,A,False)
    
    #Test Swap_Vacancy
    from MAST.structopt.moves.Swap_Vacancy import Swap_Vacancy
    ind = indiv.duplicate()
    nind = Swap_Vacancy(ind,A,False)
    
    #Test Swap
    from MAST.structopt.moves.Swap import Swap
    ind = indiv.duplicate()
    nind = Swap(ind,A,False)
    
    #Test ZP_Rotation_Fixed
    from MAST.structopt.moves.ZP_Rotation_Fixed import ZP_Rotation_Fixed
    ind = indiv.duplicate()
    nind = ZP_Rotation_Fixed(ind,A,False)
    
    #Test ZP_Rotation
    from MAST.structopt.moves.ZP_Rotation import ZP_Rotation
    ind = indiv.duplicate()
    nind = ZP_Rotation(ind,A,False)

def test_cx(ind1, ind2):
    #Test cxTP
    from MAST.structopt.crossover.cxTP import cxTP
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = cxTP(i1,i2,A,False)
    
    #Test cxTPA
    from MAST.structopt.crossover.cxTPA import cxTPA
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = cxTPA(i1,i2,A,False)
    
    #Test cxTPC
    from MAST.structopt.crossover.cxTPC import cxTPC
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = cxTPC(i1,i2,A,False)
    
    #Test NewClus
    from MAST.structopt.crossover.NewClus import NewClus
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = NewClus(i1,i2,A,False)
    
    #Test randalloybox
    from MAST.structopt.crossover.randalloybox import randalloybox
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = randalloybox(i1,i2,A,False)
    
    #Test rotct_rand_clus
    from MAST.structopt.crossover.rotct_rand_clus import rotct_rand_clus
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = rotct_rand_clus(i1,i2,A,False)
    
    #Test rotct_rand
    from MAST.structopt.crossover.rotct_rand import rotct_rand
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = rotct_rand(i1,i2,A,False)
    
    #Test rotct
    from MAST.structopt.crossover.rotct import rotct
    i1 = ind1.duplicate()
    i2 = ind2.duplicate()
    nc1, nc2 = rotct(i1,i2,A,False)

def main():
    from MAST.structopt import Optimizer
    from MAST.structopt import tools
    
    #Test for cluster structures
    parameters = {'structure':'Cluster','optimizer_type': 'GA','atomlist':[('Au', 10, 0, 0), ('Cu', 10, 0, 0)]}
    A = Optimizer(parameters)
    A.algorithm_initialize()
    A.calc = tools.setup_calculator(A)
    offspring = A.generation_set([])
    indiv = offspring[0].duplicate()
    ind1 = offspring[1].duplicate()
    ind2 = offspring[2].duplicate()
    test_mutations(indiv)
    test_cx(ind1,ind2)
    A.close_output()
    
    #Test for defect structures
    parameters = {'structure':'Defect','optimizer_type': 'GA','atomlist':[('Au', 3, 0, 0), ('Cu', 3, 0, 0)],'supercell':(3,3,3),'SolidFile':'BFe.xyz','SolidCell':[8.61, 8.61, 8.61]}
    A = Optimizer(parameters)
    A.algorithm_initialize()
    A.calc = tools.setup_calculator(A)
    offspring = A.generation_set([])
    indiv = offspring[0].duplicate()
    ind1 = offspring[1].duplicate()
    ind2 = offspring[2].duplicate()
    test_mutations(indiv)
    test_cx(ind1,ind2)
    A.close_output()
    
    #Test for crystal structures
    parameters = {'structure':'Crystal','optimizer_type': 'GA','atomlist':[('Au', 10, 0, 0), ('Cu', 10, 0, 0)]}
    A = Optimizer(parameters)
    A.algorithm_initialize()
    A.calc = tools.setup_calculator(A)
    offspring = A.generation_set([])
    indiv = offspring[0].duplicate()
    ind1 = offspring[1].duplicate()
    ind2 = offspring[2].duplicate()
    test_mutations(indiv)
    test_cx(ind1,ind2)
    A.close_output()

