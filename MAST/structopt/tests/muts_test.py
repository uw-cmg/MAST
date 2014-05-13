def mutation_test(indiv, A):
    '''Mutation function unit tests
    Key to test is:
        - imports and executes without errors
        - produces distinct individual
    '''
    print 'Beginning unit testing of mutations'
    try:
		from MAST.structopt.moves.ase_minimization import ase_minimization
		A.ase_min_fmax=0.001
		A.ase_min_maxsteps=2500
		ind = indiv.duplicate()
		nind = ase_minimization(ind,A)
		print 'ase_minimization mutation successful'
	except Exception, e:
        print 'ERROR: ase_minimization mutation test FAILED'
        print e
        pass
    try:
		from MAST.structopt.moves.atoms_add import atoms_add
		ind = indiv.duplicate()
		nind = atoms_add(ind,A)
		print 'Atoms_Add mutation successful'
	except Exception, e:
        print 'ERROR: Atoms_Add mutation test FAILED'
        print e
        pass
    try:
		from MAST.structopt.moves.atoms_remove import atoms_remove
		ind = indiv.duplicate()
		nind = atoms_remove(ind,A)
		print 'Atoms_Remove mutation successful'
	except Exception, e:
        print 'ERROR: Atoms_Remove mutation test FAILED'
        print e
        pass
    try:
		from MAST.structopt.moves.basin_hop_la import basin_hop_la
		A.finddefect=False
		A.bh_steps=10
		A.bh_temp=1000*8.617385692256675e-05
		ind = indiv.duplicate()
		nind = basin_hop_la(ind,A)
		print 'basin_hop_la mutation successful'
	except Exception, e:
        print 'ERROR: basin_hop_la mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.basin_hop_permute import basin_hop_permute
        A.finddefect=False
        A.bh_steps=10
        A.bh_temp=1000*8.617385692256675e-05
        ind = indiv.duplicate()
        nind = basin_hop_permute(ind,A)
        print 'basin_hop_permute mutation successful'
	except Exception, e:
        print 'ERROR: basin_hop_permute mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.basin_hop_ra_atoms import basin_hop_ra_atoms
        A.finddefect=False
        A.bh_steps=10
        A.bh_temp=1000*8.617385692256675e-05
        ind = indiv.duplicate()
        nind = basin_hop_ra_atoms(ind,A)
        print 'basin_hop_ra_atoms mutation successful'
	except Exception, e:
        print 'ERROR: basin_hop_ra_atoms mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.basin_hop_rotate import basin_hop_rotate
        A.finddefect=False
        A.bh_steps=10
        A.bh_temp=1000*8.617385692256675e-05
        ind = indiv.duplicate()
        nind = basin_hop_rotate(ind,A)
        print 'basin_hop_rotate mutation successful'
	except Exception, e:
        print 'ERROR: basin_hop_rotate mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.cell_relax_lammps import cell_relax_lammps
        ind = indiv.duplicate()
        nind = cell_relax_lammps(ind,A)
        print 'cell_relax_lammps mutation successful'
	except Exception, e:
        print 'ERROR: cell_relax_lammps mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.cell_shape import cell_shape
        ind = indiv.duplicate()
        nind = cell_shape(ind,A)
        print 'cell_shape mutation successful'
	except Exception, e:
        print 'ERROR: cell_shape mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.lattice_alteration_crystal import lattice_alteration_crystal
        ind = indiv.duplicate()
        nind = lattice_alteration_crystal(ind,A)
        print 'lattice_alteration_crystal mutation successful'
	except Exception, e:
        print 'ERROR: lattice_alteration_crystal mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.lattice_alteration_group import lattice_alteration_group
        ind = indiv.duplicate()
        nind = lattice_alteration_group(ind,A)
        print 'lattice_alteration_group mutation successful'
	except Exception, e:
        print 'ERROR: lattice_alteration_group mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.lattice_alteration_nn import lattice_alteration_nn
        ind = indiv.duplicate()
        nind = lattice_alteration_nn(ind,A)
        print 'lattice_alteration_nn mutation successful'
	except Exception, e:
        print 'ERROR: lattice_alteration_nn mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.lattice_alteration_rdrd import lattice_alteration_rdrd
        ind = indiv.duplicate()
        nind = lattice_alteration_rdrd(ind,A)
        print 'lattice_alteration_rdrd mutation successful'
	except Exception, e:
        print 'ERROR: lattice_alteration_rdrd mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.lattice_alteration_small import lattice_alteration_small
        ind = indiv.duplicate()
        nind = lattice_alteration_small(ind,A)
        print 'lattice_alteration_small mutation successful'
	except Exception, e:
        print 'ERROR: lattice_alteration_small mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.lattice_alteration import lattice_alteration
        ind = indiv.duplicate()
        nind = lattice_alteration(ind,A)
        print 'lattice_alteration mutation successful'
	except Exception, e:
        print 'ERROR: lattice_alteration mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.move_la import move_la
        ind = indiv.duplicate()
        nind = move_la(ind,A)
        print 'move_la mutation successful'
	except Exception, e:
        print 'ERROR: move_la mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.permutation_bulk import permutation_bulk
        ind = indiv.duplicate()
        nind = permutation_bulk(ind,A)
        print 'permutation_bulk mutation successful'
	except Exception, e:
        print 'ERROR: permutation_bulk mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.permutation_crystal_multi import permutation_crystal_multi
        ind = indiv.duplicate()
        nind = permutation_crystal_multi(ind,A)
        print 'permutation_crystal_multi mutation successful'
	except Exception, e:
        print 'ERROR: permutation_crystal_multi mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.permutation_crystal import permutation_crystal
        ind = indiv.duplicate()
        nind = permutation_crystal(ind,A)
        print 'permutation_crystal mutation successful'
	except Exception, e:
        print 'ERROR: permutation_crystal mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.permutation import permutation
        ind = indiv.duplicate()
        nind = permutation(ind,A)
        print 'permutation mutation successful'
	except Exception, e:
        print 'ERROR: permutation mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.quench import quench
        ind = indiv.duplicate()
        nind = quench(ind,A)
        print 'quench mutation successful'
	except Exception, e:
        print 'ERROR: quench mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.random_replacement import random_replacement
        ind = indiv.duplicate()
        nind = random_replacement(ind,A)
        print 'random_replacement mutation successful'
	except Exception, e:
        print 'ERROR: random_replacement mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.rotation_geo import rotation_geo
        ind = indiv.duplicate()
        nind = rotation_geo(ind,A)
        print 'rotation_geo mutation successful'
	except Exception, e:
        print 'ERROR: rotation_geo mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.rotation import rotation
        ind = indiv.duplicate()
        nind = rotation(ind,A)
        print 'rotation mutation successful'
	except Exception, e:
        print 'ERROR: rotation mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.scale_size import scale_size
        ind = indiv.duplicate()
        nind = scale_size(ind,A)
        print 'scale_size mutation successful'
	except Exception, e:
        print 'ERROR: scale_size mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.swap_int_local import swap_int_local
        ind = indiv.duplicate()
        nind = swap_int_local(ind,A)
        print 'swap_int_local mutation successful'
	except Exception, e:
        print 'ERROR: swap_int_local mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.swap_int import swap_int
        ind = indiv.duplicate()
        nind = swap_int(ind,A)
        print 'swap_int mutation successful'
	except Exception, e:
        print 'ERROR: swap_int mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.swap_vacancy import swap_vacancy
        ind = indiv.duplicate()
        nind = swap_vacancy(ind,A)
        print 'swap_vacancy mutation successful'
	except Exception, e:
        print 'ERROR: swap_vacancy mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.swap import swap
        ind = indiv.duplicate()
        nind = swap(ind,A)
        print 'swap mutation successful'
	except Exception, e:
        print 'ERROR: swap mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.zp_rotation_fixed import zp_rotation_fixed
        ind = indiv.duplicate()
        nind = zp_rotation_fixed(ind,A)
        print 'zp_rotation_fixed mutation successful'
	except Exception, e:
        print 'ERROR: zp_rotation_fixed mutation test FAILED'
        print e
        pass
    try:
        from MAST.structopt.moves.zp_rotation import zp_rotation
        ind = indiv.duplicate()
        nind = zp_rotation(ind,A)
        print 'zp_rotation mutation successful'
	except Exception, e:
        print 'ERROR: ZP_Rotation mutation test FAILED'
        print e
        pass
    print 'Mutation Function Testing Complete'
    return
