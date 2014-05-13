import os
try:
    from mpi4py import MPI
except:
    pass
from MAST.structopt.inp_out.write_individual import write_individual

def write_optimizer(Optimizer, optfile, restart=True):
    """Function to write out an Optimizer class object
    Input:
        Optimizer = Optimizer class object to be written out
        optfile = String or fileobject that Optimizer class is written to
    Output:
        No output returned but optimizer object is written to specified file
    """
    if isinstance(optfile, str):
        optfile=open(optfile, 'w')
    convergencesettings = ['convergence_scheme',
        'maxgen',
        'reqrep',
        'convergence',
        'overrideconvergence']
    parametersettings = ['parallel',
        'filename',
        'loggername',
        'atomlist',
        'structure',
        'natoms',
        'optimizer_type',
        'nindiv',
        'genealogy',
        'output_format',
        'allenergyfile',
        'best_inds_list',
        'number_of_bests',
        'indiv_defect_write',
        'vacancy_output',
        'lattice_concentration',
        'postprocessing',
        'genealogytree',
        'seed',
        'forcing',
        'debug',
        'algorithm_type',
        'migration_intervals',
        'migration_percent',
        'fingerprinting',
        'fpbin',
        'fpcutoff',
        'bulkfp',
        'fixed_region',
        'rattle_atoms',
        'constrain_position',
        'restart_ints',
        'restart',
        'r_ab',
        'size',
        'generate_flag',
        'sf',
        'supercell',
        'solidfile',
        'solidcell',
        'evalsolid',
        'finddefects',
        'trackvacs',
        'trackswaps',
        'random_loc_start',
        'random_vac_start',
        'purebulkenpa',
        'natomsbulk',
        'surfacefile',
        'surftopthick',
        'surfacecell',
        'cell_shape_options',
        'alloy',
        'calc_method',
        'vaspcalc',
        'pair_style',
        'bopcutoff',
        'buckcutoff',
        'buckparameters',
        'ps_name',
        'pair_coeff',
        'ps_other',
        'pot_file',
        'lammps_keep_files',
        'lammps_thermo_steps',
        'lammps_min',
        'lammps_min_style',
        'large_box_size',
        'ase_min',
        'ase_min_fmax',
        'ase_min_maxsteps',
        'cxpb',
        'cx_scheme',
        'selection_scheme',
        'mutpb',
        'mutation_options',
        'bh_temp',
        'bh_steps',
        'mutant_add',
        'quench_max_temp',
        'quench_min_temp',
        'quench_n_steps_1',
        'quench_n_steps_2',
        'quench_step_size',
        'isolate_mutation',
        'fitness_scheme',
        'energy_cutoff_factor',
        'stem_coeff',
        'stem_parameters' ,
        'stem_keep_files',
        'constrain_swaps',
        'swaplist',
        'natural_selection_scheme',
        'tournsize',
        'fusslimit',
        'metropolis_temp',
        'tolerance',
        'predator',
        'adaptbegin',
        'adaptmultiplier',
        'demin',
        'mark',
        'restart_optimizer']
    attributelist = ['genrep',
        'minfit',
        'generation',
        'Runtimes',
        'Evaluations',
        'CXs',
        'Muts',
        'cxattempts',
        'mutattempts',
        'optimizerfile']
    Optimizer.restart_optimizer = restart
    #Write Optimizer convergence parameters
    for one in convergencesettings:
        try:
            optpar = eval('Optimizer.{0}'.format(one))
            if isinstance(optpar,str):
                if '\n' in optpar:
                    optfile.write("'{0}':{1},\n".format(one,repr(optpar)))
                else:
                    optfile.write("'{0}':'{1}',\n".format(one,optpar))
            else:
                optfile.write("'{0}':{1},\n".format(one,optpar))
        except:
            optfile.write("'{0}':{1},\n".format(one,None))
            print 'Cannot write parameter: {0}'.format(one)
    #Write Optimizer parameters
    for one in parametersettings:
        try:
            optpar = eval('Optimizer.{0}'.format(one))
            if isinstance(optpar,str):
                if '\n' in optpar:
                    optfile.write("'{0}':{1},\n".format(one,repr(optpar)))
                else:
                    optfile.write("'{0}':'{1}',\n".format(one,optpar))
            else:
                optfile.write("'{0}':{1},\n".format(one,optpar))
        except:
            optfile.write("'{0}':{1},\n".format(one,None))
            print 'Cannot write parameter: {0}'.format(one)
    #Write Optimizer attributes
    for one in attributelist:
        try:
            optpar = eval('Optimizer.{0}'.format(one))
            if isinstance(optpar,str):
                if '\n' in optpar:
                    optifile.write("'{0}':{1},\n".format(one,repr(optpar)))
                else:
                    optfile.write("'{0}':'{1}',\n".format(one,optpar))
            else:
                optfile.write("'{0}':{1},\n".format(one,optpar))
        except:
            optfile.write("'{0}':{1},\n".format(one,None))
            print 'Cannot write attribute: {0}'.format(one)
    #Write Optimizer output files
    flist = []
    try:
        for i in range(len(Optimizer.files)):
            flist.append(Optimizer.files[i].name)
    except:
        flist = []
    optfile.write("'files':{0},\n".format(flist))
    iflist = []
    try:
        if Optimizer.ifiles:
            for i in range(len(Optimizer.ifiles)):
                iflist.append(Optimizer.ifiles[i].name)
        else:
            iflist = None
    except:
        iflist = None
    optfile.write("'ifiles':{0},\n".format(iflist))
    outputfiles = ['tenergyfile','fpfile','fpminfile','debugfile','Genealogyfile','output','summary']
    for one in outputfiles:
        try:
            optfile.write("'{0}':'{1}',\n".format(one,eval('Optimizer.{0}.name'.format(one))))
        except:
            optfile.write("'{0}':{1},\n".format(one,None))
            #print 'Cannot write filename: {0}'.format(one)
    #Write population and bests
    try:
        rank = MPI.COMM_WORLD.Get_rank()
    except:
        rank = 0
    if (len(Optimizer.population) > 0) and (len(Optimizer.BESTS) > 0):
        if isinstance(Optimizer.population[0], str):
            optfile.write("'population':{0},\n".format(Optimizer.population))
            optfile.write("'BESTS':{0}".format(Optimizer.BESTS))
        else:
            fpath = os.path.join(os.getcwd(), '{0}-rank{1}'.format(Optimizer.filename,rank))
            path = os.path.join(fpath,'Restart-files')
            if not os.path.exists(path):
                os.mkdir(path)
            popfiles = []
            for index in range(len(Optimizer.population)):
                fname = os.path.join(path,'Reload-indiv{0:02d}.txt'.format(index))
                write_individual(Optimizer.population[index],fname)
                popfiles.append(fname)
            optfile.write("'population':{0},\n".format(popfiles))
            bestfiles = []
            for index in range(len(Optimizer.BESTS)):
                fname = os.path.join(path,'Reload-bests{0:02d}.txt'.format(index))
                write_individual(Optimizer.BESTS[index],fname)
                bestfiles.append(fname)
            optfile.write("'BESTS':{0}".format(bestfiles))
    else:
        optfile.write("'population':{0},\n".format([]))
        optfile.write("'BESTS':{0}".format([]))
    optfile.close()
    return