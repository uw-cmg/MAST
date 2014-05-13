import os
import copy
import random
import time
import logging
from MAST.structopt.inp_out.loggerUtils import initialize_logger
from MAST.structopt.tools.check_atomlist_concentration import check_atomlist_concentration
import pdb
try:
    from mpi4py import MPI
except:
    pass

def read_parameter_input(input, logger):
    """Function to convert input string, file, or dictionary to a dictionary to contain 
        the parameters for use by the optimizer class.
        input:
            input : Can be dictionary, string, or filename
            logfile : Name of logfile to write to. Default is None
        output:
            parameters : dictionary for defining parameters for optimizer with defaults
    """
    parameters = None
    if isinstance(input,dict):
        # If supplied input is already a dictionary set parameters equal to that dictionary
        parameters=input
    else:
        #Check to see if input is a filename
        if os.path.exists(input):
            try:
                #Check if input file is already formatted as a python dictionary
                parameters = eval('{%s}' % open(input).read())
            except:
                f = open(input,'r')
                lines = f.readlines()
                f.close()
        else:
            #lines = input
            #Split lines in input string
            lines = input.rstrip().split('\n')
        if parameters==None:
            parameters = {}
            for one in lines:
                if one != '\n':
                    p, v = one.strip().split('=')
                    try:
                        #Convert numbers to floats
                        parameters[p.strip()] = float(v.strip())
                    except:
                        try:
                            #Convert list, tuple, and booleans
                            if p.strip() != 'calc_method':
                                parameters[p.strip()] = eval(v.strip())
                            else:
                                parameters[p.strip()] = v.strip()
                        except:
                            try:
                                #Leave remaining as strings
                                parameters[p.strip()] = v.strip()
                            except:
                                print 'Trouble with input line: ', one
    if 'parallel' not in parameters:
        parameters['parallel'] = False
    else:
        parameters['parallel'] = bool(parameters['parallel'])
    if 'filename' not in parameters:
        parameters['filename']='Output'
    else:
        parameters['filename']=str(parameters['filename'])
    try:
        rank = MPI.COMM_WORLD.Get_rank()
    except:
        rank = 0
    if logger:
        if 'loggername' not in parameters:
            parameters['loggername']= '{0}-rank{1}-{2}.log'.format(parameters['filename'],rank,time.strftime("%Y_%m%d_%H%M%S"))
            logger = initialize_logger(parameters['loggername'])
        else:
            parameters['loggername'] = str(parameters['loggername'])
            logger = initialize_logger(parameters['loggername'])
            #logging.getLogger(parameters['loggername'])
    else:
        parameters['loggername'] = None
        logger = dummy_logger_no_write()
    if ('atomlist' not in parameters):
        #Stop program if atom list parameter not in input
        logger.critical("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx',Concentration,Mass,Chemical Potential)]")
        logger.critical("Current parameters include:\n" + repr(parameters))
        raise RuntimeError("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx',Concentration,Mass,Chemical Potential)]")
    else:
        if isinstance(parameters['atomlist'],str):
            parameters['atomlist'] = eval(parameters['atomlist'].strip())
        if not isinstance(parameters['atomlist'],list):
            logger.critical('Something is wrong with atomlist parameter: {0}'.format(parameters['atomlist']))
            raise RuntimeError("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx',Concentration,Mass,Chemical Potential)]: {0}".format(parameters['atomlist']))
    if 'structure' not in parameters:
        #Stop program if structure parameter not in input
        logger.critical("Input file/dictionary must include a structure for the simulation as 'structure':'Cluster/Crystal/Defect'")
        logger.debug("Current parameters include:\n"+repr(parameters))
        raise RuntimeError("Input file/dictionary must include a structure for the simulation as 'structure':'Cluster/Crystal/Defect'")
    for one in parameters['atomlist']:
        if len(one) != 4:
            #Stop program if atom list parameter not properly formatted
            logger.critical('Format of atom list not correct. Must be [(symbol,concentration,mass,potential)]')
            logger.debug('Issue in section : {0}'.format(one))
            logger.debug('Current atomlist is formatted as : {0}'.format(parameters['atomlist']))
            raise RuntimeError('Format of atom list not correct. Must be [(symbol,concentration,mass,potential)]')
    if 'natoms' not in parameters:
        parameters['natoms'] = int(sum([abs(c) for ind,c,m,u in parameters['atomlist']]))
        if rank==0:
            logger.warning('Number of atoms in simulation not set')
            logger.warning('Assuming natoms = {0}'.format(parameters['natoms']))
    else:
        parameters['natoms'] = int(parameters['natoms'])
    parameters['atomlist'] = check_atomlist_concentration(parameters['atomlist'],parameters['natoms'],parameters['loggername'])
    if 'optimizer_type' not in parameters:
        parameters['optimizer_type'] = 'Random'
        if rank==0:
            logger.info('optimizer_type not set.  Default values set to Random.')
    
    #Set parameter defaults based on optimizer structure
    if parameters['optimizer_type'] == 'GA':
        nindiv = 10
        genealogy = True
        nbests = 100
        algtype = 'lambda+mu'
        cxpb = 0.8
        mutpb = 0.15
        natselectscheme = 'tournament'
    else:
        nindiv = 1
        genealogy = False
        nbests = 100
        algtype = parameters['optimizer_type']
        cxpb = 0.0
        mutpb = 1.0
    if parameters['optimizer_type'] == 'SA':
        natselectscheme = 'metropolis'
        predator = 'adapting'
    elif parameters['optimizer_type'] == 'BH':
        natselectscheme = 'metropolis'
    else:
        natselectscheme = 'best'
    if 'nindiv' not in parameters:
        parameters['nindiv'] = nindiv
        if rank==0:
            logger.info('Setting number of individuals in population (nindiv) = {0}'.format(parameters['nindiv']))
    else:
        parameters['nindiv'] = int(parameters['nindiv'])
    
    #Parameters for output
    if 'genealogy' not in parameters:
        parameters['genealogy'] = genealogy
        if rank==0:
            logger.info('Setting genealogy = {0}'.format(parameters['genealogy']))
    else:
        parameters['genealogy'] = bool(parameters['genealogy'])
    if 'output_format' not in parameters:
        parameters['output_format'] = 'fitness'
        if rank==0:
            logger.info('Setting output format = {0}'.format(parameters['output_format']))
    else:
        parameters['output_format'] = str(parameters['output_format']).lower()
    if 'allenergyfile' not in parameters:
        parameters['allenergyfile'] = False
        if rank==0:
            logger.info('Setting allenergyfile = {0}'.format(parameters['allenergyfile']))
    else:
        parameters['allenergyfile'] = bool(parameters['allenergyfile'])
    if 'best_inds_list' not in parameters:
        parameters['best_inds_list'] = True
        if rank==0:
            logger.info('Setting best_inds_list = {0}'.format(parameters['best_inds_list']))
    else:
        parameters['best_inds_list'] = bool(parameters['best_inds_list'])
    if 'number_of_bests' not in parameters:
        parameters['number_of_bests'] = nbests
        if rank==0:
            logger.info('Setting number_of_bests = {0}'.format(parameters['number_of_bests']))
    else:
        parameters['number_of_bests'] = int(parameters['number_of_bests'])
    if 'indiv_defect_write' not in parameters:
        parameters['indiv_defect_write'] = False
        if rank==0:
            logger.info('Setting indiv_defect_write = {0}'.format(parameters['indiv_defect_write']))
    else:
        parameters['indiv_defect_write'] = bool(parameters['indiv_defect_write'])
    if 'vacancy_output' not in parameters:
        parameters['vacancy_output'] = False
        if rank==0:
            logger.info('Setting vacancy_output = {0}'.format(parameters['vacancy_output']))
    else:
        parameters['vacancy_output'] = bool(parameters['vacancy_output'])
    if 'restart_optimizer' not in parameters:
        parameters['restart_optimizer'] = False
    else:
        parameters['restart_optimizer'] = bool(parameters['restart_optimizer'])
    #Parameters for post-processing
    if 'lattice_concentration' not in parameters:
        parameters['lattice_concentration'] = False
        if rank==0:
            logger.info('Setting lattice_concentration = {0}'.format(parameters['lattice_concentration']))
    else:
        parameters['lattice_concentration'] = bool(parameters['lattice_concentration'])
    if 'postprocessing' not in parameters:
        parameters['postprocessing'] = False
        if rank==0:
            logger.info('Setting postprocessing = {0}'.format(parameters['postprocessing']))
    else:
        parameters['postprocessing'] = bool(parameters['postprocessing'])
    if 'genealogytree' not in parameters:
        parameters['genealogytree'] = False
        if rank==0:
            logger.info('Setting genealogytree = {0}'.format(parameters['genealogytree']))
    else:
        parameters['genealogytree'] = bool(parameters['genealogytree'])
    
    #Parameters for general algorithm
    if 'seed' not in parameters:
        parameters['seed']=random.randint(0,10)
        if rank==0:
            logger.info('Setting Random number seed (seed) to {0}'.format(parameters['seed']))
    else:
        parameters['seed'] = int(parameters['seed'])
    if 'forcing' not in parameters:
        parameters['forcing'] = 'Concentration'
        if rank==0:
            logger.info('Setting forcing = {0}'.format(parameters['forcing']))
            logger.info('Assuming forcing concentration control')
    else:
        parameters['forcing'] = str(parameters['forcing'])
    if 'debug' not in parameters:
        parameters['debug'] = ['None']
        if rank==0:
            logger.info('Setting debug = {0}'.format(parameters['debug']))
    else:
        parameters['debug'] = list(parameters['debug'])
        if 'None' not in parameters['debug']: 
            print '***** DEBUGGING RUN *****'
    if 'algorithm_type' not in parameters:
        parameters['algorithm_type'] = algtype
        if rank==0:
            logger.info('Setting algorithm type = {0}'.format(parameters['algorithm_type']))
    else:
        parameters['algorithm_type'] = str(parameters['algorithm_type'])
    if 'migration_intervals' not in parameters:
        parameters['migration_intervals'] = 5
        if rank==0:
            logger.info('Setting migration_intervals = '.format(parameters['migration_intervals']))
    else:
        parameters['migration_intervals'] = int(parameters['migration_intervals'])
    if 'migration_percent' not in parameters:
        parameters['migration_percent'] = 0.05
        if rank==0:
            logger.info('Setting migration_percent = {0}'.format(parameters['migration_percent']))
    else:
        parameters['migration_percent'] = float(parameters['migration_percent'])
    if 'fingerprinting' not in parameters:
        parameters['fingerprinting'] = False
        if rank==0:
            logger.info('Setting fingerprinting = {0}'.format(parameters['fingerprinting']))
    else:
        parameters['fingerprinting'] = bool(parameters['fingerprinting'])
    if 'fpbin' not in parameters:
        parameters['fpbin'] = 0.25
        if rank==0:
            logger.info('Setting fingerprint bin to {0}'.format(parameters['fpbin']))
    else:
        parameters['fpbin'] = float(parameters['fpbin'])
    if 'fpcutoff' not in parameters:
        parameters['fpcutoff'] = 15.0
        if rank==0:
            logger.info('Setting fingerprint cutoff distance to {0}'.format(parameters['fpcutoff']))
    else:
        parameters['fpcutoff'] = float(parameters['fpcutoff'])
    parameters['bulkfp'] = None
    if 'fixed_region' not in parameters:
        parameters['fixed_region'] = False
        if rank==0:
            logger.info('Setting fixed_region = {0}'.format(parameters['fixed_region']))
    else:
        parameters['fixed_region'] = bool(parameters['fixed_region'])
    if 'rattle_atoms' not in parameters:
        parameters['rattle_atoms'] = False
        if rank==0:
            logger.info('Setting rattle_atoms = {0}'.format(parameters['rattle_atoms']))
    else:
        parameters['rattle_atoms'] = bool(parameters['rattle_atoms'])
    if 'constrain_position' not in parameters:
        parameters['constrain_position'] = False
        if rank==0:
            logger.info('Setting constrain_position = {0}'.format(parameters['constrain_position']))
    else:
        parameters['constrain_position'] = bool(parameters['constrain_position'])
    if 'restart' not in parameters:
        parameters['restart'] = False
        if rank==0:
            logger.info('Setting restart = {0}'.format(parameters['restart']))
    else:
        parameters['restart'] = bool(parameters['restart'])
    if 'restart_ints' not in parameters:
        parameters['restart_ints'] = 0
        if rank==0:
            if parameters['restart']:
                logger.info('Setting restart_ints = {0}'.format(parameters['restart_ints']))
    else:
        parameters['restart_ints'] = int(parameters['restart_ints'])
    
    # Parameters to generate the population and individual
    if 'r_ab' not in parameters:
        parameters['r_ab'] = 2.5
        if rank==0:
            logger.info('Setting r_ab = {0}'.format(parameters['r_ab']))
    else:
        parameters['r_ab'] = float(parameters['r_ab'])
    if 'size' not in parameters:
        parameters['size'] = parameters['natoms']**0.33333*parameters['r_ab']
        if rank==0:
            logger.info('Setting size to (natoms)^(1/3)*r_ab = {0}'.format(parameters['size']))
    else:
        parameters['size'] = float(parameters['size'])
    if 'generate_flag' not in parameters:
        parameters['generate_flag']='box'
        if rank==0:
            logger.info('Setting default generation scheme = {0}'.format(parameters['generate_flag']))
    else:
        parameters['generate_flag']=str(parameters['generate_flag']).lower()
    parameters['solidbulk'] = None
    if 'sf' not in parameters:
        parameters['sf'] = 1.75
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting size factor for Defect (sf) = {0}'.format(parameters['sf']))
    else:
        parameters['sf'] = float(parameters['sf'])
    if 'supercell' not in parameters:
        parameters['supercell'] = (1,1,1)
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting supercell for Defect (supercell) = {0}'.format(parameters['supercell']))
    else:
        parameters['supercell'] = tuple(parameters['supercell'])
    if 'solidfile' not in parameters:
        if parameters['structure'] == 'Defect':
            logger.critical('Must provide a file for bulk solid if running a defect simulation.')
            raise RuntimeError('Error: Bulk for Defect not specified. Enter name of file for bulk structure as solidfile parameter')
        else:
            parameters['solidfile'] = None
    else:
        parameters['solidfile'] = str(parameters['solidfile'])
    if 'solidcell' not in parameters:
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.warning('Warning cell size for Bulk Solid not specified assuming distance between 1st and 3rd atom')
        parameters['solidcell'] = None
    else:
        try:
            parameters['solidcell'] = eval('numpy.'+parameters['solidcell'])
        except:
            try:
                parameters['solidcell'] = list(parameters['solidcell'])
            except:
                if rank==0:
                    if parameters['structure'] == 'Defect':
                        logger.warning('Warning cell size for Bulk Solid not recognized assuming distance between 1st and 3rd atom')
                        logger.debug('solidcell input as: {0}'.format(parameters['solidcell']))
                parameters['solidcell'] = None 
    if 'evalsolid' not in parameters:
        parameters['evalsolid'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Not evaluating Solid')
    else:
        parameters['evalsolid'] = bool(parameters['evalsolid'])
    if 'finddefects' not in parameters:
        parameters['finddefects'] = True
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting finddefects = {0}'.format(parameters['finddefects']))
    else:
        parameters['finddefects'] = bool(parameters['finddefects'])
    if 'trackvacs' not in parameters:
        parameters['trackvacs'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting trackvacs = {0}'.format(parameters['trackvacs']))
    else:
        parameters['trackvacs'] = bool(parameters['trackvacs'])
    if 'trackswaps' not in parameters:
        parameters['trackswaps'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting trackswaps = {0}'.format(parameters['trackswaps']))
    else:
        parameters['trackswaps'] = bool(parameters['trackswaps'])
    if 'random_loc_start' not in parameters:
        parameters['random_loc_start'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting random_loc_start = {0}'.format(parameters['random_loc_start']))
    else:
        parameters['random_loc_start'] = bool(parameters['random_loc_start'])
    if 'random_vac_start' not in parameters:
        parameters['random_vac_start'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting random_vac_start = {0}'.format(parameters['random_vac_start']))
    else:
        parameters['random_vac_start'] = bool(parameters['random_vac_start'])
    if 'purebulkenpa' not in parameters:
        if parameters['structure'] == 'Defect':
            parameters['purebulkenpa'] = None
        else:
            parameters['purebulkenpa'] = 0
    if 'natomsbulk' not in parameters:
        if parameters['structure'] == 'Defect':
            parameters['natomsbulk'] = None
        else:
            parameters['natomsbulk']=0
    if 'surfacefile' not in parameters:
        if parameters['structure'] == 'Surface':
            logger.critical('Must provide a file for bulk solid if running a defect simulation.')
            raise RuntimeError('Error: Bulk for Defect not specified. Enter name of file for bulk structure as solidfile parameter')
        else:
            parameters['surfacefile'] = None
    if 'surfacecell' not in parameters:
        parameters['surfacecell'] = None
    if 'surftopthick' not in parameters:
        parameters['surftopthick'] = 0        
    if 'cell_shape_options' not in parameters:
        parameters['cell_shape_options'] = ['cubic', 'hexagonal', 'triclinic', \
        'monoclinic', 'orthorhombic', 'tetragonal']
    else:
        parameters['cell_shape_options'] = list(parameters['cell_shape_options'])
        if parameters['cell_shape_options'] == 'all':
            parameters['cell_shape_options']=['cubic', 'hexagonal', 'triclinic', \
            'monoclinic', 'orthorhombic', 'tetragonal']
    if rank==0:
        if parameters['structure'] == 'Crystal':
            logger.info('Assuming following cell shape options: {0}'.format(parameters['cell_shape_options']))
    if 'alloy' not in parameters:
        parameters['alloy'] = True
        if rank==0:
            logger.info('Setting alloy = {0}'.format(parameters['alloy']))
    else:
        parameters['alloy'] = bool(parameters['alloy'])
    
    #Parameters to Evaluate an Individual
    if 'calc_method' not in parameters:
        if 'vaspcalc' in parameters:
            parameters['calc_method'] = 'VASP'
            if rank==0:
                logger.info('Setting calculation method to VASP. calc_method = {0}'.format(parameters['calc_method']))
        elif 'pair_style' in parameters:
            parameters['calc_method'] = 'LAMMPS'
            if rank==0:
                logger.info('Setting calculation method to LAMMPS. calc_method = {0}'.format(parameters['calc_method']))
        elif 'lammps_min' in parameters:
            parameters['calc_method'] = 'LAMMPS'
            if rank==0:
                logger.info('Setting calculation method to LAMMPS. calc_method = {0}'.format(parameters['calc_method']))
        else:
            parameters['calc_method'] = 'LennardJones'
            if rank==0:
                logger.info('Calculation method not specified assuming ase LennardJones. calc_method = {0}'.format(parameters['calc_method']))
    else:
        parameters['calc_method'] = str(parameters['calc_method'])
    if 'vaspcalc' not in parameters:
        parameters['vaspcalc'] = "Vasp()"
        if rank==0:
            if parameters['calc_method'] == 'VASP':
                logger.info('Setting vaspcalc = {0}'.format(parameters['vaspcalc']))
    if 'pair_style' not in parameters:
        if rank==0:
            if parameters['calc_method'] == 'LAMMPS':
                logger.info('No pair_style for LAMMPS specified. Assuming Lennard-Jones')
        parameters['pair_style'] = None
    if 'bopcutoff' not in parameters:
        if parameters['pair_style'] == 'bop':
            logger.critical('BOP potential requires a cutoff distance.')
            logger.debug('Parameters = {0}'.format(parameters))
            raise RuntimeError("ERROR:Cutoff distance must be supplied with bop potential! Specify with bopcutoff=N")
        else:
            parameters['bopcutoff'] = 0
    else:
        parameters['bopcutoff'] = float(parameters['bopcutoff'])
    if 'buckcutoff' not in parameters:
        parameters['buckcutoff'] = 1
        if rank==0:
            if parameters['pair_style'] == 'buck':
                logger.warning('Setting buckingham potential cutoff = ' +repr(parameters['buckcutoff']))
    else:
        parameters['buckcutoff'] = int(parameters['buckcutoff'])
    if 'buckparameters' not in parameters:
        parameters['buckparameters'] = ['* * 100.00 1.5 200.0']
        if rank == 0:
            if parameters['pair_style'] == 'buck':
                logger.warning('Setting Buckingham parameters = {0}'.format(parameters['buckparameters']))
    if 'ps_name' not in parameters:
        parameters['ps_name'] = None
        if rank == 0:
            if parameters['pair_style'] == 'other':
                logger.info('Setting ps_name = {0}'.format(parameters['ps_name']))
    if 'pair_coeff' not in parameters:
        if parameters['pair_style'] == 'other':
            logger.critical('Must provide pair_coeff for use with pair_style = other')
            raise RuntimeError("ERROR:Coefficients for structure=other potential not specifed. Use pair_coeff")
        else:
            parameters['pair_coeff'] = None
    if 'ps_other' not in parameters:
        parameters['ps_other'] = None
        if rank==0:
            if parameters['pair_style'] == 'other':
                logger.info('Setting ps_other = {0}'.format(parameters['ps_other']))
    if 'pot_file' not in parameters:
        parameters['pot_file'] = None
        if rank == 0:
            if parameters['calc_method'] == 'LAMMPS':
                logger.warning('No potential file for LAMMPS specified.')
    else:
        parameters['pot_file'] = str(parameters['pot_file'])
    if 'lammps_keep_files' not in parameters:
        parameters['lammps_keep_files'] = False
        if rank == 0:
            if parameters['calc_method'] == 'LAMMPS':
                logger.info('Setting lammps_keep_files = {0}'.format(parameters['lammps_keep_files']))
    else:
        parameters['lammps_keep_files'] = bool(parameters['lammps_keep_files'])
    if 'lammps_thermo_steps' not in parameters:
        parameters['lammps_thermo_steps'] = 1
        if rank == 0:
            if parameters['calc_method'] == 'LAMMPS':
                logger.info('Setting lammps_thermo_steps = {0}'.format(parameters['lammps_thermo_steps']))
    else:
        parameters['lammps_thermo_steps'] = int(parameters['lammps_thermo_steps'])
    if ('lammps_min' in parameters) and (parameters['lammps_min'] != None):
        parameters['ase_min'] = False
    else:
        parameters['lammps_min'] = None
        if rank ==0:
            if parameters['calc_method'] == 'LAMMPS':
                logger.info('No Local minimization implemented in LAMMPS')
    if 'lammps_min_style' not in parameters:
        if parameters['lammps_min']:
            parameters['lammps_min_style'] = 'cg'
            if rank ==0:
                if parameters['calc_method'] == 'LAMMPS':
                    logger.info('Setting lammps_min_style = {0}'.format('lammps_min_style')+'. Using LAMMPS conjugate gradient local minimizer')
        else:
            parameters['lammps_min_style'] = None
    if 'large_box_size' not in parameters:
        parameters['large_box_size']=500.0
        if rank == 0:
            if parameters['structure']=='Cluster':
                logger.info('Setting large_box_size to {0}'.format(parameters['large_box_size']))
    else:
        parameters['large_box_size']=float(parameters['large_box_size'])
    if 'ase_min' not in parameters:
        if parameters['calc_method']=='LennardJones':
            parameters['ase_min'] = True
        else:
            parameters['ase_min'] = False
        if rank==0:
            logger.info('Setting ase_min = {0}'.format(parameters['ase_min']))
    else:
        parameters['ase_min'] = bool(parameters['ase_min'])
    if 'ase_min_fmax' not in parameters:
        parameters['ase_min_fmax'] = 0.01
        if rank ==0:
            if parameters['ase_min']:
                logger.info('Setting ASE maximum force value (ase_min_fmax) = {0}'.format(parameters['ase_min_fmax']))
    else:
        parameters['Sovlerfmax'] = float(parameters['ase_min_fmax'])
    if 'ase_min_maxsteps' not in parameters:
        parameters['ase_min_maxsteps'] = 2500
        if rank ==0:
            if parameters['ase_min']:
                logger.info('Setting maximum number of steps for BFGS solver (ase_min_maxsteps) = {0}'.format(parameters['ase_min_maxsteps']))
    else:
        parameters['SovlerMxSteps'] = int(parameters['ase_min_maxsteps'])
    
    #Parameters for Crossovers
    if 'cxpb' not in parameters:
        parameters['cxpb'] = cxpb
        if rank == 0:
            logger.info('Setting crossover probability (cxpb) = {0}'.format(parameters['cxpb']))
    else:
        parameters['cxpb'] = float(parameters['cxpb'])
    if 'cx_scheme' not in parameters:
        parameters['cx_scheme'] = 'cxtp'
        if rank == 0:
            logger.info('Assuming two-point crossover.  Setting cx_scheme = {0}'.format(parameters['cx_scheme']))
    else:
        parameters['cx_scheme'] = str(parameters['cx_scheme']).lower()
    if 'selection_scheme' not in parameters:
        parameters['selection_scheme'] = 'tournament2'
        if rank == 0:
            logger.info('Setting selection_scheme = {0}'.format(parameters['selection_scheme']))
    else:
        parameters['selection_scheme'] = str(parameters['selection_scheme']).lower()
    
    #Parameters for Mutations
    if 'mutpb' not in parameters:
        parameters['mutpb'] = mutpb
        if rank == 0:
            logger.info('Setting mutation probability (mutpb) = {0}'.format(parameters['mutpb']))
    else:
        parameters['mutpb'] = float(parameters['mutpb'])
    if 'mutation_options' not in parameters:
        if parameters['structure']=='Cluster':
            parameters['mutation_options']=['lattice_alteration','rotation',\
            'permutation','scale_size']
        elif parameters['structure']=='Crystal':
            parameters['mutation_options']=['lattice_alteration','rotation',\
            'permutation','scale_size', 'cell_shape', 'lammps_box_relax']
        elif parameters['structure']=='Defect':
            parameters['mutation_options']=['lattice_alteration','rotation','permutation']
        if rank == 0:
            logger.info('Setting mutations options = {0}'.format(parameters['mutation_options']))
    else:
        parameters['mutation_options'] = list(parameters['mutation_options'])
        for i in range(len(parameters['mutation_options'])):
            parameters['mutation_options'][i] = parameters['mutation_options'][i].lower()
    BHFlag=False
    for one in parameters['mutation_options']:
        if 'basin_hop' in one:
            BHFlag=True
    if 'bh_steps' not in parameters:
        parameters['bh_steps']=100
        if rank == 0:
            if BHFlag==True:
                logger.warning('Max steps not specified for Basin Hop mutation, setting bh_steps = {0}'.format(parameters['bh_steps']))
    else:
        parameters['bh_steps'] = int(parameters['bh_steps'])
    if 'bh_temp' not in parameters:
        parameters['bh_temp'] = 1000*8.617385692256675e-05
        if rank == 0:
            if BHFlag==True:
                logger.warning('Temperature not set for Basin Hop mutation, setting bh_temp in kT = {0}'.format(parameters['bh_temp']))
    else:
        parameters['bh_temp'] = float(parameters['bh_temp'])
    if 'mutant_add' not in parameters:
        parameters['mutant_add'] = False
        if rank == 0:
            logger.info('Setting mutant_add = {0}'.format(parameters['mutant_add']))
    else:
        parameters['mutant_add'] = bool(parameters['mutant_add'])
    if 'quench_max_temp' not in parameters:
        parameters['quench_max_temp'] = 1000
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_max_temp = {0}'.format(parameters['quench_max_temp']))
    else:
        parameters['quench_max_temp'] = int(parameters['quench_max_temp'])
    if 'quench_min_temp' not in parameters:
        parameters['quench_min_temp'] = 2
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Seting quench_min_temp = {0}'.format(parameters['quench_min_temp']))
    else:
        parameters['quench_min_temp'] = int(parameters['quench_min_temp'])
    if 'quench_step_size' not in parameters:
        parameters['quench_step_size'] = 0.01
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_step_size = {0}'.format(parameters['quench_step_size']))
    else:
        parameters['quench_step_size'] = float(parameters['quench_step_size'])
    if 'quench_n_steps_1' not in parameters:
        parameters['quench_n_steps_1'] = 10000
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_n_steps_1 = {0}'.format(parameters['quench_n_steps_1']))
    else:
        parameters['quench_n_steps_1'] = int(parameters['quench_n_steps_1'])
    if 'quench_n_steps_2' not in parameters:
        parameters['quench_n_steps_2'] = parameters['quench_n_steps_1']*2
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_n_steps_2 = {0}'.format(parameters['quench_n_steps_2']))
    else:
        parameters['quench_n_steps_2'] = int(parameters['quench_n_steps_2'])
    if 'isolate_mutation' not in parameters:
        parameters['isolate_mutation'] = False
        if rank == 0:
            logger.info('Setting isolate_mutation flag = {0}'.format(parameters['isolate_mutation']))
    else:
        parameters['isolate_mutation'] = bool(parameters['isolate_mutation'])
    
    #Parameters for Selection
    if 'fitness_scheme' not in parameters:
        parameters['fitness_scheme'] = 'totalenfit'
        if rank == 0:
            logger.info('Setting fitness_scheme = {0}'.format(parameters['fitness_scheme']))
    else:
        parameters['fitness_scheme'] = str(parameters['fitness_scheme']).lower()
    if 'energy_cutoff_factor' not in parameters:
        parameters['energy_cutoff_factor'] = 10.0
        if rank ==0:
            logger.info('Setting energy_cutoff_factor = {0}'.format(parameters['energy_cutoff_factor']))
    else:
        parameters['energy_cutoff_factor'] = float(parameters['energy_cutoff_factor'])
    if 'stem_parameters' not in parameters:
        if 'stem' in parameters['fitness_scheme']:
            logger.critical('Must provide stem_parameters for STEM_Cost calculation')
            raise RuntimeError("STEM parameters not specified.  Cannot simulate image files")
        else:
            parameters['stem_parameters'] = {}
    else:
        stemparams = dict(parameters['stem_parameters'])
    if 'stem_keep_files' not in parameters:
        parameters['stem_keep_files'] = True
        if rank == 0:
            if 'stem' in parameters['fitness_scheme']:
                logger.info('Setting stem_keep_files = {0}'.format(parameters['stem_keep_files']))
    else:
        parameters['stem_keep_files'] = parameters['stem_keep_files']
    if 'stem_coeff' not in parameters:
        parameters['stem_coeff'] = None
        if rank == 0:
            if 'stem' in parameters['fitness_scheme']:
                logger.info('Setting stem_coeff with first individual')
    else:
        try:
            parameters['stem_coeff'] = float(parameters['stem_coeff'])
        except:
            if rank == 0:
                if 'stem' in parameters['fitness_scheme']:
                    logger.warning('Trouble reading stem_coeff input. stem_coeff = {0}'.format(parameters['stem_coeff']))
            parameters['stem_coeff'] = None
    if 'stem' in parameters['fitness_scheme']:
        #Initialize function for experimental image
        from MAST.structopt.tools.StemCalc import ConvStem
        logger.info('Initializing ConvStem Calculator')
        parameters['stemcalc'] = ConvStem(parameters=stemparams, tmp_dir='/'+os.getcwd()+'/ConvStemImages/', keep_files=parameters['stem_keep_files'])
    else:
        parameters['stemcalc'] = None
    if 'constrain_swaps' not in parameters:
        if 'IntSwap' in parameters['mutation_options']:
            parameters['swaplist'] = None
            if rank == 0:
                logger.info('Setting swaplist for IntSwap = None')
        elif parameters['fitness_scheme'] == 'chempotswap':
            parameters['swaplist'] = None
            if rank == 0:
                logger.info('Setting swaplist for chempotswap = None')
        else:
            parameters['swaplist'] = False
            parameters['constrain_swaps'] = False
            if rank == 0:
                logger.info('Setting swaplist = False')
    else:
        parameters['swaplist'] = parameters['constrain_swaps']
    if 'natural_selection_scheme' not in parameters:
        parameters['natural_selection_scheme'] = natselectscheme
        if rank ==0:
            logger.info('Setting natural_selection_scheme = {0}'.format(parameters['natural_selection_scheme']))
    else:
        parameters['natural_selection_scheme'] = str(parameters['natural_selection_scheme']).lower()
    if 'tournsize' not in parameters:
        parameters['tournsize'] = 3
        if 'tournament' in parameters['selection_scheme'] or 'tournament' in parameters['natural_selection_scheme']:
            if rank ==0:
                logger.info('Setting Tournament size (tournsize) = {0}'.format(parameters['tournsize']))
    else:
        parameters['tournsize']=int(parameters['tournsize'])
    if 'fusslimit' not in parameters:
        parameters['fusslimit'] = 10.0
        if rank ==0:
            logger.info('Setting FUSS limit (fusslimit) = {0}'.format(parameters['fusslimit']))
    else:
        parameters['fusslimit'] = float(parameters['fusslimit'])
    if 'metropolis_temp' not in parameters:
        parameters['metropolis_temp'] = 30.0
        if rank == 0:
            logger.info('Setting metropolis_temp = {0}'.format(parameters['metropolis_temp']))
    else:
        parameters['metropolis_temp'] = float(parameters['metropolis_temp'])
    if 'mark' not in parameters:
        parameters['mark'] = None
    #Parameters for Convergence
    if 'convergence_scheme' not in parameters:
        parameters['convergence_scheme'] = 'max_gen'
        if rank ==0:
            logger.info('Setting convergence scheme (convergence_scheme) = {0}'.format(parameters['convergence_scheme']))
    else:
        parameters['convergence_scheme'] = str(parameters['convergence_scheme']).lower()
    if 'maxgen' not in parameters:
        parameters['maxgen'] = 5
        if rank ==0:
            logger.info('Setting Max Number of generations (maxgen) = {0}'.format(parameters['maxgen']))
    else:
        parameters['maxgen'] = int(parameters['maxgen'])
    if 'reqrep' not in parameters:
        parameters['reqrep'] = 10
        if rank ==0:
            if 'rep' in parameters['convergence_scheme']:
                logger.info('Setting max number of energy repetitions (reqrep) = {0}'.format(parameters['reqrep']))
    else:
        parameters['reqrep'] = int(parameters['reqrep'])
    if 'tolerance' not in parameters:
        parameters['tolerance'] = 0.001
        if rank ==0:
            if 'rep' in parameters['convergence_scheme']:
                logger.info('Setting energy tolerance (tolerance) = {0}'.format(parameters['tolerance']))
    else:
        parameters['tolerance'] = float(parameters['tolerance'])
    if 'predator' not in parameters:
        parameters['predator'] = 'mutation_dups'
        if rank == 0:
            logger.info('Setting predator = {0}'.format(parameters['predator']))
    else:
        parameters['predator'] = str(parameters['predator']).lower()
    if 'adaptbegin' not in parameters:
        parameters['adaptbegin'] = 0.75
        if rank == 0:
            if parameters['predator'] == 'adapting':
                logger.info('Setting adaptation predator to begin (adaptbegin) at genrep*{0}'.format(parameters['adaptbegin']))
    else:
        parameters['adaptbegin'] = float(parameters['adaptbegin'])
    if 'adaptmultiplier' not in parameters:
        parameters['adaptmultiplier'] = 3.0
        if rank == 0:
            if parameters['predator'] == 'adapting':
                logger.info('Setting adaptation predator multiplier (adaptmultiplier) = {0}'.format(parameters['adaptmulitplier']))
    else:
        parameters['adaptmultiplier'] = float(parameters['adaptmultiplier'])
    if 'demin' not in parameters:
        parameters['demin'] = 0.005
        if rank == 0:
            logger.info('Setting cutoff convergence energy (demin) = {0}'.format(parameters['demin']))
    else:
        parameters['demin'] = float(parameters['demin'])
    
    return parameters

class dummy_logger():
    def __init__(self):
        return
    def critical(self,message):
        print 'CRITICAL: {0}'.format(message)
        return
    def debug(self,message):
        print 'DEBUG: {0}'.format(message)
        return
    def warning(self,message):
        print 'WARNING: {0}'.format(message)
        return
    def info(self,message):
        print 'MESSAGE: {0}'.format(message)
        return

class dummy_logger_no_write():
    def __init__(self):
        return
    def critical(self,message):
        return
    def debug(self,message):
        return
    def warning(self,message):
        return
    def info(self,message):
        return