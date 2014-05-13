import time
import os

def write_parameters(Optimizer):
    """
    Subprogram to write input parameters to output file
    Inputs:
        Optimizer = structopt Optimizer class object with output file data
    Outputs:
        None. Data is written directly to output files provided by Optimizer object
    """
    Optimizer.output.write('-----StructOpt - Global Structure Optimization Run----\n')
    if 'None' not in Optimizer.debug: Optimizer.output.write('*******DEBUGGING RUN*******\n')
    Optimizer.output.write(Optimizer.__version__ + '\n')
    localtime = time.asctime( time.localtime(time.time()) )
    Optimizer.output.write('Local time : ' + repr(localtime) + '\n')

    section_break = '\n--------------------\n'
    
    Optimizer.output.write(section_break + 'General Algorithm Information\n\n')
    if Optimizer.parallel:
        Optimizer.output.write('Algorithm Performed in parallel\n')
    else:
        Optimizer.output.write('Algorithm Performed in serial\n')
    Optimizer.output.write('Algorithm type : ' + Optimizer.algorithm_type + '\n')
    if 'Island_Method' in Optimizer.algorithm_type:
        Optimizer.output.write('    Migration Intervals : ' + repr(Optimizer.migration_intervals) + '\n')
        Optimizer.output.write('    Migration Percent : ' + repr(Optimizer.migration_percent) + '\n')
    Optimizer.output.write('Type of structure to optimize (structure) = ' + repr(Optimizer.structure) + '\n')
    Optimizer.output.write('Random number seed = ' + repr(Optimizer.seed) + '\n')
    if Optimizer.forcing=='Concentration':
        Optimizer.output.write('Concentration in individual held constant through crossover\n')
    elif Optimizer.forcing=='Energy_bias':
        Optimizer.output.write('Concentration in individual held constant through energy bias\n')
    elif Optimizer.forcing=='Chem_pot':
        Optimizer.output.write('Concentration in individual allowed to change based on Chemical Potentials\n')
    else:
        Optimizer.output.write('Concentration not held constant')
    Optimizer.output.write('Structure Fingerprinting (fingerprinting) : ' + repr(Optimizer.fingerprinting) + '\n')
    if Optimizer.fingerprinting:
        Optimizer.output.write('Fingerprint bin size (fpbin) : ' + repr(Optimizer.fpbin) + '\n')
        Optimizer.output.write('Fingerprint cutoff distance (fpcutoff) : ' + repr(Optimizer.fpcutoff) + '\n')
    if Optimizer.fixed_region: Optimizer.output.write('Fixed Bulk calculation \n')
    if Optimizer.constrain_position: Optimizer.output.write('Constrained position calculation \n')

    Optimizer.output.write(section_break + 'Population Generation Parameters\n\n')
    if Optimizer.restart: 
        Optimizer.output.write('Restarting Old Run\n')
        if Optimizer.restart_ints !=0: Optimizer.output.write('Assuming specified restart ints = ' + repr(Optimizer.restart_ints) + '\n')
    Optimizer.output.write('Number of structures in population (nindiv) : '  + repr(Optimizer.nindiv) + '\n')
    Optimizer.output.write('Starting structure configuration (atomlist) : '+ repr(Optimizer.atomlist) + '\n')
    Optimizer.output.write('Number of atoms in initial structure (natoms) : ' + repr(Optimizer.natoms) + '\n')
    Optimizer.output.write('Starting structure size (size) : ' + repr(Optimizer.size) + '\n')
    Optimizer.output.write('Average distance between atoms (r_ab) : ' + repr(Optimizer.r_ab) + '\n')
    Optimizer.output.write('Generation flag type = '+repr(Optimizer.generate_flag)+'\n')
    if Optimizer.structure=='Defect':
        Optimizer.output.write('File for bulk solid for Defect run (solidfile) : ' + repr(Optimizer.solidfile) + '\n')
        Optimizer.output.write('Cell size for bulk solid for Defect run (solidcell) : ' + repr(Optimizer.solidcell) + '\n')
        Optimizer.output.write('Size of bulk solid supercell (supercell) : ' + repr(Optimizer.supercell) + '\n')
        Optimizer.output.write('Evaluation of energy of bulk solid (evalsolid) : ' + repr(Optimizer.evalsolid) + '\n')
        Optimizer.output.write('Scale factor for surrounding material to include in Defect run (sf) : ' + repr(Optimizer.sf) + '\n')
        Optimizer.output.write('Start from random location option (Random_loc_start) : ' + repr(Optimizer.random_loc_start) + '\n')
        Optimizer.output.write('Start from random location option (random_vac_start) : ' + repr(Optimizer.random_vac_start) + '\n')
        Optimizer.output.write('Interstitial finding scheme is (finddefects) : ' + repr(Optimizer.finddefects) + '\n')
        Optimizer.output.write('Vacancy finding scheme is (trackvacs) : ' + repr(Optimizer.trackvacs) + '\n')
        Optimizer.output.write('Substitutions finding scheme is (trackswaps) : ' + repr(Optimizer.trackswaps) + '\n')
        if Optimizer.alloy:
            Optimizer.output.write('Allow free exchange between bulk and cluster (alloy) : ' + repr(Optimizer.alloy) + '\n')
    elif Optimizer.structure=='Crystal':
        Optimizer.output.write('Crystal cell shape options (cell_shape_options) : ' + repr(Optimizer.cell_shape_options) + '\n')

    Optimizer.output.write(section_break + 'Crossover Configuration Parameters\n\n')
    Optimizer.output.write('Crossover probability (cxpb) : ' + repr(Optimizer.cxpb) + '\n')
    Optimizer.output.write('Crossover scheme (cx_scheme) : ' + repr(Optimizer.cx_scheme) + '\n')
    Optimizer.output.write('Parent selection scheme (selection_scheme) : ' + repr(Optimizer.selection_scheme) + '\n')
    if 'tournament' in Optimizer.selection_scheme: Optimizer.output.write('    Tournament Size (tournsize) : ' + repr(Optimizer.tournsize) + '\n')
    if 'fuss' in Optimizer.selection_scheme:
        if Optimizer.selection_scheme != 'fuss1': 
            Optimizer.output.write('    fuss limit (fusslimit) : ' + repr(Optimizer.fusslimit) + ' eV \n')
        if Optimizer.natural_selection_scheme == 'fussf':
            Optimizer.output.write('Fingerprint bin size (fpbin) : ' + repr(Optimizer.fpbin) + '\n')
            Optimizer.output.write('Fingerprint cutoff distance (fpcutoff) : ' + repr(Optimizer.fpcutoff) + '\n')

    Optimizer.output.write(section_break + 'Mutation Configuration Parameters\n\n')
    Optimizer.output.write('Mutation Probability (mutpb) : ' + repr(Optimizer.mutpb) + '\n')
    Optimizer.output.write('Mutation Options (mutation_options) : ' + repr(Optimizer.mutation_options) + '\n')
    if Optimizer.swaplist==True or Optimizer.swaplist==None: Optimizer.output.write('Allow Constrained Swapping of atoms is True\n')
    BHFlag=False
    for one in Optimizer.mutation_options:
        if 'basin_hop' in one:
            BHFlag=True
    if BHFlag==True:
        Optimizer.output.write('Basin Hop mutation max steps (bh_steps) : ' + repr(Optimizer.bh_steps) + '\n')
        Optimizer.output.write('Basing Hop mutation temperature (bh_temp) : ' + repr(Optimizer.bh_temp) + '\n')
    if 'quench' in Optimizer.mutation_options:
        Optimizer.output.write('Quench mutation max temperature (K)' + repr(Optimizer.quench_max_temp)+'\n')
        Optimizer.output.write('Quench mutation min temperature (K)' +repr(Optimizer.quench_min_temp)+'\n')
        Optimizer.output.write('Quench mutation step size' +repr(Optimizer.quench_step_size)+'\n')
        Optimizer.output.write('Quench mutation number of steps 1'+repr(Optimizer.quench_n_steps_1)+'\n')
        Optimizer.output.write('Quench mutation number of steps 2'+repr(Optimizer.quench_n_steps_2)+'\n')
    if Optimizer.structure=='Defect':
        Optimizer.output.write('Setting mutation isolation to only region 1 = '+repr(Optimizer.isolate_mutation)+'\n')

    Optimizer.output.write(section_break + 'Individual Evaluation Parameters\n\n')
    if Optimizer.calc_method=='VASP':
        Optimizer.output.write('\n----VASP INPUT------\n')
        Optimizer.output.write(Optimizer.vaspcalc + '\n')
    elif Optimizer.calc_method=='LAMMPS':
        Optimizer.output.write('\n----LAMMPS Input------\n')
        Optimizer.output.write('The potential used is '+repr(Optimizer.pot_file)+'\n')
        Optimizer.output.write('The potential style is '+repr(Optimizer.pair_style)+'\n')
    elif Optimizer.calc_method=='Mixed':
        Optimizer.output.write('\n----Mixed Run------\n')
        if Optimizer.pot_file != None:
            Optimizer.output.write('The potential used is '+repr(Optimizer.pot_file)+'\n')
        if Optimizer.pair_style !=None:
            Optimizer.output.write('The potential style is ' + repr(Optimizer.pair_style) + '\n')
        else:
            Optimizer.output.write('A Leonard Jones Potential is assumed for LAMMPS\n')
        Optimizer.output.write('VASP calculations are performed with: \n' + Optimizer.vaspcalc + '\n')
    if Optimizer.ase_min==True:
        Optimizer.output.write('Using ASE BFGS Minimizer \n')
        Optimizer.output.write('fmax = ' + repr(Optimizer.ase_min_fmax) + '\n')
        Optimizer.output.write('Solver_Max_Steps = ' + repr(Optimizer.ase_min_maxsteps) + '\n')
    elif Optimizer.lammps_min != None:
        Optimizer.output.write('Lammps minimizer used ' + repr(Optimizer.lammps_min) + repr(Optimizer.lammps_min_style) + '\n')
    elif Optimizer.calc_method=='VASP':
        Optimizer.output.write('See VASP input to see local minimizer options\n')
    else:
        Optimizer.output.write('No Local Minimizer Used')
    if Optimizer.structure =='Cluster':
        Optimizer.output.write('LAMMPS box size for periodic evaluation = '+repr(Optimizer.large_box_size)+'\n')

    Optimizer.output.write(section_break+'Selection Parameters\n\n')
    Optimizer.output.write('Fitness Scheme (fitness_scheme) : ' + repr(Optimizer.fitness_scheme) + '\n')
    if 'stem' in Optimizer.fitness_scheme:
        Optimizer.output.write('Parameters for Convoluation AutoSTEM calculation:\n'+repr(Optimizer.stemcalc.parameters)+'\n')
        Optimizer.output.write('Keep STEM simulated images = ' + repr(Optimizer.stemcalc.keep_files) + '\n')
        Optimizer.output.write('STEM alpha coefficient = ' + repr(Optimizer.stem_coeff) + '\n')
    Optimizer.output.write('Natural Selection Scheme (natural_selection_scheme): ' + repr(Optimizer.natural_selection_scheme) + '\n')
    if 'tournament' in Optimizer.natural_selection_scheme: Optimizer.output.write('    Tournament Size (tournsize) = ' + repr(Optimizer.tournsize) + '\n')
    if 'fuss' in Optimizer.natural_selection_scheme:
        if Optimizer.natural_selection_scheme != 'fuss1': 
            Optimizer.output.write('FUSS limit (fusslim) : ' + repr(Optimizer.fusslimit) + ' eV \n')
        if Optimizer.natural_selection_scheme == 'fussf':
            Optimizer.output.write('Fingerprint bin size (fpbin) : ' + repr(Optimizer.fpbin) + '\n')
            Optimizer.output.write('Fingerprint cutoff distance (fpcutoff) : ' + repr(Optimizer.fpcutoff) + '\n')
    if 'metropolis' in Optimizer.natural_selection_scheme: Optimizer.output.write('Starting Temperature for Metropolis Selection = '+repr(Optimizer.metropolis_temp)+'kT\n')

    Optimizer.output.write(section_break+'Convergence Parameters\n\n')
    Optimizer.output.write('Population convergence scheme (CONVERGENCE_SCHEME) : ' + repr(Optimizer.convergence_scheme) + '\n')
    Optimizer.output.write('Maximum number of generation (maxgen) = ' + repr(Optimizer.maxgen) + '\n')
    if 'rep' in Optimizer.convergence_scheme:
        Optimizer.output.write('Required generation of repetition of energy for convergence (reqrep) = ' + repr(Optimizer.reqrep) + '\n')
        Optimizer.output.write('Energy cutoff for repetition (tolerance) : ' + repr(Optimizer.tolerance) + '\n')
    Optimizer.output.write('Duplicate structure convergence control scheme (predator) '+ repr(Optimizer.predator) + '\n')
    Optimizer.output.write('Minimum energy difference for duplicate consideration (demin) : ' + repr(Optimizer.demin) + '\n')
    if Optimizer.predator=='adapting':
        Optimizer.output.write('Fraction of repeated generations at which to begin adaptation (adaptbegin) : ' + repr(Optimizer.adaptbegin) + '\n')
        Optimizer.output.write('Multiplier for adaptation exponential (adaptmultiplier) : ' + repr(Optimizer.adaptmultiplier) + '\n')

    Optimizer.output.write(section_break+'Output and Post-Processing Parameters\n\n')
    Optimizer.output.write('Structure atoms filename (filename) : ' + os.getcwd() + '/' + Optimizer.filename + '\n')
    Optimizer.output.write('Summary filename : ' + repr(Optimizer.summary.name) + '\n')
    Optimizer.output.write('Format for summary file (output_format) : ' + repr(Optimizer.output_format) + '\n')
    if Optimizer.genealogy: Optimizer.output.write('Genealogy File name : ' + repr(Optimizer.Genealogyfile.name) + '\n')
    if Optimizer.allenergyfile: Optimizer.output.write('All Energy file written \n')
    if Optimizer.lattice_concentration: Optimizer.output.write('Lattice Concentration File written \n')
    if Optimizer.indiv_defect_write: Optimizer.output.write('Individual cluster files will also be written \n')
    Optimizer.output.write('Vacancies output to final structures : ' + repr(Optimizer.vacancy_output) + '\n')
    Optimizer.summary.write(Optimizer.filename+'\n')

    Optimizer.output.flush()
    return