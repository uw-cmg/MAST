from MAST.utility import InputOptions
import numpy as np
#MAST INPUT OPTIONS
inputoptions = InputOptions()
 
###################################
#Section: ingredients
####################################
 
ingredients = dict()
ingredients['neb_to_neb'] = dict()
ingredients['neb_to_neb']['mast_nodes'] = '1'
ingredients['neb_to_neb']['mast_write_method'] = 'write_neb'
ingredients['neb_to_neb']['lclimb'] = 'True'
ingredients['neb_to_neb']['mast_update_children_method'] = 'give_neb_structures_to_neb'
ingredients['neb_to_neb']['spring'] = '-5'
ingredients['neb_to_neb']['images'] = '1'
ingredients['neb_to_neb']['ibrion'] = '1'
ingredients['neb_to_neb']['mast_kpoints'] = list()
ingredients['neb_to_neb']['mast_kpoints'].append(1)
ingredients['neb_to_neb']['mast_kpoints'].append(1)
ingredients['neb_to_neb']['mast_kpoints'].append(1)
ingredients['neb_to_neb']['mast_kpoints'].append('G')
ingredients['neb_to_neb']['potim'] = '0.5'
ingredients['phononparse'] = dict()
ingredients['phononparse']['lnosym'] = '.True.'
ingredients['phononparse']['mast_exec'] = '$MAST_INSTALL_PATH/bin/phon_henry'
ingredients['phononparse']['ldrift'] = '.False.'
ingredients['phononparse']['temperature'] = '273'
ingredients['phononparse']['mast_multiplyencut'] = '1.25'
ingredients['phononparse']['lfree'] = '.True.'
ingredients['phononparse']['nd'] = '3'
ingredients['phononparse']['qa'] = '11'
ingredients['phononparse']['qc'] = '11'
ingredients['phononparse']['qb'] = '11'
ingredients['phononparse']['mast_program'] = 'phon'
ingredients['phononparse']['ptemp'] = '10 110'
ingredients['phononparse']['lsuper'] = '.False.'
ingredients['phonon_to_phononparse'] = dict()
ingredients['phonon_to_phononparse']['mast_write_method'] = 'write_phonon_multiple'
ingredients['phonon_to_phononparse']['mast_complete_method'] = 'complete_subfolders'
ingredients['phonon_to_phononparse']['mast_ready_method'] = 'ready_subfolders'
ingredients['phonon_to_phononparse']['mast_run_method'] = 'run_subfolders'
ingredients['phonon_to_phononparse']['mast_update_children_method'] = 'give_phonon_multiple_forces_and_displacements'
ingredients['phonon_to_phononparse']['icharg'] = '1'
ingredients['phonon_to_phononparse']['ibrion'] = '5'
ingredients['phonon_to_phononparse']['nfree'] = '2'
ingredients['phonon_to_phononparse']['potim'] = '0.01'
ingredients['phonon_to_phononparse']['istart'] = '1'
ingredients['singlerun_to_neb'] = dict()
ingredients['singlerun_to_neb']['mast_update_children_method'] = 'give_structure_and_energy_to_neb'
ingredients['singlerun_to_neb']['lwave'] = 'True'
ingredients['singlerun_to_neb']['ibrion'] = '-1'
ingredients['singlerun_to_neb']['lcharge'] = 'True'
ingredients['singlerun_to_neb']['nsw'] = '0'
ingredients['nebstat_to_nebphonon'] = dict()
ingredients['nebstat_to_nebphonon']['mast_write_method'] = 'write_neb_subfolders'
ingredients['nebstat_to_nebphonon']['mast_complete_method'] = 'complete_neb_subfolders'
ingredients['nebstat_to_nebphonon']['mast_ready_method'] = 'ready_neb_subfolders'
ingredients['nebstat_to_nebphonon']['mast_run_method'] = 'run_neb_subfolders'
ingredients['nebstat_to_nebphonon']['mast_update_children_method'] = 'give_saddle_structure'
ingredients['nebstat_to_nebphonon']['ibrion'] = '-1'
ingredients['nebstat_to_nebphonon']['nsw'] = '0'
ingredients['global'] = dict()
ingredients['global']['mast_exec'] = '//share/apps/vasp5.2_cNEB'
ingredients['global']['ibrion'] = '2'
ingredients['global']['mast_nodes'] = '1'
ingredients['global']['mast_write_method'] = 'write_singlerun'
ingredients['global']['mast_complete_method'] = 'complete_singlerun'
ingredients['global']['mast_xc'] = 'PW91'
ingredients['global']['isif'] = '2'
ingredients['global']['mast_ready_method'] = 'ready_singlerun'
ingredients['global']['mast_queue'] = 'default'
ingredients['global']['lcharg'] = 'False'
ingredients['global']['mast_run_method'] = 'run_singlerun'
ingredients['global']['prec'] = 'Accurate'
ingredients['global']['mast_kpoints'] = list()
ingredients['global']['mast_kpoints'].append(2)
ingredients['global']['mast_kpoints'].append(2)
ingredients['global']['mast_kpoints'].append(2)
ingredients['global']['mast_kpoints'].append('M')
ingredients['global']['mast_update_children_method'] = 'give_structure'
ingredients['global']['mast_multiplyencut'] = '1.5'
ingredients['global']['mast_program'] = 'vasp'
ingredients['global']['mast_ppn'] = '1'
ingredients['global']['ismear'] = '1'
ingredients['global']['lwave'] = 'False'
ingredients['global']['sigma'] = '0.2'
ingredients['global']['nsw'] = '191'
ingredients['volrelax_to_singlerun'] = dict()
ingredients['volrelax_to_singlerun']['isif'] = '3'
ingredients['neb_to_nebstat'] = dict()
ingredients['neb_to_nebstat']['mast_nodes'] = '1'
ingredients['neb_to_nebstat']['mast_write_method'] = 'write_neb'
ingredients['neb_to_nebstat']['lclimb'] = 'True'
ingredients['neb_to_nebstat']['mast_update_children_method'] = 'give_neb_structures_to_neb'
ingredients['neb_to_nebstat']['spring'] = '-5'
ingredients['neb_to_nebstat']['ibrion'] = '1'
ingredients['neb_to_nebstat']['images'] = '1'
ingredients['neb_to_nebstat']['potim'] = '0.5'
ingredients['induce_defect'] = dict()
ingredients['induce_defect']['mast_ready_method'] = 'ready_defect'
ingredients['induce_defect']['mast_run_method'] = 'run_defect'
ingredients['induce_defect']['mast_write_method'] = 'no_setup'
ingredients['induce_defect']['mast_complete_method'] = 'complete_structure'
ingredients['singlerun_to_phonon'] = dict()
ingredients['singlerun_to_phonon']['mast_multiplyencut'] = '1.25'
ingredients['singlerun_to_phonon']['lcharge'] = 'True'
ingredients['singlerun_to_phonon']['mast_update_children_method'] = 'give_structure_and_restart_files'
ingredients['singlerun_to_phonon']['ibrion'] = '-1'
ingredients['singlerun_to_phonon']['lwave'] = 'True'
ingredients['singlerun_to_phonon']['nsw'] = '0'
inputoptions.options['ingredients'] = ingredients
 
###################################
#Section: recipe
####################################
 
recipe = dict()
recipe['recipe_file'] = '/home/tam/tammast/recipe_templates/treetest.txt'
inputoptions.options['recipe'] = recipe
 
###################################
#Section: neb
####################################
 
neb = dict()
neb['images'] = 1
neb['neblines'] = dict()
neb['neblines']['vac1-vac2'] = list()
neb['neblines']['vac1-vac2'].append(['Al', ' 0 0 0', ' 0.5 0.5 0.0'])
inputoptions.options['neb'] = neb
 
###################################
#Section: phonon
####################################
 
phonon = dict()
phonon['phonon'] = dict()
phonon['phonon']['perfect'] = dict()
phonon['phonon']['perfect']['phonon_center_site'] = '0.5 0.5 0'
phonon['phonon']['perfect']['phonon_center_radius'] = '1'
phonon['phonon']['vac2'] = dict()
phonon['phonon']['vac2']['phonon_center_site'] = '0.0 0.0 0'
phonon['phonon']['vac2']['phonon_center_radius'] = '1'
phonon['phonon']['vac1'] = dict()
phonon['phonon']['vac1']['phonon_center_site'] = '0.5 0.5 0'
phonon['phonon']['vac1']['phonon_center_radius'] = '1'
phonon['phonon']['vac1-vac2'] = dict()
phonon['phonon']['vac1-vac2']['phonon_center_site'] = '0.25 0.25 0'
phonon['phonon']['vac1-vac2']['phonon_center_radius'] = '1'
inputoptions.options['phonon'] = phonon
 
###################################
#Section: mast
####################################
 
mast = dict()
mast['system_name'] = 'smalldemo_Al'
mast['scratch_directory'] = '/home/tam/tammast/SCRATCH'
mast['timestamp'] = 'Wed Aug 21 12:45:39 2013'
mast['input_stem'] = '/home/tam/tammast/SCRATCH/test_small_new_template_20130821T124539_'
mast['program'] = 'vasp'
mast['working_directory'] = '/home/tam/tammast/SCRATCH/smalldemo_Al_Al_None_test_small_new_template_20130821T124539'
inputoptions.options['mast'] = mast
 
###################################
#Section: defects
####################################
 
defects = dict()
defects['coord_type'] = 'fractional'
defects['num_defects'] = 2
defects['defects'] = dict()
defects['defects']['vac2'] = dict()
defects['defects']['vac2']['threshold'] = 0.0001
defects['defects']['vac2']['charge'] = list()
defects['defects']['vac2']['charge'].append(0)
defects['defects']['vac2']['coord_type'] = 'fractional'
defects['defects']['vac2']['subdefect_1'] = dict()
defects['defects']['vac2']['subdefect_1']['symbol'] = 'Al'
defects['defects']['vac2']['subdefect_1']['type'] = 'vacancy'
defects['defects']['vac2']['subdefect_1']['coordinates']=np.array([0.5, 0.5, 0.0])
defects['defects']['vac1'] = dict()
defects['defects']['vac1']['threshold'] = 0.0001
defects['defects']['vac1']['charge'] = list()
defects['defects']['vac1']['charge'].append(0)
defects['defects']['vac1']['coord_type'] = 'fractional'
defects['defects']['vac1']['subdefect_1'] = dict()
defects['defects']['vac1']['subdefect_1']['symbol'] = 'Al'
defects['defects']['vac1']['subdefect_1']['type'] = 'vacancy'
defects['defects']['vac1']['subdefect_1']['coordinates']=np.array([0.0, 0.0, 0.0])
inputoptions.options['defects'] = defects
 
###################################
#Section: structure
####################################
 
structure = dict()
structure['element_map'] = dict()
structure['element_map']['X1'] = 'Al'
structure['primitive'] = False
structure['atom_list'] = list()
structure['atom_list'].append('Al')
structure['atom_list'].append('Al')
structure['atom_list'].append('Al')
structure['atom_list'].append('Al')
structure['coordinates']=np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
structure['lattice']=np.array([[3.5, 0.0, 0.0], [0.0, 3.5, 0.0], [0.0, 0.0, 3.5]])
structure['posfile'] = None
structure['coord_type'] = 'fractional'
structure['symmetry_only'] = False
structure['spacegroup'] = None
#ignoring created structure
inputoptions.options['structure'] = structure
 
###############################################
#Structure Creation Section
################################################
 
iopscoords=inputoptions.get_item('structure','coordinates')
iopslatt=inputoptions.get_item('structure','lattice')
iopsatoms=inputoptions.get_item('structure','atom_list')
iopsctype=inputoptions.get_item('structure','coord_type')
from MAST.utility import MAST2Structure
structure = MAST2Structure(lattice=iopslatt,
    coordinates=iopscoords, atom_list=iopsatoms,
    coord_type=iopsctype)
inputoptions.update_item('structure','structure',structure)
 
###############################################
#MAST Command Section
################################################
 
from MAST.mast import MAST
mast_obj = MAST()
mast_obj.start_from_input_options(inputoptions)
