from MAST.utility import InputOptions
import numpy as np
#MAST INPUT OPTIONS
inputoptions = InputOptions()
 
###################################
#Section: ingredients
####################################
 
ingredients = dict()
ingredients['neblowmesh'] = dict()
ingredients['neblowmesh']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['neblowmesh']['mast_multiplyencut'] = '1.25'
ingredients['neblowmesh']['mast_nodes'] = '1'
ingredients['neblowmesh']['lclimb'] = 'true'
ingredients['neblowmesh']['mast_xc'] = 'pw91'
ingredients['neblowmesh']['isif'] = '2'
ingredients['neblowmesh']['spring'] = '-5'
ingredients['neblowmesh']['lcharg'] = 'false'
ingredients['neblowmesh']['prec'] = 'accurate'
ingredients['neblowmesh']['images'] = '1'
ingredients['neblowmesh']['ibrion'] = '1'
ingredients['neblowmesh']['mast_ppn'] = '1'
ingredients['neblowmesh']['mast_queue'] = 'default'
ingredients['neblowmesh']['mast_kpoints'] = list()
ingredients['neblowmesh']['mast_kpoints'].append(2)
ingredients['neblowmesh']['mast_kpoints'].append(2)
ingredients['neblowmesh']['mast_kpoints'].append(2)
ingredients['neblowmesh']['mast_kpoints'].append('m')
ingredients['neblowmesh']['ismear'] = '1'
ingredients['neblowmesh']['lwave'] = 'false'
ingredients['neblowmesh']['sigma'] = '0.2'
ingredients['neblowmesh']['potim'] = '0.5'
ingredients['neblowmesh']['nsw'] = '191'
ingredients['global'] = dict()
ingredients['global']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['global']['ibrion'] = '2'
ingredients['global']['mast_nodes'] = '1'
ingredients['global']['mast_xc'] = 'pw91'
ingredients['global']['isif'] = '2'
ingredients['global']['mast_queue'] = 'default'
ingredients['global']['lcharg'] = 'false'
ingredients['global']['prec'] = 'accurate'
ingredients['global']['mast_kpoints'] = list()
ingredients['global']['mast_kpoints'].append(1)
ingredients['global']['mast_kpoints'].append(1)
ingredients['global']['mast_kpoints'].append(1)
ingredients['global']['mast_kpoints'].append('g')
ingredients['global']['mast_multiplyencut'] = '1.25'
ingredients['global']['mast_ppn'] = '1'
ingredients['global']['ismear'] = '1'
ingredients['global']['lwave'] = 'false'
ingredients['global']['sigma'] = '0.2'
ingredients['global']['nsw'] = '191'
ingredients['neb'] = dict()
ingredients['neb']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['neb']['mast_multiplyencut'] = '1.25'
ingredients['neb']['mast_nodes'] = '1'
ingredients['neb']['lclimb'] = 'true'
ingredients['neb']['mast_xc'] = 'pw91'
ingredients['neb']['lcharg'] = 'false'
ingredients['neb']['spring'] = '-5'
ingredients['neb']['mast_kpoints'] = list()
ingredients['neb']['mast_kpoints'].append(1)
ingredients['neb']['mast_kpoints'].append(1)
ingredients['neb']['mast_kpoints'].append(1)
ingredients['neb']['mast_kpoints'].append('g')
ingredients['neb']['prec'] = 'accurate'
ingredients['neb']['isif'] = '2'
ingredients['neb']['ibrion'] = '1'
ingredients['neb']['mast_ppn'] = '1'
ingredients['neb']['mast_queue'] = 'default'
ingredients['neb']['images'] = '1'
ingredients['neb']['ismear'] = '1'
ingredients['neb']['lwave'] = 'false'
ingredients['neb']['sigma'] = '0.2'
ingredients['neb']['potim'] = '0.5'
ingredients['neb']['nsw'] = '191'
ingredients['sortnebstatic'] = dict()
ingredients['sortnebstatic']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['sortnebstatic']['mast_multiplyencut'] = '1.25'
ingredients['sortnebstatic']['mast_nodes'] = '1'
ingredients['sortnebstatic']['mast_xc'] = 'pw91'
ingredients['sortnebstatic']['lcharg'] = 'false'
ingredients['sortnebstatic']['mast_queue'] = 'default'
ingredients['sortnebstatic']['prec'] = 'accurate'
ingredients['sortnebstatic']['isif'] = '2'
ingredients['sortnebstatic']['ibrion'] = '2'
ingredients['sortnebstatic']['mast_ppn'] = '1'
ingredients['sortnebstatic']['mast_kpoints'] = list()
ingredients['sortnebstatic']['mast_kpoints'].append(1)
ingredients['sortnebstatic']['mast_kpoints'].append(1)
ingredients['sortnebstatic']['mast_kpoints'].append(1)
ingredients['sortnebstatic']['mast_kpoints'].append('g')
ingredients['sortnebstatic']['ismear'] = '1'
ingredients['sortnebstatic']['lwave'] = 'false'
ingredients['sortnebstatic']['sigma'] = '0.2'
ingredients['sortnebstatic']['nsw'] = '191'
ingredients['static'] = dict()
ingredients['static']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['static']['ibrion'] = '-1'
ingredients['static']['mast_nodes'] = '1'
ingredients['static']['mast_xc'] = 'pw91'
ingredients['static']['lcharg'] = 'false'
ingredients['static']['mast_queue'] = 'default'
ingredients['static']['prec'] = 'accurate'
ingredients['static']['isif'] = '2'
ingredients['static']['mast_multiplyencut'] = '1'
ingredients['static']['mast_ppn'] = '1'
ingredients['static']['mast_kpoints'] = list()
ingredients['static']['mast_kpoints'].append(1)
ingredients['static']['mast_kpoints'].append(1)
ingredients['static']['mast_kpoints'].append(1)
ingredients['static']['mast_kpoints'].append('g')
ingredients['static']['ismear'] = '1'
ingredients['static']['lwave'] = 'false'
ingredients['static']['sigma'] = '0.2'
ingredients['static']['nsw'] = '0'
ingredients['optimizelowmeshdefect'] = dict()
ingredients['optimizelowmeshdefect']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['optimizelowmeshdefect']['mast_multiplyencut'] = '1.25'
ingredients['optimizelowmeshdefect']['mast_nodes'] = '1'
ingredients['optimizelowmeshdefect']['mast_xc'] = 'pw91'
ingredients['optimizelowmeshdefect']['lcharg'] = 'false'
ingredients['optimizelowmeshdefect']['mast_queue'] = 'default'
ingredients['optimizelowmeshdefect']['prec'] = 'accurate'
ingredients['optimizelowmeshdefect']['isif'] = '2'
ingredients['optimizelowmeshdefect']['ibrion'] = '2'
ingredients['optimizelowmeshdefect']['mast_ppn'] = '1'
ingredients['optimizelowmeshdefect']['mast_kpoints'] = list()
ingredients['optimizelowmeshdefect']['mast_kpoints'].append(1)
ingredients['optimizelowmeshdefect']['mast_kpoints'].append(1)
ingredients['optimizelowmeshdefect']['mast_kpoints'].append(1)
ingredients['optimizelowmeshdefect']['mast_kpoints'].append('g')
ingredients['optimizelowmeshdefect']['ismear'] = '1'
ingredients['optimizelowmeshdefect']['lwave'] = 'false'
ingredients['optimizelowmeshdefect']['sigma'] = '0.2'
ingredients['optimizelowmeshdefect']['nsw'] = '191'
ingredients['optimizelowmesh'] = dict()
ingredients['optimizelowmesh']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['optimizelowmesh']['mast_multiplyencut'] = '1.25'
ingredients['optimizelowmesh']['mast_nodes'] = '1'
ingredients['optimizelowmesh']['mast_xc'] = 'pw91'
ingredients['optimizelowmesh']['lcharg'] = 'false'
ingredients['optimizelowmesh']['mast_queue'] = 'default'
ingredients['optimizelowmesh']['prec'] = 'accurate'
ingredients['optimizelowmesh']['isif'] = '3'
ingredients['optimizelowmesh']['ibrion'] = '2'
ingredients['optimizelowmesh']['mast_ppn'] = '1'
ingredients['optimizelowmesh']['mast_kpoints'] = list()
ingredients['optimizelowmesh']['mast_kpoints'].append(1)
ingredients['optimizelowmesh']['mast_kpoints'].append(1)
ingredients['optimizelowmesh']['mast_kpoints'].append(1)
ingredients['optimizelowmesh']['mast_kpoints'].append('g')
ingredients['optimizelowmesh']['ismear'] = '1'
ingredients['optimizelowmesh']['lwave'] = 'false'
ingredients['optimizelowmesh']['sigma'] = '0.2'
ingredients['optimizelowmesh']['nsw'] = '191'
ingredients['staticdefect'] = dict()
ingredients['staticdefect']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['staticdefect']['ibrion'] = '-1'
ingredients['staticdefect']['mast_nodes'] = '1'
ingredients['staticdefect']['mast_xc'] = 'pw91'
ingredients['staticdefect']['lcharg'] = 'false'
ingredients['staticdefect']['mast_queue'] = 'default'
ingredients['staticdefect']['prec'] = 'accurate'
ingredients['staticdefect']['isif'] = '2'
ingredients['staticdefect']['mast_multiplyencut'] = '1'
ingredients['staticdefect']['mast_ppn'] = '1'
ingredients['staticdefect']['mast_kpoints'] = list()
ingredients['staticdefect']['mast_kpoints'].append(1)
ingredients['staticdefect']['mast_kpoints'].append(1)
ingredients['staticdefect']['mast_kpoints'].append(1)
ingredients['staticdefect']['mast_kpoints'].append('g')
ingredients['staticdefect']['ismear'] = '1'
ingredients['staticdefect']['lwave'] = 'false'
ingredients['staticdefect']['sigma'] = '0.2'
ingredients['staticdefect']['nsw'] = '0'
ingredients['optimizedefect'] = dict()
ingredients['optimizedefect']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['optimizedefect']['mast_multiplyencut'] = '1.25'
ingredients['optimizedefect']['mast_nodes'] = '1'
ingredients['optimizedefect']['mast_xc'] = 'pw91'
ingredients['optimizedefect']['lcharg'] = 'false'
ingredients['optimizedefect']['mast_queue'] = 'default'
ingredients['optimizedefect']['prec'] = 'accurate'
ingredients['optimizedefect']['isif'] = '2'
ingredients['optimizedefect']['ibrion'] = '2'
ingredients['optimizedefect']['mast_ppn'] = '1'
ingredients['optimizedefect']['mast_kpoints'] = list()
ingredients['optimizedefect']['mast_kpoints'].append(1)
ingredients['optimizedefect']['mast_kpoints'].append(1)
ingredients['optimizedefect']['mast_kpoints'].append(1)
ingredients['optimizedefect']['mast_kpoints'].append('g')
ingredients['optimizedefect']['ismear'] = '1'
ingredients['optimizedefect']['lwave'] = 'false'
ingredients['optimizedefect']['sigma'] = '0.2'
ingredients['optimizedefect']['nsw'] = '191'
ingredients['optimize'] = dict()
ingredients['optimize']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['optimize']['mast_multiplyencut'] = '1.25'
ingredients['optimize']['mast_nodes'] = '1'
ingredients['optimize']['mast_xc'] = 'pw91'
ingredients['optimize']['lcharg'] = 'false'
ingredients['optimize']['mast_queue'] = 'default'
ingredients['optimize']['prec'] = 'accurate'
ingredients['optimize']['isif'] = '3'
ingredients['optimize']['ibrion'] = '2'
ingredients['optimize']['mast_ppn'] = '1'
ingredients['optimize']['mast_kpoints'] = list()
ingredients['optimize']['mast_kpoints'].append(1)
ingredients['optimize']['mast_kpoints'].append(1)
ingredients['optimize']['mast_kpoints'].append(1)
ingredients['optimize']['mast_kpoints'].append('g')
ingredients['optimize']['ismear'] = '1'
ingredients['optimize']['lwave'] = 'false'
ingredients['optimize']['sigma'] = '0.2'
ingredients['optimize']['nsw'] = '191'
ingredients['inducedefect'] = dict()
ingredients['inducedefect']['mast_exec'] = '//share/apps/vasp5.2_cneb'
ingredients['inducedefect']['mast_multiplyencut'] = '1.25'
ingredients['inducedefect']['mast_nodes'] = '1'
ingredients['inducedefect']['mast_xc'] = 'pw91'
ingredients['inducedefect']['lcharg'] = 'false'
ingredients['inducedefect']['mast_queue'] = 'default'
ingredients['inducedefect']['prec'] = 'accurate'
ingredients['inducedefect']['isif'] = '2'
ingredients['inducedefect']['ibrion'] = '2'
ingredients['inducedefect']['mast_ppn'] = '1'
ingredients['inducedefect']['mast_kpoints'] = list()
ingredients['inducedefect']['mast_kpoints'].append(1)
ingredients['inducedefect']['mast_kpoints'].append(1)
ingredients['inducedefect']['mast_kpoints'].append(1)
ingredients['inducedefect']['mast_kpoints'].append('g')
ingredients['inducedefect']['ismear'] = '1'
ingredients['inducedefect']['lwave'] = 'false'
ingredients['inducedefect']['sigma'] = '0.2'
ingredients['inducedefect']['nsw'] = '191'
inputoptions.options['ingredients'] = ingredients
 
###################################
#Section: neb
####################################
 
neb = dict()
neb['images'] = 1
neb['neblines'] = dict()
neb['neblines']['vac1-vac2'] = list()
neb['neblines']['vac1-vac2'].append(['fe', ' 0 0 0', ' 0 .5 0.5'])
inputoptions.options['neb'] = neb
 
###################################
#Section: recipe
####################################
 
recipe = dict()
recipe['recipe_file'] = '/home/tam/bin/git/MAST4pymatgen/recipe_templates/neb.txt'
inputoptions.options['recipe'] = recipe
 
###################################
#Section: mast
####################################
 
mast = dict()
mast['system_name'] = 'smallfcc'
mast['scratch_directory'] = '//home/tam/bin/git/MAST4pymatgen/scratch2'
mast['program'] = 'vasp'
mast['input_stem'] = '//home/tam/bin/git/MAST4pymatgen/scratch2/neb_test_noint_20130620T092336_'
inputoptions.options['mast'] = mast
 
###################################
#Section: defects
####################################
 
defects = dict()
defects['coord_type'] = 'fractional'
defects['num_defects'] = 2
defects['defects'] = dict()
defects['defects']['defect_vac2'] = dict()
defects['defects']['defect_vac2']['subdefect_1'] = dict()
defects['defects']['defect_vac2']['subdefect_1']['symbol'] = 'fe'
defects['defects']['defect_vac2']['subdefect_1']['type'] = 'vacancy'
defects['defects']['defect_vac2']['subdefect_1']['coordinates']=np.array([0.0, 0.5, 0.5])
defects['defects']['defect_vac1'] = dict()
defects['defects']['defect_vac1']['subdefect_1'] = dict()
defects['defects']['defect_vac1']['subdefect_1']['symbol'] = 'fe'
defects['defects']['defect_vac1']['subdefect_1']['type'] = 'vacancy'
defects['defects']['defect_vac1']['subdefect_1']['coordinates']=np.array([0.0, 0.0, 0.0])
inputoptions.options['defects'] = defects
 
###################################
#Section: structure
####################################
 
structure = dict()
structure['primitive'] = False
structure['atom_list'] = None
structure['coordinates'] = None
structure['lattice'] = None
structure['posfile'] = 'smallfcc_poscar'
structure['coord_type'] = 'cartesian'
structure['symmetry_only'] = False
structure['spacegroup'] = None
#ignoring created structure
inputoptions.options['structure'] = structure
 
###############################################
#Structure Creation Section
################################################
 
from pymatgen.io.vaspio import Poscar
structure = Poscar.from_file('smallfcc_poscar').structure
inputoptions.set_item('structure','structure',structure)
 
###############################################
#Buffet Command Section
################################################
 
from MAST.recipe.recipeinput import RecipeInput
from MAST.buffet.buffetmanager import Buffet
from MAST.mast import MAST
mast_obj = MAST()
mast_obj.start_from_input_options(inputoptions)
