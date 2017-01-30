"""Tests for Inputparser"""

from MAST.parsers.inputparser import InputParser

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import InputOptions
from MAST.utility import MASTFile
from MAST.utility import MASTError
from pymatgen.io.vasp import Poscar
import numpy as np

testname="inputparser_test"
testdir = dirutil.get_test_dir(testname)

class TestInputparser(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        s_dict = dict({\
                'mast' : myip.parse_mast_section,
                'structure' : myip.parse_structure_section,
                'ingredients' : myip.parse_ingredients_section,
                'defects' : myip.parse_defects_section,
                'recipe' : myip.parse_recipe_section,
                'neb' : myip.parse_neb_section,
                'chemical_potentials' : myip.parse_chemical_potentials_section,
                'summary' : myip.parse_summary_section,
                'personal_recipe' : myip.parse_personal_recipe_section,
                'scaling' : myip.parse_scaling_section,
                               })
        self.assertItemsEqual(myip.section_parsers, s_dict)
        #self.testclass.__init__(**kwargs)

    def test_parse(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        myoptions = myip.parse()
        print myoptions.options.keys()
        self.assertItemsEqual(myoptions.options.keys(),['mast','structure','defects','recipe','ingredients','neb'])
        #ofile = MASTFile()
        #ofile.data = str(myoptions)
        #ofile.to_file(os.path.join(os.getcwd(),"long_input_options_output.txt"))
        #self.testclass.parse()
        
    def test_parse_defect_formation(self):
        myip = InputParser(inputfile="defect_formation_energy.inp")
        myoptions = myip.parse()
        print myoptions.options.keys()
        self.assertItemsEqual(myoptions.options.keys(),['mast','structure','defects','recipe','ingredients','chemical_potentials','summary'])
        #ofile = MASTFile()
        #ofile.data = str(myoptions)
        #ofile.to_file(os.path.join(os.getcwd(),"long_input_options_output.txt"))
        #self.testclass.parse()

    def test_parse_mast_section(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/mast_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_mast_section('mast',cleanlines,myoptions)
        mdict=dict()
        mdict['mast']=dict()
        mdict['mast']['system_name'] = "SystemName!"
        self.assertItemsEqual(myoptions.options.keys(), ['mast'])
        self.assertItemsEqual(myoptions.options['mast'].keys(), ['system_name','mast_auto_correct'])
        self.assertEqual(myoptions.options['mast']['system_name'],mdict['mast']['system_name']) 
        self.assertEqual(myoptions.get_item('mast', 'mast_auto_correct'), "True")
        #print myoptions.options['mast']
        #self.testclass.parse_mast_section(section_name, section_content, options)

    def test_parse_structure_section_noposfile(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/structure_noposfile_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_structure_section('structure',cleanlines,myoptions)
        #print myoptions
        self.assertEqual(myoptions.options['structure']['atom_list'],['Al','Mg','Al','Fe'])
        self.assertEqual(myoptions.options['structure']['coord_type'],'fractional')
        self.assertTrue(np.array_equal(myoptions.options['structure']['lattice'],np.array([[3.5,1,2],[3,3.5,4],[5,6,3.5]])))
        self.assertTrue(np.array_equal(myoptions.options['structure']['coordinates'],np.array([[0,0,0],[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]])))
        self.assertItemsEqual(myoptions.options['structure']['element_map'],dict({'X1':'Al','X2':'Mg','X3':'Fe'}))

        #mdict=dict()
        #mdict['structure'] = dict()
        #mdict['structure']['posfile'] = os.path.join(testdir, "POSCAR_startpos")
        #self.assertItemsEqual(myoptions.options.keys(), ['structure'])
        #self.assertEqual(myoptions.options['structure']['posfile'],mdict['structure']['posfile']) 
        #self.testclass.parse_structure_section(section_name, section_content, options)
    
    def test_parse_structure_section_posfile(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/structure_posfile_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_structure_section('structure',cleanlines,myoptions)
        #print myoptions
        mdict=dict()
        mdict['structure'] = dict()
        mdict['structure']['posfile']="POSCAR_startpos"
        #mdict['structure']['posfile'] = os.path.join(testdir, "POSCAR_startpos")
        self.assertItemsEqual(myoptions.options.keys(), ['structure'])
        self.assertEqual(myoptions.options['structure']['posfile'],mdict['structure']['posfile']) 
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/structure_posfile_lines_nosuch.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        self.assertRaises(MASTError, myip.parse_structure_section,'structure',cleanlines,myoptions)
        #self.testclass.parse_structure_section(section_name, section_content, options)

    def test_parse_defects_section(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/defects_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_defects_section('defects',cleanlines,myoptions)
        print myoptions.options['defects']
        #Test coordinates first, to avoid dictionary-equal error on comparing arryas
        self.assertTrue(np.array_equal(myoptions.options['defects']['defects']['defect_3']['subdefect_1']['coordinates'],np.array([0.,0.,0.],'float')))
        self.assertTrue(np.array_equal(myoptions.options['defects']['defects']['Vac@Al-Sub@Fe']['subdefect_1']['coordinates'],np.array([0.,0.,0.],'float')))
        self.assertTrue(np.array_equal(myoptions.options['defects']['defects']['Vac@Al-Sub@Fe']['subdefect_2']['coordinates'],np.array([0.5,0.5,0.],'float')))
        self.assertTrue(np.array_equal(myoptions.options['defects']['defects']['AntiGe@Fe']['subdefect_1']['coordinates'],np.array([0.1,0.3,0.4],'float')))
        self.assertTrue(np.array_equal(myoptions.options['defects']['defects']['MgInt']['subdefect_1']['coordinates'],np.array([0.,0.2,0.3],'float')))
        self.assertItemsEqual(myoptions.options['defects']['defects']['Vac@Al-Sub@Fe']['phonon'].keys(),['host1','sub2'])
        #Have to blank out the numpy arrays so that the rest of the diff'ing will work
        myoptions.options['defects']['defects']['defect_3']['subdefect_1']['coordinates']='blank'
        myoptions.options['defects']['defects']['Vac@Al-Sub@Fe']['subdefect_1']['coordinates']='blank'
        myoptions.options['defects']['defects']['Vac@Al-Sub@Fe']['subdefect_2']['coordinates']='blank'
        myoptions.options['defects']['defects']['AntiGe@Fe']['subdefect_1']['coordinates']='blank'
        myoptions.options['defects']['defects']['MgInt']['subdefect_1']['coordinates']='blank'
        mdict=dict()
        mdict['coord_type']='fractional'
        mdict['num_defects']=4
        mdict['defects']=dict()
        mdict['defects']['defect_3']=dict()
        mdict['defects']['defect_3']['threshold']=0.0001
        mdict['defects']['defect_3']['charge']=[0] #charge of 0
        mdict['defects']['defect_3']['coord_type']='fractional'
        mdict['defects']['defect_3']['subdefect_1']=dict()
        mdict['defects']['defect_3']['subdefect_1']['symbol']='Al'
        mdict['defects']['defect_3']['subdefect_1']['type']='vacancy'
        mdict['defects']['defect_3']['subdefect_1']['coordinates']='blank' #np.array([0.,0.,0.],'float')
        mdict['defects']['defect_3']['phonon']={}
        mdict['defects']['Vac@Al-Sub@Fe']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['threshold']=0.0001
        mdict['defects']['Vac@Al-Sub@Fe']['charge']=[-2,-1,0,1,2,3,4,5]
        mdict['defects']['Vac@Al-Sub@Fe']['coord_type']='fractional'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']['symbol']='Al'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']['type']='vacancy'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']['coordinates']='blank' #np.array([0.,0.,0.],'float')
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']['symbol']='Fe'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']['type']='substitution'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']['coordinates']='blank' #np.array([0.5,0.5,0.],'float')
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']={}
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['host1']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['host1']['phonon_center_site']='0.1 0.2 0.3'
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['host1']['phonon_center_radius']=3.0
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['host1']['threshold']=0.07
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['sub2']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['sub2']['phonon_center_site']='0.2 0.4 0.0'
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['sub2']['phonon_center_radius']=1.0
        mdict['defects']['Vac@Al-Sub@Fe']['phonon']['sub2']['threshold']=0.1
        mdict['defects']['AntiGe@Fe']=dict()
        mdict['defects']['AntiGe@Fe']['threshold']=0.0001
        mdict['defects']['AntiGe@Fe']['charge']=[3]
        mdict['defects']['AntiGe@Fe']['coord_type']='fractional'
        mdict['defects']['AntiGe@Fe']['subdefect_1']=dict()
        mdict['defects']['AntiGe@Fe']['subdefect_1']['symbol']='Ge'
        mdict['defects']['AntiGe@Fe']['subdefect_1']['type']='antisite'
        mdict['defects']['AntiGe@Fe']['subdefect_1']['coordinates']='blank' #np.array([0.1,0.3,0.4],'float')
        mdict['defects']['AntiGe@Fe']['phonon']={}
        mdict['defects']['MgInt']=dict()
        mdict['defects']['MgInt']['threshold']=0.0001
        mdict['defects']['MgInt']['charge']=[0,1,2,3]
        mdict['defects']['MgInt']['coord_type']='fractional'
        mdict['defects']['MgInt']['subdefect_1']=dict()
        mdict['defects']['MgInt']['subdefect_1']['symbol']='Mg'
        mdict['defects']['MgInt']['subdefect_1']['type']='interstitial'
        mdict['defects']['MgInt']['subdefect_1']['coordinates']='blank' #np.array([0.,0.2,0.3],'float')
        mdict['defects']['MgInt']['phonon']={}
        maxdiff=self.maxDiff
        self.maxDiff=None
        self.assertEqual(myoptions.options['defects'],mdict) 
        self.maxDiff=maxdiff
        #self.testclass.parse_defects_section(section_name, section_content, options)

    def test_parse_recipe_section(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        print testdir
        minput = MASTFile("%s/recipe_neb.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line)
        myoptions = InputOptions()
        myip.parse_recipe_section('recipe',cleanlines,myoptions)
        #print myoptions
        mdict=dict()
        #mdict['recipe_file'] = os.path.join(os.getenv("MAST_RECIPE_PATH"),'defects_test.txt')
        mdict['recipe_file'] = cleanlines
        self.assertEqual(myoptions.options['recipe'],mdict)
        #self.testclass.parse_recipe_section(section_name, section_content, options)

    def test_parse_ingredients_section(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/ingredient_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_ingredients_section('ingredients',cleanlines,myoptions)
        print myoptions
        mdict=dict()
        mdict['global']=dict()
        mdict['global']['mast_complete_method']='complete_singlerun'
        mdict['global']['mast_update_children_method']='give_structure'
        mdict['global']['mast_queue']='default'
        mdict['global']['mast_run_method']='run_singlerun'
        mdict['global']['lcharg']='False'
        mdict['global']['mast_multiplyencut']='1.5'
        mdict['global']['mast_ppn']='1'
        mdict['global']['ismear']='1'
        mdict['global']['nsw']='191'
        mdict['global']['mast_exec']='mpiexec //share/apps/vasp5.2_cNEB'
        mdict['global']['mast_nodes']='1'
        mdict['global']['mast_xc']='PW91'
        mdict['global']['prec']='Accurate'
        mdict['global']['mast_kpoints']=[2, 2, 2, 'M']
        mdict['global']['mast_write_method']='write_singlerun'
        mdict['global']['ldauu']='1 3 5'
        mdict['global']['ibrion']='2'
        mdict['global']['mast_ready_method']='ready_singlerun'
        mdict['global']['isif']='2'
        mdict['global']['mast_program']='vasp'
        mdict['global']['lwave']='False'
        mdict['global']['sigma']='0.2'
        mdict['phononingredient']=dict()
        mdict['phononingredient']['mast_exec']='phon several words WiTh CaPs'
        mdict['phononingredient']['ibrion']='5'
        mdict['phononingredient']['mast_program']='phon'
        maxdiff=self.maxDiff
        self.maxDiff=None
        self.assertEqual(myoptions.options['ingredients'],mdict)
        self.maxDiff=maxdiff
        #self.testclass.parse_ingredients_section(section_name, section_content, options)

    def test_parse_neb_section(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/neb_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_neb_section('neb',cleanlines,myoptions)
        print myoptions
        mdict=dict()
        mdict['nebs']=dict()
        mdict['nebs']['vac1-vac2']=dict()
        mdict['nebs']['vac1-vac2']['images']=3
        mdict['nebs']['vac1-vac2']['lines']=list()
        mdict['nebs']['vac1-vac2']['phonon']=dict()
        mdict['nebs']['vac1-vac2']['lines'].append(list(['Cr',' 0.5 0.5 0.0',' 0.0 0.0 0.0']))
        mdict['nebs']['vac1-vac2']['lines'].append(list(['Fe',' 0.2 0.1 0.3',' 0.2 0.1 0.4']))
        mdict['nebs']['vac1-vac3']=dict()
        mdict['nebs']['vac1-vac3']['images']=5
        mdict['nebs']['vac1-vac3']['lines']=list()
        mdict['nebs']['vac1-vac3']['lines'].append(list(['Al',' 0.1 0.1 0.0',' 0.3 0.4 0.1']))
        mdict['nebs']['vac1-vac3']['phonon']=dict()
        mdict['nebs']['vac1-vac3']['phonon']['solvent1']=dict()
        mdict['nebs']['vac1-vac3']['phonon']['solvent1']['phonon_center_site']='0.3 0.2 0.1'
        mdict['nebs']['vac1-vac3']['phonon']['solvent1']['phonon_center_radius']=3.0
        mdict['nebs']['vac1-vac3']['phonon']['solvent1']['threshold']=0.1
        mdict['nebs']['vac1-vac3']['phonon']['solute1']=dict()
        mdict['nebs']['vac1-vac3']['phonon']['solute1']['phonon_center_site']='0.2 0.1 0.5'
        mdict['nebs']['vac1-vac3']['phonon']['solute1']['phonon_center_radius']=0.1
        mdict['nebs']['vac1-vac3']['phonon']['solute1']['threshold']=0.05
        self.assertEqual(myoptions.options['neb'],mdict)
        #self.testclass.parse_neb_section(section_name, section_content, options)

    def test_parse_chemical_potentials_section(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/chemical_potential_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_chemical_potentials_section('chemical_potentials',cleanlines,myoptions)
        print myoptions
        mdict=dict()
        mdict['As rich']=dict()
        mdict['As rich']['As']=3.5
        mdict['As rich']['Ga']=4.5
        mdict['Ga rich']=dict()
        mdict['Ga rich']['As']=4.5
        mdict['Ga rich']['Ga']=3.5
        self.assertEqual(myoptions.options['chemical_potentials'],mdict)
        #self.testclass.parse_chemical_potentials_section(section_name, section_content, options)

    def test_parse_phonon_section(self):
        raise SkipTest
        #Phonon section in input file is now obsolete!
        myip = InputParser(inputfile="neb_with_phonons.inp")
        minput = MASTFile("%s/phonon_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_phonon_section('phonon',cleanlines,myoptions)
        print myoptions
        mdict=dict()
        mdict['perfect']=dict()
        mdict['perfect']['phonon_center_site']='0.5 0.5 0'
        mdict['perfect']['phonon_center_radius']='1'
        mdict['vac1']=dict()
        mdict['vac1']['phonon_center_site']='0.5 0.5 0'
        mdict['vac1']['phonon_center_radius']='1'
        mdict['vac2']=dict()
        mdict['vac2']['phonon_center_site']='0.0 0.0 0'
        mdict['vac2']['phonon_center_radius']='1'
        mdict['vac1-vac2']=dict()
        mdict['vac1-vac2']['phonon_center_site']='0.25 0.25 0'
        mdict['vac1-vac2']['phonon_center_radius']='1'
        self.assertEqual(myoptions.options['phonon']['phonon'],mdict)
        #self.testclass.parse_phonon_section(section_name, section_content, options)

    def test_perform_element_mapping(self):
        myip = InputParser(inputfile="element_mapping.inp")
        myoptions = myip.parse()
        self.assertItemsEqual(myoptions.options.keys(),['mast','structure','defects','recipe','ingredients','neb','chemical_potentials'])
        self.assertEqual(myoptions.options['defects']['defects']['group1']['subdefect_1']['symbol'],'Xe')
        self.assertEqual(myoptions.options['defects']['defects']['group1']['subdefect_2']['symbol'],'Ar')
        self.assertEqual(myoptions.options['defects']['defects']['group2']['subdefect_1']['symbol'],'Ar')
        self.assertEqual(myoptions.options['defects']['defects']['group2']['subdefect_2']['symbol'],'Xe')
        self.assertEqual(myoptions.options['neb']['nebs']['group1-group2']['lines'][0][0],'Ar')
        self.assertEqual(myoptions.options['neb']['nebs']['group1-group2']['lines'][1][0],'Xe')

        #self.testclass.perform_element_mapping(input_options)

    def test_validate_execs(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        myoptions = myip.parse()
        myoptions.options['ingredients']['global'].pop('mast_exec')
        myoptions.options['ingredients']['diffcoeff'].pop('mast_exec')
        self.assertRaises(MASTError,myip.validate_execs,myoptions)
        #self.testclass.validate_execs(input_options)

    def test_set_structure_from_inputs(self):
        myip = InputParser(inputfile="neb_with_phonons.inp")
        myoptions = myip.parse()
        compare_struct = Poscar.from_file("nebphononsposcar").structure
        self.assertEqual(myoptions.options['structure']['structure'],compare_struct)
        #self.testclass.set_structure_from_inputs(input_options)

