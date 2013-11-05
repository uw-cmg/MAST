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
import numpy as np

testname="inputtest_notupdated"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestInputparser(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        myip = InputParser(inputfile="long_input.inp")
        s_dict = dict({\
                'mast' : myip.parse_mast_section,
                'structure' : myip.parse_structure_section,
                'ingredients' : myip.parse_ingredients_section,
                'defects' : myip.parse_defects_section,
                'recipe' : myip.parse_recipe_section,
                'neb' : myip.parse_neb_section,
                'chemical_potentials' : myip.parse_chemical_potentials_section,
                'phonon' : myip.parse_phonon_section,
                               })
        self.assertItemsEqual(myip.section_parsers, s_dict)
        #self.testclass.__init__(**kwargs)

    def test_parse(self):
        myip = InputParser(inputfile="long_input.inp")
        myoptions = myip.parse()
        #ofile = MASTFile()
        #ofile.data = str(myoptions)
        #ofile.to_file(os.path.join(os.getcwd(),"long_input_options_output.txt"))
        #self.testclass.parse()

    def test_parse_mast_section(self):
        myip = InputParser(inputfile="long_input.inp")
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
        self.assertItemsEqual(myoptions.options['mast'].keys(), ['system_name'])
        self.assertEqual(myoptions.options['mast']['system_name'],mdict['mast']['system_name']) 
        #print myoptions.options['mast']
        #self.testclass.parse_mast_section(section_name, section_content, options)

    def test_parse_structure_section_noposfile(self):
        myip = InputParser(inputfile="long_input.inp")
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
        myip = InputParser(inputfile="long_input.inp")
        minput = MASTFile("%s/structure_posfile_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_structure_section('structure',cleanlines,myoptions)
        #print myoptions
        mdict=dict()
        mdict['structure'] = dict()
        mdict['structure']['posfile'] = os.path.join(testdir, "POSCAR_startpos")
        self.assertItemsEqual(myoptions.options.keys(), ['structure'])
        self.assertEqual(myoptions.options['structure']['posfile'],mdict['structure']['posfile']) 
        myip = InputParser(inputfile="long_input.inp")
        minput = MASTFile("%s/structure_posfile_lines_nosuch.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        self.assertRaises(MASTError, myip.parse_structure_section,'structure',cleanlines,myoptions)
        #self.testclass.parse_structure_section(section_name, section_content, options)

    def test_parse_defects_section(self):
        myip = InputParser(inputfile="long_input.inp")
        minput = MASTFile("%s/defects_lines.txt" % testdir)
        cleanlines = list()
        for line in minput.data:
            cleanlines.append(line.strip())
        myoptions = InputOptions()
        myip.parse_defects_section('defects',cleanlines,myoptions)
        print myoptions.options['defects']
        mdict=dict()
        mdict['coord_type']='fractional'
        mdict['num_defects']=4
        mdict['defects']=dict()
        mdict['defects']['defect_3']=dict()
        mdict['defects']['defect_3']['threshold']=0.0001
        mdict['defects']['defect_3']['charge']=[3]
        mdict['defects']['defect_3']['coord_type']='fractional'
        mdict['defects']['defect_3']['subdefect_1']=dict()
        mdict['defects']['defect_3']['subdefect_1']['symbol']='Al'
        mdict['defects']['defect_3']['subdefect_1']['type']='vacancy'
        mdict['defects']['defect_3']['subdefect_1']['coordinates']=np.array([0.,0.,0.],'float')
        mdict['defects']['Vac@Al-Sub@Fe']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['threshold']=0.0001
        mdict['defects']['Vac@Al-Sub@Fe']['charge']=[-2,-1,0,1,2,3,4,5]
        mdict['defects']['Vac@Al-Sub@Fe']['coord_type']='fractional'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']['symbol']='Al'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']['type']='vacancy'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_1']['coordinates']=np.array([0.,0.,0.],'float')
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']=dict()
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']['symbol']='Fe'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']['type']='substitution'
        mdict['defects']['Vac@Al-Sub@Fe']['subdefect_2']['coordinates']=np.array([0.5,0.5,0.],'float')
        

        mdict['defects']['defect_3']['subdefect_1']=dict()
        mdict['defects']['defect_3']['subdefect_1']=dict()

        mdict['system_name'] = "SystemName!"
        self.assertItemsEqual(myoptions.options['defects'],mdict) 
        #self.testclass.parse_defects_section(section_name, section_content, options)

    def test_parse_recipe_section(self):
        raise SkipTest
        #self.testclass.parse_recipe_section(section_name, section_content, options)

    def test_parse_ingredients_section(self):
        raise SkipTest
        #self.testclass.parse_ingredients_section(section_name, section_content, options)

    def test_parse_neb_section(self):
        raise SkipTest
        #self.testclass.parse_neb_section(section_name, section_content, options)

    def test_parse_chemical_potentials_section(self):
        raise SkipTest
        #self.testclass.parse_chemical_potentials_section(section_name, section_content, options)

    def test_parse_phonon_section(self):
        raise SkipTest
        #self.testclass.parse_phonon_section(section_name, section_content, options)

    def test_perform_element_mapping(self):
        raise SkipTest
        #self.testclass.perform_element_mapping(input_options)

    def test_validate_execs(self):
        raise SkipTest
        #self.testclass.validate_execs(input_options)

    def test_set_structure_from_inputs(self):
        raise SkipTest
        #self.testclass.set_structure_from_inputs(input_options)

