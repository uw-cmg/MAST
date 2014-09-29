"""Tests for Recipeplan"""

from MAST.recipe.recipeplan import RecipePlan

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
import shutil
from MAST.utility import dirutil
from MAST.utility import MASTFile
from MAST.utility import MASTError
testname="recipeplan_test"
testdir = dirutil.get_test_dir(testname)
old_control = os.getenv("MAST_CONTROL")
old_scratch = os.getenv("MAST_SCRATCH")
old_recipe = os.getenv("MAST_RECIPE_PATH")

class TestRecipeplan(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir("test_control"):
            os.mkdir("test_control")
        if not os.path.isdir("recipedir"):
            os.mkdir("recipedir")
        if not os.path.isdir("recipedir/ing2a"):
            os.mkdir("recipedir/ing2a")
        if not os.path.isdir("recipedir/ing2b"):
            os.mkdir("recipedir/ing2b")
        if not os.path.isdir("recipedir/ing3"):
            os.mkdir("recipedir/ing3")
        os.environ['MAST_CONTROL'] = testdir + "/test_control"
        os.environ['MAST_RECIPE_PATH'] = testdir
        os.environ['MAST_SCRATCH'] = testdir

    def tearDown(self):
        removelist=list()
        removelist.append("recipedir/status.txt")
        removelist.append("recipedir/ing2b/INCAR")
        removelist.append("recipedir/ing2b/POSCAR")
        removelist.append("recipedir/ing2b/POTCAR")
        removelist.append("recipedir/ing2b/KPOINTS")
        removelist.append("recipedir/ing2b/submit.sh")
        removelist.append("recipedir/ing2b/CHGCAR")
        removelist.append("recipedir/ing2b/WAVECAR")
        removelist.append("recipedir/ing2a/POSCAR")
        removelist.append("recipedir/printed.txt")
        removelist.append("test_control/submitlist")
        removelist.append("test_control/set_platform")
        for myfile in removelist:
            try:
                os.remove(myfile)
            except OSError:
                pass
        os.environ['MAST_CONTROL'] = old_control
        os.environ['MAST_RECIPE_PATH'] = old_recipe
        os.environ['MAST_SCRATCH'] = old_scratch

    def test___init__(self):
        rp = RecipePlan("recipedir")
        self.assertEquals(rp.working_directory,"recipedir")
        self.assertEquals(rp.status,"I")
        #self.assertEquals(rp.name,"test_recipe")
        #raise SkipTest
        #self.testclass.__init__(name, working_directory)

    def test_write_ingredient(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        #metad = MASTFile("files/metadata_single")
        #metad.to_file("%s/metadata.txt" % ingdir)
        rp = RecipePlan("recipedir")
        rp.ingredients['ing2b'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        rp.ingred_input_options['ing2b']=dict()
        rp.ingred_input_options['ing2b']['name']="recipedir/ing2b"
        rp.ingred_input_options['ing2b']['program_keys']=kdict
        rp.ingred_input_options['ing2b']['structure']=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.write_methods['ing2b']=[['write_singlerun']]
        rp.write_ingredient('ing2b')
        self.assertTrue(os.path.isfile('recipedir/ing2b/INCAR'))
        self.assertTrue(os.path.isfile('recipedir/ing2b/POSCAR'))
        self.assertTrue(os.path.isfile('recipedir/ing2b/POTCAR'))
        self.assertTrue(os.path.isfile('recipedir/ing2b/KPOINTS'))
        self.assertTrue(os.path.isfile('recipedir/ing2b/submit.sh'))
        #self.testclass.write_ingredient(iname)

    def test_complete_ingredient(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        #metad = MASTFile("files/metadata_single")
        #metad.to_file("%s/metadata.txt" % ingdir)
        rp = RecipePlan("recipedir")
        rp.ingredients['ing1'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        rp.ingred_input_options['ing1']=dict()
        rp.ingred_input_options['ing1']['name']="%s/recipedir/ing1" % testdir
        rp.ingred_input_options['ing1']['program_keys']=kdict
        rp.ingred_input_options['ing1']['structure']=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.complete_methods['ing1']=[['complete_singlerun']]
        self.assertTrue(rp.complete_ingredient('ing1'))
        #self.testclass.complete_ingredient(iname)

    def test_ready_ingredient(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        #metad = MASTFile("files/metadata_single")
        #metad.to_file("%s/metadata.txt" % ingdir)
        rp = RecipePlan("recipedir")
        rp.ingredients['ing2b'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        rp.ingred_input_options['ing2b']=dict()
        rp.ingred_input_options['ing2b']['name']="recipedir/ing2b"
        rp.ingred_input_options['ing2b']['program_keys']=kdict
        rp.ingred_input_options['ing2b']['structure']=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.write_methods['ing2b']=[['write_singlerun']]
        rp.write_ingredient('ing2b')
        rp.ready_methods['ing2b']=[['ready_singlerun']]
        self.assertTrue(rp.ready_ingredient('ing2b'))
        #self.testclass.ready_ingredient(iname)

    def test_run_ingredient(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        #metad = MASTFile("files/metadata_single")
        #metad.to_file("%s/metadata.txt" % ingdir)
        rp = RecipePlan("recipedir")
        rp.ingredients['ing2b'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        rp.ingred_input_options['ing2b']=dict()
        rp.ingred_input_options['ing2b']['name']="recipedir/ing2b"
        rp.ingred_input_options['ing2b']['program_keys']=kdict
        rp.ingred_input_options['ing2b']['structure']=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.write_methods['ing2b']=[['write_singlerun']]
        rp.write_ingredient('ing2b')
        rp.run_methods['ing2b']=[['run_singlerun']]
        rp.run_ingredient('ing2b')
        mysubmit = MASTFile("test_control/submitlist")
        self.assertEquals(mysubmit.data[0],"recipedir/ing2b\n")
        #self.testclass.run_ingredient(iname)

    def test_update_children(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        #metad = MASTFile("files/metadata_single")
        #metad.to_file("%s/metadata.txt" % ingdir)
        rp = RecipePlan("%s/recipedir" % testdir)
        rp.ingredients['ing1'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        rp.ingred_input_options['ing1']=dict()
        rp.ingred_input_options['ing1']['name']="%s/ing1" % rp.working_directory
        rp.ingred_input_options['ing1']['program_keys']=kdict
        rp.ingred_input_options['ing1']['structure']=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.update_methods['ing1']=dict()
        rp.update_methods['ing1']['ing2a']=[['give_structure']]
        rp.update_methods['ing1']['ing2b']=[['give_structure_and_restart_files']]
        rp.update_children('ing1')
        self.assertTrue(os.path.isfile("recipedir/ing2a/POSCAR"))
        self.assertTrue(os.path.isfile("recipedir/ing2b/POSCAR"))
        #CHGCAR softlink only sent to second child
        self.assertFalse(os.path.exists("recipedir/ing2a/CHGCAR"))
        self.assertTrue(os.path.exists("recipedir/ing2b/CHGCAR"))
        #self.testclass.update_children(iname)

    def test_fast_forward_check_complete(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing1/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = labela\n")
        metad.to_file("recipedir/ing2a/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = labelb\n")
        metad.to_file("recipedir/ing2b/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing3/metadata.txt")
        rp = RecipePlan("recipedir")
        rp.ingredients['ing1'] = "I"
        rp.ingredients['ing2a'] = "I"
        rp.ingredients['ing2b'] = "I"
        rp.ingredients['ing3'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        my_struc = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.ingred_input_options['ing1']=dict()
        rp.ingred_input_options['ing1']['name']="%s/recipedir/ing1" % testdir
        rp.ingred_input_options['ing1']['program_keys']=kdict
        rp.ingred_input_options['ing1']['structure']=my_struc
        rp.complete_methods['ing1']=[['complete_singlerun']]
        rp.update_methods['ing1']=dict()
        rp.update_methods['ing1']['ing2a']=[['give_structure']]
        rp.update_methods['ing1']['ing2b']=[['give_structure']]
        rp.ingred_input_options['ing2a']=dict()
        rp.ingred_input_options['ing2a']['name']="%s/recipedir/ing2a" % testdir
        rp.ingred_input_options['ing2a']['program_keys']=kdict
        rp.ingred_input_options['ing2a']['structure']=my_struc
        rp.complete_methods['ing2a']=[['complete_singlerun']]
        rp.ready_methods['ing2a']=[['ready_structure']]
        rp.ingred_input_options['ing2b']=dict()
        rp.ingred_input_options['ing2b']['name']="%s/recipedir/ing2b" % testdir
        rp.ingred_input_options['ing2b']['program_keys']=kdict
        rp.ingred_input_options['ing2b']['structure']=my_struc
        rp.complete_methods['ing2b']=[['complete_singlerun']]
        rp.ready_methods['ing2b']=[['ready_structure']]
        rp.ingred_input_options['ing3']=dict()
        rp.ingred_input_options['ing3']['name']="%s/recipedir/ing3" % testdir
        rp.ingred_input_options['ing3']['program_keys']=kdict
        rp.ingred_input_options['ing3']['structure']=my_struc
        rp.complete_methods['ing3']=[['complete_singlerun']]
        rp.ready_methods['ing3']=[['ready_structure']]
        rp.fast_forward_check_complete()
        #self.assertTrue(rp.complete_ingredient('ing1'))
        self.assertEquals(rp.ingredients, {'ing1':'C','ing2a':'I','ing2b':'I','ing3':'I'})
        self.assertTrue(rp.ready_ingredient('ing2a'))
        self.assertTrue(rp.ready_ingredient('ing2b'))
        #self.testclass.fast_forward_check_complete()

    def test_check_if_have_parents(self):
        rp = RecipePlan("recipedir")
        rp.ingredients['ing1'] = "I"
        rp.ingredients['ing2a'] = "I"
        rp.ingredients['ing2b'] = "I"
        rp.ingredients['ing3'] = "I"
        rp.parents_to_check['ing3']=['ing2a','ing2b']
        rp.parents_to_check['ing2a']=['ing1']
        rp.parents_to_check['ing2b']=[]
        rp.parents_to_check['ing1']=[]
        rp.check_if_have_parents()
        self.assertEquals(rp.ingredients,{'ing1':'S','ing2a':'W','ing2b':'S','ing3':'W'})
        #self.testclass.check_if_have_parents()

    def test_check_if_ready_to_proceed_are_complete(self):
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing1/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing2a/metadata.txt")
        rp = RecipePlan("recipedir")
        rp.ingredients['ing1'] = "P"
        rp.ingredients['ing2a'] = "I"
        rp.ingredients['ing2b'] = "I"
        rp.ingredients['ing3'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        my_struc = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.ingred_input_options['ing1']=dict()
        rp.ingred_input_options['ing1']['name']="%s/recipedir/ing1" % testdir
        rp.ingred_input_options['ing1']['program_keys']=kdict
        rp.ingred_input_options['ing1']['structure']=my_struc
        rp.complete_methods['ing1']=[['complete_singlerun']]
        rp.update_methods['ing1']=dict()
        rp.update_methods['ing1']['ing2a']=[['give_structure']]
        rp.update_methods['ing1']['ing2b']=[['give_structure']]
        rp.ingred_input_options['ing2a']=dict()
        rp.ingred_input_options['ing2a']['name']="%s/recipedir/ing2a" % testdir
        rp.ingred_input_options['ing2a']['program_keys']=kdict
        rp.ingred_input_options['ing2a']['structure']=my_struc
        rp.complete_methods['ing2a']=[['complete_singlerun']]
        rp.ready_methods['ing2a']=[['ready_structure']]
        rp.check_if_ready_to_proceed_are_complete()
        self.assertTrue(rp.ready_ingredient('ing2a'))
        self.assertEquals
        self.assertEquals(rp.ingredients,{'ing1':'C','ing2a':'I','ing2b':'I','ing3':'I'})
        #self.testclass.check_if_ready_to_proceed_are_complete()

    def test_check_if_parents_are_complete(self):
        rp = RecipePlan("recipedir")
        rp.ingredients['ing1'] = "C"
        rp.ingredients['ing2a'] = "W"
        rp.ingredients['ing2b'] = "I"
        rp.ingredients['ing3'] = "I"
        rp.parents_to_check['ing3']=['ing2a','ing2b']
        rp.parents_to_check['ing2a']=['ing1']
        rp.parents_to_check['ing2b']=[]
        rp.parents_to_check['ing1']=[]
        rp.check_if_parents_are_complete()
        self.assertEquals(rp.ingredients,{'ing1':'C','ing2a':'S','ing2b':'I','ing3':'I'})
        #self.testclass.check_if_parents_are_complete()

    def test_run_staged_ingredients(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        #metad = MASTFile("files/metadata_single")
        #metad.to_file("%s/metadata.txt" % ingdir)
        rp = RecipePlan("recipedir")
        rp.ingredients['ing1']="C"
        rp.ingredients['ing2a'] = "W"
        rp.ingredients['ing2b'] = "S"
        rp.ingredients['ing3'] = "W"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        rp.ingred_input_options['ing2b']=dict()
        rp.ingred_input_options['ing2b']['name']="recipedir/ing2b"
        rp.ingred_input_options['ing2b']['program_keys']=kdict
        rp.ingred_input_options['ing2b']['structure']=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.write_methods['ing2b']=[['write_singlerun']]
        rp.write_ingredient('ing2b')
        rp.ready_methods['ing2b']=[['ready_singlerun']]
        rp.run_methods['ing2b']=[['run_singlerun']]
        rp.complete_methods['ing2b']=[['complete_singlerun']]
        rp.run_staged_ingredients()
        mysubmit = MASTFile("test_control/submitlist")
        self.assertEquals(mysubmit.data[0],"recipedir/ing2b\n")
        self.assertEquals(rp.ingredients,{'ing1':'C','ing2a':'W','ing2b':'P','ing3':'W'})
        #self.testclass.run_staged_ingredients()

    def test_check_recipe_status(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing1/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = labela\n")
        metad.to_file("recipedir/ing2a/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = labelb\n")
        metad.to_file("recipedir/ing2b/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing3/metadata.txt")
        rp = RecipePlan("%s/recipedir" % testdir)
        rp.ingredients['ing1'] = "I"
        rp.ingredients['ing2a'] = "I"
        rp.ingredients['ing2b'] = "I"
        rp.ingredients['ing3'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        my_struc = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.ingred_input_options['ing1']=dict()
        rp.ingred_input_options['ing1']['name']="%s/recipedir/ing1" % testdir
        rp.ingred_input_options['ing1']['program_keys']=kdict
        rp.ingred_input_options['ing1']['structure']=my_struc
        rp.complete_methods['ing1']=[['complete_singlerun']]
        rp.run_methods['ing1']=[['run_singlerun']]
        rp.update_methods['ing1']=dict()
        rp.update_methods['ing1']['ing2a']=[['give_structure']]
        rp.update_methods['ing1']['ing2b']=[['give_structure']]
        rp.ingred_input_options['ing2a']=dict()
        rp.ingred_input_options['ing2a']['name']="%s/recipedir/ing2a" % testdir
        rp.ingred_input_options['ing2a']['program_keys']=kdict
        rp.ingred_input_options['ing2a']['structure']=my_struc
        rp.complete_methods['ing2a']=[['complete_singlerun']]
        rp.ready_methods['ing2a']=[['ready_structure']]
        rp.run_methods['ing2a']=[['run_singlerun']]
        rp.ingred_input_options['ing2b']=dict()
        rp.ingred_input_options['ing2b']['name']="%s/recipedir/ing2b" % testdir
        rp.ingred_input_options['ing2b']['program_keys']=kdict
        rp.ingred_input_options['ing2b']['structure']=my_struc
        rp.complete_methods['ing2b']=[['complete_singlerun']]
        rp.ready_methods['ing2b']=[['ready_structure']]
        rp.run_methods['ing2b']=[['run_singlerun']]
        rp.ingred_input_options['ing3']=dict()
        rp.ingred_input_options['ing3']['name']="%s/recipedir/ing3" % testdir
        rp.ingred_input_options['ing3']['program_keys']=kdict
        rp.ingred_input_options['ing3']['structure']=my_struc
        rp.complete_methods['ing3']=[['complete_singlerun']]
        rp.ready_methods['ing3']=[['ready_structure']]
        rp.run_methods['ing3']=[['run_singlerun']]
        rp.parents_to_check['ing3']=['ing2a','ing2b']
        rp.parents_to_check['ing2a']=['ing1']
        rp.parents_to_check['ing2b']=[]
        rp.parents_to_check['ing1']=[]
        rp.check_recipe_status()
        mystatus = MASTFile("%s/recipedir/status.txt" % testdir)
        status_compare = MASTFile("%s/files/status_current.txt" % testdir)
        self.assertEqual(mystatus.data, status_compare.data)
        
        #self.testclass.check_recipe_status(verbose=1)

    def test_print_status(self):
        rp = RecipePlan("%s/recipedir" % testdir)
        rp.ingredients['ing1'] = "S"
        rp.ingredients['ing2a'] = "W"
        rp.ingredients['ing2b'] = "C"
        rp.ingredients['ing3'] = "P"
        rp.print_status()
        mystatus = MASTFile("recipedir/status.txt")
        status_compare = MASTFile("files/status_artificial.txt")
        self.assertEqual(mystatus.data, status_compare.data)
        #self.testclass.print_status(verbose=1)

    def test___repr__(self):
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("recipedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing1/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = labela\n")
        metad.to_file("recipedir/ing2a/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = labelb\n")
        metad.to_file("recipedir/ing2b/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("recipedir/ing3/metadata.txt")
        rp = RecipePlan("recipedir")
        rp.ingredients['ing1'] = "I"
        rp.ingredients['ing2a'] = "I"
        rp.ingredients['ing2b'] = "I"
        rp.ingredients['ing3'] = "I"
        kdict=dict()
        kdict['mast_program']='vasp'
        kdict['mast_xc']='pw91'
        kdict['mast_kpoints']=[1,2,3,"G"]
        my_struc = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        rp.ingred_input_options['ing1']=dict()
        rp.ingred_input_options['ing1']['name']="recipedir/ing1"
        rp.ingred_input_options['ing1']['program_keys']=kdict
        rp.ingred_input_options['ing1']['structure']=my_struc
        rp.complete_methods['ing1']=[['complete_singlerun']]
        rp.run_methods['ing1']=[['run_singlerun']]
        rp.ready_methods['ing1']=[['ready_singlerun']]
        rp.write_methods['ing1']=[['write_singlerun']]
        rp.update_methods['ing1']=dict()
        rp.update_methods['ing1']['ing2a']=[['give_structure']]
        rp.update_methods['ing1']['ing2b']=[['give_structure']]
        rp.ingred_input_options['ing2a']=dict()
        rp.ingred_input_options['ing2a']['name']="recipedir/ing2a"
        rp.ingred_input_options['ing2a']['program_keys']=kdict
        rp.ingred_input_options['ing2a']['structure']=my_struc
        rp.complete_methods['ing2a']=[['complete_singlerun']]
        rp.ready_methods['ing2a']=[['ready_structure']]
        rp.run_methods['ing2a']=[['run_singlerun']]
        rp.write_methods['ing2a']=[['write_singlerun']]
        rp.update_methods['ing2a']=dict()
        rp.update_methods['ing2a']['ing3']=[['give_structure_and_restart_files']]
        rp.ingred_input_options['ing2b']=dict()
        rp.ingred_input_options['ing2b']['name']="recipedir/ing2b"
        rp.ingred_input_options['ing2b']['program_keys']=kdict
        rp.ingred_input_options['ing2b']['structure']=my_struc
        rp.complete_methods['ing2b']=[['complete_singlerun']]
        rp.ready_methods['ing2b']=[['ready_structure']]
        rp.run_methods['ing2b']=[['run_singlerun']]
        rp.write_methods['ing2b']=[['write_singlerun']]
        rp.update_methods['ing2b']=dict()
        rp.update_methods['ing2b']['ing3']=[['give_structure_and_restart_files']]
        rp.ingred_input_options['ing3']=dict()
        rp.ingred_input_options['ing3']['name']="recipedir/ing3"
        rp.ingred_input_options['ing3']['program_keys']=kdict
        rp.ingred_input_options['ing3']['structure']=my_struc
        rp.complete_methods['ing3']=[['complete_singlerun']]
        rp.ready_methods['ing3']=[['ready_structure']]
        rp.run_methods['ing3']=[['run_singlerun']]
        rp.write_methods['ing3']=[['write_singlerun']]
        rp.update_methods['ing3']=dict()
        rp.parents_to_check['ing3']=['ing2a','ing2b']
        rp.parents_to_check['ing2a']=['ing1']
        rp.parents_to_check['ing2b']=[]
        rp.parents_to_check['ing1']=[]
        mylines=rp.__repr__()
        printed =MASTFile()
        printed.data = mylines
        printed.to_file("recipedir/printed.txt")
        compare_print=MASTFile("files/printout.txt")
        printedread = MASTFile("recipedir/printed.txt")
        self.assertEqual(printedread.data, compare_print.data)

        #self.testclass.__repr__()

    def test_get_statuses_from_file(self):
        rp = RecipePlan("recipedir")
        mystatus = MASTFile("files/status_random.txt")
        self.assertRaises(MASTError, rp.get_statuses_from_file)
        mystatus.to_file("recipedir/status.txt")
        rp.ingredients['ing1']="I"
        rp.ingredients['ing2a']="I"
        rp.ingredients['ing2b']="I"
        self.assertRaises(MASTError,rp.get_statuses_from_file)
        rp.ingredients['ing3']="I"
        rp.get_statuses_from_file()
        statusdict=dict()
        statusdict={'ing1':'alpha','ing2a':'beta','ing2b':'gamma','ing3':'delta'}
        self.assertEquals(rp.ingredients, statusdict)
        #self.testclass.get_statuses_from_file()

