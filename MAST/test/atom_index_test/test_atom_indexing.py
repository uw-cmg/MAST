"""Tests for Atom Indexing"""

from MAST.parsers.inputparser import InputParser
from MAST.recipe.recipesetup import RecipeSetup
from MAST.ingredients.pmgextend.atom_index import AtomIndex
from MAST.utility import MASTError
import unittest
import numpy as np
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import InputOptions
from MAST.controllers.mastinput import MASTInput
from MAST.utility import MASTFile
import shutil
testname="atom_index_test"
testdir = dirutil.get_test_dir(testname)

class TestAtomIndexing(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir(os.path.join(testdir,'workdir')):
            os.mkdir('workdir')

    def tearDown(self):
        pass
        #if os.path.isdir(os.path.join(testdir,'workdir')):
        #    shutil.rmtree('workdir')
        #if os.path.isdir(os.path.join(testdir,'structure_index_files')):
        #    shutil.rmtree('structure_index_files')

    def test_interstitial_neb_setup(self):
        raise SkipTest #SD
        wdir=os.path.join(testdir,'workdir')
        tdir=os.path.join(testdir,'interstitial_neb_files')
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        #print rwdir
        fnames = os.listdir(os.path.join(rwdir, "structure_index_files"))
        for fname in fnames:
            f1=MASTFile(os.path.join(tdir,"structure_index_files",fname))
            f2=MASTFile(os.path.join(rwdir,"structure_index_files",fname))
            self.assertEqual(f1.data,f2.data)
        #iopt=InputParser(inputfile='neb_pathfinder.inp').parse()
        #rfile = iopt.get_item('personal_recipe', 'personal_recipe_list')
        #print rfile
        #struc = iopt.get_item('structure','structure')
        #myrs=RecipeSetup(recipeFile=rfile, workingDirectory=wdir, inputOptions=iopt, structure=struc)
        #myrs.start()
        #self.assertEqual(type(myrs.metafile),MAST.utility.metadata.Metadata)
        return

    def test__init__(self):
        raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=mysid)
        self.assertEqual(myai.atomcount, 1)
        return
    
    def test__init__2(self):
        raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='multidefect.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=mysid)
        self.assertEqual(myai.atomcount, 1)
        return
    
    def test_write_defected(self):
        raise SkipTest
        tdir=os.path.join(testdir,'multidefect_files')
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='multidefect.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        fnames = os.listdir(os.path.join(rwdir, "structure_index_files"))
        for fname in fnames:
            f1=MASTFile(os.path.join(tdir,"structure_index_files",fname))
            f2=MASTFile(os.path.join(rwdir,"structure_index_files",fname))
            self.assertEqual(f1.data,f2.data)
        return

    def test_find_orig_frac_coord_in_atom_indices(self):
        raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        test_sid = os.path.join(testdir, "find_coord_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        orig_coord=np.array([0.0,0.0,0.5],'float')
        self.assertRaises(MASTError, myai.find_orig_frac_coord_in_atom_indices,
                orig_coord,"")
        print "subtest1 ok"
        findtest2 = myai.find_orig_frac_coord_in_atom_indices(orig_coord,"",
                scaling_label="",find_multiple=True)
        self.assertListEqual(findtest2, ["0000000000000x12","0000000000000xE2"])
        print "subtest2 ok"
        findtest3 = myai.find_orig_frac_coord_in_atom_indices(orig_coord,"He")
        self.assertEqual(findtest3, "0000000000000xE2")
        print "subtest3 ok"
        findtest4 = myai.find_orig_frac_coord_in_atom_indices(orig_coord,"Al",
            find_multiple=True, tol=0.001)
        self.assertListEqual(findtest4, ["0000000000000x12","0000000000000xTOL"])
        print "subtest4 ok"
        return
    
    def test_find_any_frac_coord_in_atom_indices(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        test_sid = os.path.join(testdir, "find_coord_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        relaxed_coord=np.array([-0.00501174,-0.00501174, 0.50707951],'float')
        findtest1 = myai.find_any_frac_coord_in_atom_indices(relaxed_coord,"",
                scaling_label="",find_multiple=False)
        print findtest1
        self.assertRaises(MASTError, myai.find_any_frac_coord_in_atom_indices,
                relaxed_coord,"")
        print "subtest1 ok"
        findtest2 = myai.find_any_frac_coord_in_atom_indices(relaxed_coord,"",
                scaling_label="",find_multiple=True)
        self.assertListEqual(findtest2, ["0000000000000x12","0000000000000xE2"])
        print "subtest2 ok"
        findtest3 = myai.find_any_frac_coord_in_atom_indices(relaxed_coord,"He")
        self.assertEqual(findtest3, "0000000000000xE2")
        print "subtest3 ok"
        findtest4 = myai.find_any_frac_coord_in_atom_indices(relaxed_coord,"Al",
            find_multiple=True, tol=0.001)
        self.assertListEqual(findtest4, ["0000000000000x12","0000000000000xTOL"])
        print "subtest4 ok"
        return
