"""Tests for Atom Indexing"""

from MAST.parsers.inputparser import InputParser
from MAST.recipe.recipesetup import RecipeSetup
from MAST.ingredients.pmgextend.atom_index import AtomIndex
import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import InputOptions
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

    def test_general_setup(self):
        #raise SkipTest #SD
        wdir=os.path.join(testdir,'workdir')
        from MAST.controllers.mastinput import MASTInput
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        print rwdir
        #iopt=InputParser(inputfile='neb_pathfinder.inp').parse()
        #rfile = iopt.get_item('personal_recipe', 'personal_recipe_list')
        #print rfile
        #struc = iopt.get_item('structure','structure')
        #myrs=RecipeSetup(recipeFile=rfile, workingDirectory=wdir, inputOptions=iopt, structure=struc)
        #myrs.start()
        #self.assertEqual(type(myrs.metafile),MAST.utility.metadata.Metadata)
        return
