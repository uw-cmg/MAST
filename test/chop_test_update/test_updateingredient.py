"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import UpdateChildrenIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile
import shutil

testname="chop_test_update"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestUpdateChildrenIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir("writedir"):
            os.mkdir("writedir")
        if not os.path.isdir("writedir/next_ingred"):
            os.mkdir("writedir/next_ingred")
        if not os.path.isdir("writedir/single_label1"):
            os.mkdir("writedir/single_label1")
        if not os.path.isdir("writedir/neb_labelinit-labelfin"):
            os.mkdir("writedir/neb_labelinit-labelfin")

    def tearDown(self):
        for foldername in ["writedir/single_label1","writedir/next_ingred","writedir/neb_labelinit-labelfin"]:
            try:
                shutil.rmtree(foldername)
            except OSError:
                pass

    def test___init__(self):
        self.assertTrue(True)
        #self.testclass.__init__(**kwargs)

    def test__fullpath_childname(self):
        ingdir="%s/writedir/single_label1" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        myuci = UpdateChildrenIngredient(name=ingdir,program_keys=kdict)
        fullpath = myuci._fullpath_childname("next_ingred")
        self.assertEqual(fullpath, "%s/writedir/next_ingred" % testdir)
        #self.testclass._fullpath_childname(childname)

    def test_give_structure(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        myuci = UpdateChildrenIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        #self.testclass.give_structure(childname)

    def test_give_neb_structures_to_neb(self):
        ingdir="%s/writedir/neb_labelinit-labelfin" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp_neb'
        kdict['images'] = 3
        my_structure = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        myrelaxed=dict()
        for subdir in ['00','01','02','03','04']:
            os.mkdir("writedir/neb_labelinit-labelfin/%s" % subdir)
            myrelaxed[subdir] = MASTFile("files/POSCAR_%s" % subdir)
            myrelaxed[subdir].to_file("writedir/neb_labelinit-labelfin/%s/CONTCAR" % subdir)
        myuci = UpdateChildrenIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_neb_structures_to_neb("next_ingred")
        for subdir in ['01','02','03']:
            givenstr = MASTFile("%s/writedir/next_ingred/parent_structure_labelinit-labelfin_%s" % (testdir,subdir))
            self.assertEqual(myrelaxed[subdir].data, givenstr.data)
        #self.testclass.give_neb_structures_to_neb(childname)

    def test_give_saddle_structure(self):
        raise SkipTest
        #self.testclass.give_saddle_structure(childname)

    def test_give_phonon_multiple_forces_and_displacements(self):
        raise SkipTest
        #self.testclass.give_phonon_multiple_forces_and_displacements(childname)

    def test_give_phonon_single_forces_and_displacements(self):
        raise SkipTest
        #self.testclass.give_phonon_single_forces_and_displacements(childname)

    def test_give_structure_and_energy_to_neb(self):
        raise SkipTest
        #self.testclass.give_structure_and_energy_to_neb(childname)

    def test_give_structure_and_restart_files(self):
        raise SkipTest
        #self.testclass.give_structure_and_restart_files(childname)

