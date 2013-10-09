"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import RunIngredient
from MAST.ingredients.chopingredient import WriteIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile
import shutil

testname="chop_test_run"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestRunIngredient(unittest.TestCase):

    def setUp(self):
        os.environ['MAST_CONTROL'] = testdir + "/test_control"
        os.environ['MAST_RECIPE_PATH'] = testdir
        os.environ['MAST_SCRATCH'] = testdir
        os.chdir(testdir)

        if not os.path.isdir("writedir"):
            os.mkdir("writedir")
        if not os.path.isdir("test_control"):
            os.mkdir("test_control")
        if not os.path.isdir("writedir/next_ingred"):
            os.mkdir("writedir/next_ingred")
        if not os.path.isdir("writedir/single_label1"):
            os.mkdir("writedir/single_label1")
        if not os.path.isdir("writedir/neb_labelinit-labelfin_stat"):
            os.mkdir("writedir/neb_labelinit-labelfin_stat")
        if not os.path.isdir("writedir/single_phonon_label1"):
            os.mkdir("writedir/single_phonon_label1")

    def tearDown(self):
        tearlist = list()
        #tearlist.append("writedir/single_label1")
        tearlist.append("writedir/next_ingred")
        tearlist.append("writedir/neb_labelinit-labelfin_stat")
        tearlist.append("writedir/single_phonon_label1")
        for foldername in tearlist:
            try:
                shutil.rmtree(foldername)
            except OSError:
                pass
        removelist = list()
        removelist.append("test_control/submitlist")
        for ritem in removelist:
            try:
                os.remove(ritem)
            except OSError:
                pass

    def test___init__(self):
        self.assertTrue(True)
        #raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_run_singlerun(self):
        raise SkipTest
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        kdict['mast_kpoints'] = [2,2,2,"M"]
        kdict['mast_xc'] = 'pw91'
        my_structure = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywr = WriteIngredient(name=ingdir, program_keys = kdict, structure=my_structure)
        mywr.write_singlerun()
        myri = RunIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_singlerun()
        self.assertTrue(myri.checker.is_ready_to_run())
        mysubmit = MASTFile("test_control/submitlist")
        self.assertEquals(mysubmit.data[0], "%s\n" % ingdir)
        #self.testclass.run_singlerun(mode='serial')

    def test_run_neb_subfolders(self):
        ingdir="%s/writedir/neb_labelinit-labelfin_stat" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        kdict['mast_kpoints'] = [2,2,2,"M"]
        kdict['mast_xc'] = 'pw91'
        kdict['images'] = 3
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['neblines']=dict()
        kdict['neblines']['labelinit-labelfin']=neblines
        str_00 = MASTFile("files/POSCAR_00")
        str_00.to_file("%s/parent_structure_labelinit" % ingdir) 
        str_04 = MASTFile("files/POSCAR_04")
        str_04.to_file("%s/parent_structure_labelfin" % ingdir) 
        str_01 = MASTFile("files/POSCAR_01")
        str_01.to_file("%s/parent_structure_labelinit-labelfin_01" % ingdir) 
        str_02 = MASTFile("files/POSCAR_02")
        str_02.to_file("%s/parent_structure_labelinit-labelfin_02" % ingdir) 
        str_03 = MASTFile("files/POSCAR_03")
        str_03.to_file("%s/parent_structure_labelinit-labelfin_03" % ingdir) 
        en_00 = MASTFile("files/OSZICAR_00")
        en_00.to_file("%s/parent_energy_labelinit" % ingdir) 
        en_04 = MASTFile("files/OSZICAR_04")
        en_04.to_file("%s/parent_energy_labelfin" % ingdir) 
        my_structure = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywr = WriteIngredient(name=ingdir, program_keys = kdict, structure=my_structure)
        mywr.write_neb_subfolders()
        myri = RunIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_neb_subfolders()
        mysubmit = MASTFile("test_control/submitlist")
        myri.checker.keywords['name'] = "%s/00" % ingdir
        self.assertFalse(myri.checker.is_ready_to_run()) #do not run endpoints again
        myri.checker.keywords['name'] = "%s/01" % ingdir
        self.assertTrue(myri.checker.is_ready_to_run())
        myri.checker.keywords['name'] = "%s/02" % ingdir
        self.assertTrue(myri.checker.is_ready_to_run())
        myri.checker.keywords['name'] = "%s/03" % ingdir
        self.assertTrue(myri.checker.is_ready_to_run())
        myri.checker.keywords['name'] = "%s/04" % ingdir
        self.assertFalse(myri.checker.is_ready_to_run()) #do not run endpoints again
        self.assertEquals(mysubmit.data[0], "%s/01\n" % ingdir)
        self.assertEquals(mysubmit.data[1], "%s/02\n" % ingdir)
        self.assertEquals(mysubmit.data[2], "%s/03\n" % ingdir)
        #self.testclass.run_neb_subfolders()

    def test_run_subfolders(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        kdict['mast_kpoints'] = [2,2,2,"M"]
        kdict['mast_xc'] = 'pw91'
        my_structure = pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        for subfolder in ['sub1','sub2','sub3','sub4']:
            subname = "%s/%s" % (ingdir, subfolder)
            os.mkdir(subname)
            mywr = WriteIngredient(name=subname, program_keys = kdict, structure=my_structure)
            mywr.write_singlerun()
        myri = RunIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_subfolders()
        self.assertFalse(myri.checker.is_ready_to_run())
        for subfolder in ['sub1','sub2','sub3','sub4']:
            subname = "%s/%s" % (ingdir, subfolder)
            myri.checker.keywords['name'] = subname
            self.assertTrue(myri.checker.is_ready_to_run())
        mysubmit = MASTFile("test_control/submitlist")
        self.assertEquals(mysubmit.data[0], "%s/sub1\n" % ingdir)
        self.assertEquals(mysubmit.data[1], "%s/sub2\n" % ingdir)
        self.assertEquals(mysubmit.data[2], "%s/sub3\n" % ingdir)
        self.assertEquals(mysubmit.data[3], "%s/sub4\n" % ingdir)
        #self.testclass.run_subfolders()

    def test_run_defect(self):
        raise SkipTest
        #self.testclass.run_defect()

    def test_run_strain(self):
        raise SkipTest
        #self.testclass.run_strain()

