"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import ChopIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from pymatgen.io.vasp import Poscar
from MAST.utility import dirutil
from MAST.utility import MASTFile
import shutil

testname="chop_test_ready"
testdir = dirutil.get_test_dir(testname)

class TestIsReadyToRunIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

        if not os.path.isdir("writedir"):
            os.mkdir("writedir")
        if not os.path.isdir("writedir/single_label1"):
            os.mkdir("writedir/single_label1")

    def tearDown(self):
        tearlist = list()
        tearlist.append("writedir")
        for foldername in tearlist:
            try:
                shutil.rmtree(foldername)
            except OSError:
                pass

    def test___init__(self):
        self.assertTrue(True)
        #raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_ready_singlerun(self):
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
        my_structure = Poscar.from_file("files/perfect_structure").structure
        mywr = ChopIngredient(name=ingdir, program_keys = kdict, structure=my_structure)
        mywr.write_singlerun()
        myrdi = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        self.assertTrue(myrdi.ready_singlerun())
        os.remove("%s/POSCAR" % ingdir)
        self.assertFalse(myrdi.ready_singlerun())
        #self.testclass.ready_singlerun()

    def test_ready_structure(self):
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
        my_pos = Poscar.from_file("files/perfect_structure")
        my_pos.write_file("writedir/single_label1/POSCAR")
        myrdi = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_pos.structure)
        self.assertTrue(myrdi.ready_structure())
        os.remove("%s/POSCAR" % ingdir)
        self.assertFalse(myrdi.ready_structure())
        #self.testclass.ready_structure()

    def test_ready_neb_subfolders(self):
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
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images'] = 3
        my_structure = Poscar.from_file("files/perfect_structure").structure
        mywr = ChopIngredient(name=ingdir, program_keys = kdict, structure=my_structure)
        for subdir in ['00','01','02','03','04']:
            subname = "%s/%s" % (ingdir, subdir)
            os.mkdir(subname)
            mywr.keywords['name'] = subname
            mywr.checker.keywords['name'] = subname 
            if not subdir in ['00','04']:
                mywr.write_singlerun()
                mywr.write_submit_script()
        myrdi = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        self.assertTrue(myrdi.ready_neb_subfolders())
        os.remove("%s/01/POSCAR" % ingdir)
        self.assertFalse(myrdi.ready_neb_subfolders())
        #self.testclass.ready_neb_subfolders()

    def test_ready_subfolders(self):
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
        kdict['images'] = 3
        my_structure = Poscar.from_file("files/perfect_structure").structure
        mywr = ChopIngredient(name=ingdir, program_keys = kdict, structure=my_structure)
        for subdir in ['00','01','02','03','04']:
            subname = "%s/%s" % (ingdir, subdir)
            os.mkdir(subname)
            mywr.keywords['name'] = subname
            mywr.checker.keywords['name'] = subname 
            mywr.write_singlerun()
            mywr.write_submit_script()
        myrdi = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        self.assertTrue(myrdi.ready_subfolders())
        os.remove("%s/00/POSCAR" % ingdir)
        self.assertFalse(myrdi.ready_subfolders())
        #self.testclass.ready_subfolders()

