"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import ChopIngredient
import filecmp
import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile
from MAST.utility import MASTError
import shutil
import numpy as np
testname="fsstest"
testdir = dirutil.get_test_dir(testname)
old_control = os.getenv("MAST_CONTROL")
old_recipe = os.getenv("MAST_RECIPE_PATH")
old_scratch = os.getenv("MAST_SCRATCH")

class TestFSSChopIngredient(unittest.TestCase):

    def setUp(self):
        os.environ['MAST_CONTROL'] = testdir + "/test_control"
        os.environ['MAST_RECIPE_PATH'] = testdir
        os.environ['MAST_SCRATCH'] = testdir
        os.chdir(testdir)
        if not os.path.isdir("%s/test_control" % testdir):
            os.mkdir("%s/test_control" % testdir)
        if not os.path.isdir("%s/writedir" % testdir):
            os.mkdir("%s/writedir" % testdir)
        if not os.path.isdir("%s/writedir/inducedefect_label1" % testdir):
            os.mkdir("%s/writedir/inducedefect_label1" % testdir)
        if not os.path.isdir("%s/writedir/child_mislabeled" % testdir):
            os.mkdir("%s/writedir/child_mislabeled" % testdir)
        if not os.path.isdir("%s/writedir/child_scalewhat" % testdir):
            os.mkdir("%s/writedir/child_scalewhat" % testdir)
        if not os.path.isdir("%s/writedir/child_scale5" % testdir):
            os.mkdir("%s/writedir/child_scale5" % testdir)
    def tearDown(self):
        tearlist = list()
        tearlist.append("writedir")
        for foldername in tearlist:
            try:
                shutil.rmtree(foldername)
            except OSError:
                pass
        os.environ['MAST_CONTROL'] = old_control
        os.environ['MAST_RECIPE_PATH'] = old_recipe
        os.environ['MAST_SCRATCH'] = old_scratch

    def test___init__(self):
        self.assertTrue(True)
        #raise SkipTest
        #self.testclass.__init__(**kwargs)


    def test_run_supercell_defect_set(self):
        ingdir="%s/writedir/inducedefect_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = label1\n")
        metad.to_file("%s/metadata.txt" % ingdir)
        supercellfile = MASTFile("gensc/supercell_list.txt")
        supercellfile.to_file("%s/supercell_list.txt" % ingdir)
        kdict=dict()
        kdict['label1']=dict()
        kdict['label1']['subdefect1']=dict()
        kdict['label1']['subdefect1']['symbol']='Cr'
        kdict['label1']['subdefect1']['type']='interstitial'
        kdict['label1']['subdefect1']['coordinates']=np.array([0.8, 0.7, 0.6])
        kdict['label1']['subdefect2']=dict()
        kdict['label1']['subdefect2']['symbol']='Sr'
        kdict['label1']['subdefect2']['type']='antisite'
        kdict['label1']['subdefect2']['coordinates']=np.array([0.5,0.5,0.0])
        kdict['label1']['subdefect5']=dict()
        kdict['label1']['subdefect5']['symbol']='Al'
        kdict['label1']['subdefect5']['type']='vacancy'
        kdict['label1']['subdefect5']['coordinates']=np.array([0.0,0.5,0.5])
        kdict['label1']['coord_type'] = 'fractional'
        kdict['label1']['threshold'] = 0.01
        kdict['label1']['charge'] = '2'
        ddict=dict()
        ddict['mast_defect_settings']=dict()
        ddict['mast_defect_settings'].update(kdict['label1']) #single defect grouping
        ddict['mast_program'] = 'vasp'
        myri = ChopIngredient(name=ingdir,program_keys=ddict, structure=None)
        iscomplete = myri.complete_supercell_defect_set()
        self.assertFalse(iscomplete)
        myri.run_supercell_defect_set("gensc")
        subfolder_list = dirutil.immediate_subdirs(ingdir, 1)
        for subfolder in subfolder_list:
            my_defected = pymatgen.io.vaspio.Poscar.from_file("%s/%s/CONTCAR" % (ingdir, subfolder)).structure.get_sorted_structure()
            defected_compare = pymatgen.io.vaspio.Poscar.from_file("files/induced_defects/%s/CONTCAR" % subfolder).structure.get_sorted_structure()
            self.assertEquals(my_defected, defected_compare)
        iscomplete = myri.complete_supercell_defect_set()
        self.assertTrue(iscomplete)
        #Remove one of the CONTCAR files (from last folder)
        os.remove(os.path.join(ingdir, subfolder, "CONTCAR"))
        iscomplete = myri.complete_supercell_defect_set()
        self.assertFalse(iscomplete)

    def test_give_supercell_subfolder_file(self):
        ingdir="%s/writedir/induced_defects" % testdir
        shutil.copytree("files/induced_defects", "%s/writedir/induced_defects" % testdir)
        ingdirchildm="%s/writedir/child_mislabeled" % testdir
        ingdirchildw="%s/writedir/child_scalewhat" % testdir
        ingdirchild5="%s/writedir/child_scale5" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = label1\n")
        metad.to_file("%s/metadata.txt" % ingdir)
        metad.to_file("%s/metadata.txt" % ingdirchildm)
        metad.to_file("%s/metadata.txt" % ingdirchildw)
        metad.to_file("%s/metadata.txt" % ingdirchild5)
        ddict=dict()
        ddict['mast_program'] = 'vasp'
        myri = ChopIngredient(name=ingdir,program_keys=ddict, structure=None)
        self.assertRaises(MASTError, myri.give_supercell_subfolder_file,"CONTCAR","POSCAR",ingdirchildm)
        self.assertRaises(ValueError, myri.give_supercell_subfolder_file,"CONTCAR","POSCAR",ingdirchildw) 
        myri.give_supercell_subfolder_file("CONTCAR","POSCAR",ingdirchild5) 
        myri.give_supercell_subfolder_file("KPOINTS","KPOINTS",ingdirchild5)
        self.assertTrue(filecmp.cmp("%s/POSCAR" % ingdirchild5, "%s/4x4x4/CONTCAR" % ingdir))
        self.assertFalse(filecmp.cmp("%s/POSCAR" % ingdirchild5, "%s/2x2x2/CONTCAR" % ingdir))
        self.assertTrue(filecmp.cmp("%s/KPOINTS" % ingdirchild5, "%s/4x4x4/KPOINTS" % ingdir))

