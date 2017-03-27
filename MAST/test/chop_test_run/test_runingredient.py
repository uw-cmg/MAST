"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import ChopIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile
import shutil
import numpy as np
from pymatgen.io.vasp import Poscar
testname="chop_test_run"
testdir = dirutil.get_test_dir(testname)
old_control = os.getenv("MAST_CONTROL")
old_scratch = os.getenv("MAST_SCRATCH")

class TestRunIngredient(unittest.TestCase):

    def setUp(self):
        self.test_control = testdir + "/control_%s" % self._testMethodName
        os.environ["MAST_CONTROL"] = self.test_control
        os.environ['MAST_SCRATCH'] = testdir
        os.chdir(testdir)

        if not os.path.isdir("writedir"):
            os.mkdir("writedir")
        if not os.path.isdir(self.test_control):
            os.mkdir(self.test_control)
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
        tearlist.append("writedir")
        tearlist.append(self.test_control)
        for tearfolder in tearlist:
            subdirs = os.listdir(tearfolder)
            for subdir in subdirs:
                subpath = os.path.join(tearfolder, subdir)
                if os.path.isdir(subpath):
                    shutil.rmtree(subpath)
                else:
                    os.remove(subpath)
        os.environ['MAST_CONTROL'] = old_control
        os.environ['MAST_SCRATCH'] = old_scratch

    def test___init__(self):
        #raise SkipTest
        self.assertTrue(True)
        #self.testclass.__init__(**kwargs)

    def test_run_singlerun(self):
        #raise SkipTest
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
        myri = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_singlerun()
        self.assertTrue(myri.checker.is_ready_to_run())
        mysubmit = MASTFile("%s/submitlist" % self.test_control)
        self.assertEquals(mysubmit.data[0], "%s\n" % ingdir)
        #self.testclass.run_singlerun(mode='serial')

    def test_run_neb_subfolders(self):
        #raise SkipTest
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
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['lines']=neblines
        kdict['mast_neb_settings']['images']=3
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
        my_structure = Poscar.from_file("files/perfect_structure").structure
        mywr = ChopIngredient(name=ingdir, program_keys = kdict, structure=my_structure)
        mywr.write_neb_subfolders()
        myri = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_neb_subfolders()
        mysubmit = MASTFile("%s/submitlist" % self.test_control)
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
        #raise SkipTest
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
        for subfolder in ['sub1','sub2','sub3','sub4']:
            subname = "%s/%s" % (ingdir, subfolder)
            os.mkdir(subname)
            shutil.copy("files/metadata_single","%s/metadata.txt" % subname)
            mywr = ChopIngredient(name=subname, program_keys = kdict, structure=my_structure)
            mywr.write_singlerun()
        myri = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_subfolders()
        self.assertFalse(myri.checker.is_ready_to_run())
        for subfolder in ['sub1','sub2','sub3','sub4']:
            subname = "%s/%s" % (ingdir, subfolder)
            myri.checker.keywords['name'] = subname
            self.assertTrue(myri.checker.is_ready_to_run())
        mysubmit = MASTFile("%s/submitlist" % self.test_control)
        self.assertEquals(mysubmit.data[0], "%s/sub1\n" % ingdir)
        self.assertEquals(mysubmit.data[1], "%s/sub2\n" % ingdir)
        self.assertEquals(mysubmit.data[2], "%s/sub3\n" % ingdir)
        self.assertEquals(mysubmit.data[3], "%s/sub4\n" % ingdir)
        #self.testclass.run_subfolders()

    def test_run_defect(self):
        #raise SkipTest
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = label1\n")
        metad.data.append("scaling_size = [1,1,1]\n")
        metad.to_file("%s/metadata.txt" % ingdir)
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
        kdict['label1']['subdefect3']=dict()
        kdict['label1']['subdefect3']['symbol']='Cr'
        kdict['label1']['subdefect3']['type']='interstitial'
        kdict['label1']['subdefect3']['coordinates']=np.array([0.3,0.3,0.2])
        kdict['label1']['subdefect4']=dict()
        kdict['label1']['subdefect4']['symbol']='Fe'
        kdict['label1']['subdefect4']['type']='substitution'
        kdict['label1']['subdefect4']['coordinates']=np.array([0.25,0.75,0.75])
        kdict['label1']['subdefect5']=dict()
        kdict['label1']['subdefect5']['symbol']='O'
        kdict['label1']['subdefect5']['type']='vacancy'
        kdict['label1']['subdefect5']['coordinates']=np.array([0.5,0.25,0.75])
        kdict['label1']['subdefect6']=dict()
        kdict['label1']['subdefect6']['symbol']='Fe'
        kdict['label1']['subdefect6']['type']='substitution'
        kdict['label1']['subdefect6']['coordinates']=np.array([0.25,0.25,0.25])
        kdict['label1']['subdefect7']=dict()
        kdict['label1']['subdefect7']['symbol']='La'
        kdict['label1']['subdefect7']['type']='vacancy'
        kdict['label1']['subdefect7']['coordinates']=np.array([0,0,0.5])
        kdict['label1']['subdefect8']=dict()
        kdict['label1']['subdefect8']['symbol']='Ni'
        kdict['label1']['subdefect8']['type']='interstitial'
        kdict['label1']['subdefect8']['coordinates']=np.array([0.4,0.1,0.3])
        kdict['label1']['coord_type'] = 'fractional'
        kdict['label1']['threshold'] = 0.01
        kdict['label1']['charge'] = '2'
        ddict=dict()
        ddict['mast_defect_settings']=dict()
        ddict['mast_defect_settings'].update(kdict['label1']) #single defect grouping
        ddict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/POSCAR_perfect").structure
        myperf = MASTFile("files/POSCAR_perfect")
        myperf.to_file("%s/POSCAR" % ingdir)
        myri = ChopIngredient(name=ingdir,program_keys=ddict, structure=my_structure)
        myri.run_defect()
        #
        #defect = kdict['label1']
        #base_structure = my_structure.copy()
        #from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
        #for key in defect:
        #    if 'subdefect' in key:
        #        subdefect = defect[key]
        #        sxtend = StructureExtensions(struc_work1=base_structure)
        #        base_structure = sxtend.induce_defect(subdefect, defect['coord_type'], defect['threshold'])
        #        print base_structure
        #    else:
        #        pass
        #return
        #
        #
        my_defected = Poscar.from_file("%s/CONTCAR" % ingdir).structure.get_sorted_structure()
        defected_compare = Poscar.from_file("files/POSCAR_multi").structure.get_sorted_structure()
        self.assertEquals(my_defected, defected_compare)
        self.assertFalse(os.path.isfile("%s/submitlist" % self.test_control))
        #self.testclass.run_defect()

    def test_run_strain(self):
        #raise SkipTest
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        kdict['mast_strain']=" 0.98 0.92 1.03 \n"
        my_structure = Poscar.from_file("files/POSCAR_perfect").structure
        myunstrained = MASTFile("files/POSCAR_unstrained")
        myunstrained.to_file("%s/POSCAR" % ingdir)
        myri = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_strain()
        my_strained = Poscar.from_file("%s/CONTCAR" % ingdir).structure
        strained_compare = Poscar.from_file("files/POSCAR_strained").structure
        self.assertEquals(my_strained, strained_compare)
        self.assertFalse(os.path.isfile("%s/submitlist" % self.test_control))
        #self.testclass.run_strain()

    def test_run_scale_defect(self):
        #raise SkipTest
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = label1\n")
        metad.data.append("scaling_size = [2,2,2]\n")
        metad.to_file("%s/metadata.txt" % ingdir)
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
        kdict['label1']['subdefect3']=dict()
        kdict['label1']['subdefect3']['symbol']='Cr'
        kdict['label1']['subdefect3']['type']='interstitial'
        kdict['label1']['subdefect3']['coordinates']=np.array([0.3,0.3,0.2])
        kdict['label1']['subdefect4']=dict()
        kdict['label1']['subdefect4']['symbol']='Fe'
        kdict['label1']['subdefect4']['type']='substitution'
        kdict['label1']['subdefect4']['coordinates']=np.array([0.25,0.75,0.75])
        kdict['label1']['subdefect5']=dict()
        kdict['label1']['subdefect5']['symbol']='O'
        kdict['label1']['subdefect5']['type']='vacancy'
        kdict['label1']['subdefect5']['coordinates']=np.array([0.5,0.25,0.75])
        kdict['label1']['subdefect6']=dict()
        kdict['label1']['subdefect6']['symbol']='Fe'
        kdict['label1']['subdefect6']['type']='substitution'
        kdict['label1']['subdefect6']['coordinates']=np.array([0.25,0.25,0.25])
        kdict['label1']['subdefect7']=dict()
        kdict['label1']['subdefect7']['symbol']='La'
        kdict['label1']['subdefect7']['type']='vacancy'
        kdict['label1']['subdefect7']['coordinates']=np.array([0,0,0.5])
        kdict['label1']['subdefect8']=dict()
        kdict['label1']['subdefect8']['symbol']='Ni'
        kdict['label1']['subdefect8']['type']='interstitial'
        kdict['label1']['subdefect8']['coordinates']=np.array([0.4,0.1,0.3])
        kdict['label1']['coord_type'] = 'fractional'
        kdict['label1']['threshold'] = 0.01
        kdict['label1']['charge'] = '2'
        mdict=dict()
        mdict['mast_program'] = 'vasp'
        mdict['mast_defect_settings']=dict()
        mdict['mast_defect_settings'].update(kdict['label1'])
        my_structure = Poscar.from_file("files/POSCAR_perfect").structure
        myperf = MASTFile("files/POSCAR_perfect")
        myperf.to_file("%s/POSCAR" % ingdir)
        myri = ChopIngredient(name=ingdir,program_keys=mdict, structure=my_structure)
        myri.run_scale()
        #copy contcar to poscar
        os.rename("%s/CONTCAR" % ingdir, "%s/POSCAR" % ingdir)
        myri.run_defect()
        my_defected = Poscar.from_file("%s/CONTCAR" % ingdir).structure.get_sorted_structure()
        defected_compare = Poscar.from_file("files/POSCAR_scaled_defected").structure.get_sorted_structure()
        self.assertEquals(my_defected, defected_compare)
    def test_run_scale(self):
        #raise SkipTest
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = label1\n")
        metad.data.append("scaling_size = [2,2,2]\n")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/POSCAR_perfect").structure
        myperf = MASTFile("files/POSCAR_perfect")
        myperf.to_file("%s/POSCAR" % ingdir)
        myri = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myri.run_scale()
        my_scaled = Poscar.from_file("%s/CONTCAR" % ingdir).structure.get_sorted_structure()
        scaled_compare = Poscar.from_file("files/POSCAR_scaled").structure.get_sorted_structure()
        self.assertEquals(my_scaled, scaled_compare)
