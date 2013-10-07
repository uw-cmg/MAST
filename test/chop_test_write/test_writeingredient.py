"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import WriteIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile

testname="chop_test_write"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestWriteIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_no_setup(self):
        raise SkipTest
        #self.testclass.no_setup()

    def test_write_neb_from_endpoints_only(self):
        ingdir="writedir/neb_labelinit-labelfin"
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        peinit = MASTFile("files/parent_energy_labelinit")
        peinit.to_file("%s/parent_energy_labelinit" % ingdir)
        pefin = MASTFile("files/parent_energy_labelfin")
        pefin.to_file("%s/parent_energy_labelfin" % ingdir)
        psinit = MASTFile("files/parent_structure_labelinit")
        psinit.to_file("%s/parent_structure_labelinit" % ingdir)
        psfin = MASTFile("files/parent_structure_labelfin")
        psfin.to_file("%s/parent_structure_labelfin" % ingdir)
        kdict=dict()
        kdict['images']=3
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp_neb'
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['neblines']=dict()
        kdict['neblines']['labelinit-labelfin']=neblines
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = WriteIngredient(name=ingdir,program="vasp_neb",program_keys=kdict,structure=my_structure)
        mywi.write_neb()
        #self.testclass.write_neb()
    def test_write_neb_with_parent_image_structures(self):
        raise SkipTest
        metad = MASTFile("files/metadata_neb")
        metad.to_file("writedir/metadata.txt")
        peinit = MASTFile("files/parent_energy_labelinit")
        peinit.to_file("writedir/parent_energy_labelinit")
        pefin = MASTFile("files/parent_energy_labelfin")
        pefin.to_file("writedir/parent_energy_labelfin")
        psinit = MASTFile("files/parent_structure_labelinit")
        psinit.to_file("writedir/parent_structure_labelinit")
        psfin = MASTFile("files/parent_structure_labelfin")
        psfin.to_file("writedir/parent_structure_labelfin")
        psim1 = MASTFile("files/parent_structure_labelinit-labelfin_01")
        psim1.to_file("writedir/parent_structure_labelinit-labelfin_01")
        psim2 = MASTFile("files/parent_structure_labelinit-labelfin_02")
        psim2.to_file("writedir/parent_structure_labelinit-labelfin_02")
        psim3 = MASTFile("files/parent_structure_labelinit-labelfin_03")
        psim3.to_file("writedir/parent_structure_labelinit-labelfin_03")
        kdict=dict()
        kdict['images']=3
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = WriteIngredient(name="writedir/neb_labelinit-labelfin",program="vasp",program_keys=kdict,structure=my_structure)
        raise SkipTest

    def test_get_parent_structures(self):
        raise SkipTest
        #self.testclass.get_parent_structures()

    def test_get_parent_image_structures(self):
        raise SkipTest
        #self.testclass.get_parent_image_structures()

    def test_place_parent_energy_files(self):
        raise SkipTest
        #self.testclass.place_parent_energy_files()

    def test_write_neb_subfolders(self):
        raise SkipTest
        #self.testclass.write_neb_subfolders()

    def test_write_singlerun(self):
        raise SkipTest
        #self.testclass.write_singlerun()

    def test_write_phonon_multiple(self):
        raise SkipTest
        #self.testclass.write_phonon_multiple()

    def test_write_phonon_single(self):
        raise SkipTest
        #self.testclass.write_phonon_single()

    def test_get_my_phonon_params(self):
        raise SkipTest
        #self.testclass.get_my_phonon_params()

