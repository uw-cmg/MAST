"""Tests for Vaspnebchecker"""

from MAST.ingredients.checker.vaspnebchecker import VaspNEBChecker
import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
import shutil
testname="checker_test_vaspneb"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestVaspnebchecker(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        #return
        dnames = dirutil.walkdirs("childdir",1,1)
        for dname in dnames:
            shutil.rmtree(dname)
        fnames = dirutil.walkfiles("childdir",1,1)
        for fname in fnames:
            os.remove(fname)

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_get_path_to_write_neb_parent_energy(self):
        myvcneb=VaspNEBChecker(name="childdir")
        mypath = myvcneb.get_path_to_write_neb_parent_energy(3,1)
        self.assertEqual(mypath, "childdir/00/OSZICAR")
        mypath = myvcneb.get_path_to_write_neb_parent_energy(3,2)
        self.assertEqual(mypath, "childdir/04/OSZICAR")
        mypath = myvcneb.get_path_to_write_neb_parent_energy(3,'02')
        self.assertEqual(mypath, "childdir/02/OSZICAR")
        #self.testclass.get_path_to_write_neb_parent_energy(myimages, parent)

    def test_set_up_neb_folders_no_mast_coordinates(self):
        mystrs=list()
        pos=dict()
        for posstr in ['00','01','02','03','04']:
            pos[posstr] = pymatgen.io.vaspio.Poscar.from_file("structures/POSCAR_%s" % posstr)
            mystrs.append(pos[posstr].structure)
        kdict=dict()
        kdict['images']=3
        myvcneb=VaspNEBChecker(name="childdir",program_keys=kdict)
        myvcneb.set_up_neb_folders(mystrs)
        for subdir in ['00','01','02','03','04']:
            mypos = pymatgen.io.vaspio.Poscar.from_file("childdir/POSCAR_%s" % subdir)
            self.assertEqual(mypos.structure,pos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,pos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,pos[subdir].structure.sites)
            mypos = pymatgen.io.vaspio.Poscar.from_file("childdir/%s/POSCAR" % subdir)
            self.assertEqual(mypos.structure,pos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,pos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,pos[subdir].structure.sites)
        #self.testclass.set_up_neb_folders(image_structures)
    def test_set_up_neb_folders_with_mast_coordinates(self):
        kdict=dict()
        kdict['images']=3
        kdict['mast_coordinates']=["structures/POSCAR_coords_01","structures/POSCAR_coords_02","structures/POSCAR_coords_03"]
        mystrs=list()
        pos=dict()
        graftedpos=dict()
        for posstr in ['00','01','02','03','04']:
            graftedpos[posstr] = pymatgen.io.vaspio.Poscar.from_file("structures/POSCAR_grafted_%s" % posstr)
            pos[posstr] = pymatgen.io.vaspio.Poscar.from_file("structures/POSCAR_%s" % posstr)
            mystrs.append(pos[posstr].structure)
        myvcneb=VaspNEBChecker(name="childdir",program_keys=kdict)
        myvcneb.set_up_neb_folders(mystrs)
        for subdir in ['00','01','02','03','04']:
            mypos = pymatgen.io.vaspio.Poscar.from_file("childdir/POSCAR_%s" % subdir)
            self.assertEqual(mypos.structure,graftedpos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,graftedpos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,graftedpos[subdir].structure.sites)
            mypos = pymatgen.io.vaspio.Poscar.from_file("childdir/%s/POSCAR" % subdir)
            self.assertEqual(mypos.structure,graftedpos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,graftedpos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,graftedpos[subdir].structure.sites)
        #self.testclass.set_up_neb_folders(image_structures)

    def test_is_complete(self):
        raise SkipTest
        #self.testclass.is_complete()

    def test_is_ready_to_run(self):
        raise SkipTest
        #self.testclass.is_ready_to_run()

    def test__vasp_incar_setup(self):
        raise SkipTest
        #self.testclass._vasp_incar_setup(my_potcar, my_poscar)

    def test_set_up_program_input_neb(self):
        raise SkipTest
        #self.testclass.set_up_program_input_neb(image_structures)

    def test_get_energy_from_energy_file(self):
        raise SkipTest
        #self.testclass.get_energy_from_energy_file()

    def test_is_started(self):
        raise SkipTest
        #self.testclass.is_started()

