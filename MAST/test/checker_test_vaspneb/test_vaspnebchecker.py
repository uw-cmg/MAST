"""Tests for Vaspnebchecker"""

from MAST.ingredients.checker.vaspnebchecker import VaspNEBChecker
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
testname="checker_test_vaspneb"
testdir = dirutil.get_test_dir(testname)

class TestVaspnebchecker(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir('childdir'):
            os.mkdir('childdir')
        shutil.copy("files/metadata.txt","childdir")

    def tearDown(self):
        #return
        dnames = dirutil.walkdirs("childdir",1,1)
        for dname in dnames:
            shutil.rmtree(dname)
        fnames = dirutil.walkfiles("childdir",1,1)
        for fname in fnames:
            os.remove(fname)


    def test_get_path_to_write_neb_parent_energy(self):
        kdict=dict()
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        myvcneb=VaspNEBChecker(name="childdir",program_keys = kdict)
        mypath = myvcneb.get_path_to_write_neb_parent_energy(1)
        self.assertEqual(mypath, "childdir/00/OSZICAR")
        mypath = myvcneb.get_path_to_write_neb_parent_energy(2)
        self.assertEqual(mypath, "childdir/04/OSZICAR")
        mypath = myvcneb.get_path_to_write_neb_parent_energy('02')
        self.assertEqual(mypath, "childdir/02/OSZICAR")
        #self.testclass.get_path_to_write_neb_parent_energy(myimages, parent)

    def test_set_up_neb_folders_no_mast_coordinates(self):
        mystrs=list()
        pos=dict()
        for posstr in ['00','01','02','03','04']:
            pos[posstr] = Poscar.from_file("structures/POSCAR_%s" % posstr)
            mystrs.append(pos[posstr].structure)
        kdict=dict()
        kdict['images']=3
        myvcneb=VaspNEBChecker(name="childdir",program_keys=kdict)
        myvcneb.set_up_neb_folders(mystrs)
        for subdir in ['00','01','02','03','04']:
            mypos = Poscar.from_file("childdir/POSCAR_%s" % subdir)
            self.assertEqual(mypos.structure,pos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,pos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,pos[subdir].structure.sites)
            mypos = Poscar.from_file("childdir/%s/POSCAR" % subdir)
            self.assertEqual(mypos.structure,pos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,pos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,pos[subdir].structure.sites)
        #self.testclass.set_up_neb_folders(image_structures)
    def test_set_up_neb_folders_with_mast_coordinates(self):
        kdict=dict()
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['mast_coordinates']=["structures/POSCAR_coords_01","structures/POSCAR_coords_02","structures/POSCAR_coords_03"]
        mystrs=list()
        pos=dict()
        graftedpos=dict()
        for posstr in ['00','01','02','03','04']:
            graftedpos[posstr] = Poscar.from_file("structures/POSCAR_grafted_%s" % posstr)
            pos[posstr] = Poscar.from_file("structures/POSCAR_%s" % posstr)
            mystrs.append(pos[posstr].structure)
        myvcneb=VaspNEBChecker(name="childdir",program_keys=kdict)
        myvcneb.set_up_neb_folders(mystrs)
        for subdir in ['00','01','02','03','04']:
            mypos = Poscar.from_file("childdir/POSCAR_%s" % subdir)
            self.assertEqual(mypos.structure,graftedpos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,graftedpos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,graftedpos[subdir].structure.sites)
            mypos = Poscar.from_file("childdir/%s/POSCAR" % subdir)
            self.assertEqual(mypos.structure,graftedpos[subdir].structure)
            self.assertEqual(mypos.structure.lattice,graftedpos[subdir].structure.lattice)
            self.assertEqual(mypos.structure.sites,graftedpos[subdir].structure.sites)
        #self.testclass.set_up_neb_folders(image_structures)

    def test_is_complete(self):
        kdict=dict()
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        myvcneb=VaspNEBChecker(name="notready1",program_keys=kdict)
        self.assertFalse(myvcneb.is_complete())
        myvcneb=VaspNEBChecker(name="ready",program_keys=kdict)
        self.assertFalse(myvcneb.is_complete())
        myvcneb=VaspNEBChecker(name="started",program_keys=kdict)
        self.assertFalse(myvcneb.is_complete())
        myvcneb=VaspNEBChecker(name="done",program_keys=kdict)
        self.assertTrue(myvcneb.is_complete())
        #self.testclass.is_complete()

    def test_is_ready_to_run(self):
        kdict=dict()
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        myvcneb=VaspNEBChecker(name="notready1",program_keys=kdict)
        self.assertFalse(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="notready2",program_keys=kdict)
        self.assertFalse(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="notready3",program_keys=kdict)
        self.assertFalse(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="notready4",program_keys=kdict)
        self.assertFalse(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="notready5",program_keys=kdict)
        self.assertFalse(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="notready6",program_keys=kdict)
        self.assertFalse(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="ready", program_keys=kdict)
        self.assertTrue(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="started", program_keys=kdict)
        self.assertTrue(myvcneb.is_ready_to_run())
        myvcneb=VaspNEBChecker(name="done", program_keys=kdict)
        self.assertTrue(myvcneb.is_ready_to_run())
        
        #self.testclass.is_ready_to_run()

    def test__vasp_incar_setup(self):
        my_poscar = Poscar.from_file("structures/POSCAR_00")
        my_structure = my_poscar.structure
        my_poscar.write_file("childdir/POSCAR")
        kdict=dict()
        kdict['mast_xc'] = "pw91"
        kdict['mast_pp_setup']={'Cr':'Cr_pv','Fe':'Fe_sv','Ni':'Ni_pv'}
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['IBRION'] = '1'
        kdict['POTIM'] = '0.5'
        kdict['LCLIMB'] = 'True'
        kdict['SPRING'] ='-5'
        myvc = VaspNEBChecker(name="childdir",program_keys=kdict,structure=my_structure)
        mypot = myvc._vasp_potcar_setup(my_poscar)
        myvc._vasp_incar_setup(mypot, my_poscar)
        myvc._vasp_neb_incar_modify()
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR"))
        self.assertEqual(myincar, incar_compare)
        #self.testclass._vasp_incar_setup(my_potcar, my_poscar)

    def test_set_up_program_input(self):
        mystrs=list()
        pos=dict()
        for posstr in ['00','01','02','03','04']:
            pos[posstr] = Poscar.from_file("structures/POSCAR_%s" % posstr)
            mystrs.append(pos[posstr].structure)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']="pw91"
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['mast_neb_settings']['image_structures'] = mystrs
        myvcneb=VaspNEBChecker(name="childdir",program_keys=kdict)
        myvcneb.set_up_program_input()
        submitscript = MASTFile()
        submitscript.data="Submission script placeholder."
        submitscript.to_file("childdir/submit.sh")
        ozholder = MASTFile()
        ozholder.data="OSZICAR placeholder."
        submitscript.to_file("childdir/00/OSZICAR")
        submitscript.to_file("childdir/04/OSZICAR")
        self.assertTrue(myvcneb.is_ready_to_run())
        #self.testclass.set_up_program_input_neb(image_structures)

    def test_get_energy_from_energy_file(self):
        kdict=dict()
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        myvcneb=VaspNEBChecker(name="done",program_keys=kdict)
        estr = myvcneb.get_energy_from_energy_file()
        estr_compare="-99.860;-99.649;-99.538;-99.649;-99.860"
        self.assertEqual(estr,estr_compare)
        #self.testclass.get_energy_from_energy_file()

    def test_is_started(self):
        kdict=dict()
        kdict['images']=3
        myvcneb=VaspNEBChecker(name="notready1",program_keys=kdict)
        self.assertFalse(myvcneb.is_started())
        myvcneb=VaspNEBChecker(name="ready",program_keys=kdict)
        self.assertFalse(myvcneb.is_started())
        myvcneb=VaspNEBChecker(name="started",program_keys=kdict)
        self.assertTrue(myvcneb.is_started())
        myvcneb=VaspNEBChecker(name="done",program_keys=kdict)
        self.assertTrue(myvcneb.is_started())
        #self.testclass.is_started()

