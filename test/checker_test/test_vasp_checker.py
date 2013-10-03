import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
from MAST.ingredients.checker.vaspchecker import VaspChecker
import shutil
import pymatgen
import numpy as np
from MAST.utility import MASTError

testname ="checker_test"
#oldcontrol = os.getenv("MAST_CONTROL")
#oldrecipe = os.getenv("MAST_RECIPE_PATH")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldrecipe
#print oldscratch
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)


class TestVaspChecker(unittest.TestCase):
    """Test Vasp Checker
    """
    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        for fname in ["POSCAR","XDATCAR","DYNMAT","OSZICAR","DYNMAT_combined","KPOINTS"]:
            try:
                os.remove("childdir/%s" % fname)
            except OSError:
                pass

    def test_get_structure_from_file(self):
        myvc = VaspChecker(name = "structure")
        getstr = myvc.get_structure_from_file("structure/POSCAR_perfect")
        mystr = pymatgen.io.vaspio.Poscar.from_file("structure/POSCAR_perfect").structure
        self.assertEqual(mystr, getstr)

    def test_get_initial_structure(self):
        myvc = VaspChecker(name = "structure")
        nodir = myvc.get_initial_structure_from_directory()
        compare_nodir = pymatgen.io.vaspio.Poscar.from_file("structure/POSCAR").structure
        self.assertEqual(nodir, compare_nodir)
        withdir = myvc.get_initial_structure_from_directory("withdir")
        compare_withdir = pymatgen.io.vaspio.Poscar.from_file("withdir/POSCAR").structure
        self.assertEqual(withdir, compare_withdir)

    def test_get_final_structure(self):
        myvc = VaspChecker(name = "structure")
        nodir = myvc.get_final_structure_from_directory()
        compare_nodir = pymatgen.io.vaspio.Poscar.from_file("structure/CONTCAR").structure
        self.assertEqual(nodir, compare_nodir)
        withdir = myvc.get_final_structure_from_directory("withdir")
        compare_withdir = pymatgen.io.vaspio.Poscar.from_file("withdir/CONTCAR").structure
        self.assertEqual(withdir, compare_withdir)
    def test_forward_final_structure(self):
        myvc = VaspChecker(name = "structure")
        myvc.forward_final_structure_file(os.path.join(testdir,"childdir"))
        fstrp = pymatgen.io.vaspio.Poscar.from_file("structure/CONTCAR").structure
        fstrc = pymatgen.io.vaspio.Poscar.from_file("childdir/POSCAR").structure
        self.assertEqual(fstrp, fstrc)

    def test_forward_initial_structure(self):
        myvc = VaspChecker(name = "structure")
        myvc.forward_initial_structure_file(os.path.join(testdir,"childdir"))
        istrp = pymatgen.io.vaspio.Poscar.from_file("structure/POSCAR").structure
        istrc = pymatgen.io.vaspio.Poscar.from_file("childdir/POSCAR").structure
        self.assertEqual(istrp, istrc)
    
    def test_forward_dynamical_matrix(self):
        myvc = VaspChecker(name = "dynamics")
        myvc.forward_dynamical_matrix_file(os.path.join(testdir,"childdir"))
        dp = myvc.read_my_dynamical_matrix_file("dynamics")
        dc = myvc.read_my_dynamical_matrix_file("childdir")
        self.assertEqual(dp,dc)

    def test_forward_displacement(self):
        myvc = VaspChecker(name="dynamics")
        myvc.forward_displacement_file(os.path.join(testdir,"childdir"))
        dp = myvc.read_my_displacement_file("dynamics")
        dc = myvc.read_my_displacement_file("childdir")
        self.assertEqual(dp,dc)

    def test_write_displacement(self):
        myvc = VaspChecker(name="dynamics")
        dp = myvc.read_my_displacement_file("dynamics")
        myvc.write_my_displacement_file("childdir", dp)
        dw = myvc.read_my_displacement_file("childdir")
        self.assertEqual(dp,dw)

    def test_forward_energy_file(self):
        myvc = VaspChecker(name="energy")
        myvc.forward_energy_file(os.path.join(testdir, "childdir"))
        op = pymatgen.io.vaspio.Outcar("energy/OSZICAR")
        oc = pymatgen.io.vaspio.Outcar("childdir/OSZICAR")
        self.assertEqual(op.run_stats,oc.run_stats)
    def test_is_complete(self):
        vcc = VaspChecker(name="done")
        vcs = VaspChecker(name="started")
        self.assertFalse(vcs.is_complete())
        self.assertTrue(vcc.is_complete())

    def test_is_ready(self):
        vcnr1 = VaspChecker(name="notready1")
        vcnr2 = VaspChecker(name="notready2")
        vcnr3 = VaspChecker(name="notready3")
        vcnr4 = VaspChecker(name="notready4")
        vcnr5 = VaspChecker(name="notready5")
        vcr = VaspChecker(name="ready")
        self.assertFalse(vcnr1.is_ready_to_run())
        self.assertFalse(vcnr2.is_ready_to_run())
        self.assertFalse(vcnr3.is_ready_to_run())
        self.assertFalse(vcnr4.is_ready_to_run())
        self.assertFalse(vcnr5.is_ready_to_run())
        self.assertTrue(vcr.is_ready_to_run())

    def test_combine_dynamical_matrix(self):
        myvc = VaspChecker(name="dynamics_split")
        myvc.combine_dynamical_matrix_files(myvc.keywords['name'])
        shutil.move(os.path.join(testdir, "dynamics_split/DYNMAT"),os.path.join(testdir,"childdir"))
        shutil.move(os.path.join(testdir, "dynamics_split/DYNMAT_combined"),os.path.join(testdir,"childdir"))
        dynmat_compare = myvc.read_my_dynamical_matrix_file(myvc.keywords['name'],"DYNMAT_compare")
        dynmat_combined = myvc.read_my_dynamical_matrix_file(os.path.join(testdir, "childdir"),"DYNMAT_combined")
        self.assertEqual(dynmat_compare, dynmat_combined)

    def test_vasp_poscar_setup_has_poscar(self):
        compare_pos = pymatgen.io.vaspio.Poscar.from_file(os.path.join(testdir, "structure","POSCAR_perfect"))
        compare_pos.write_file(os.path.join(testdir,"childdir","POSCAR"))
        myvc = VaspChecker(name="childdir")
        mypos = myvc._vasp_poscar_setup()
        self.assertEqual(mypos.structure, compare_pos.structure)

    def test_vasp_poscar_setup_no_poscar(self):
        compare_pos = pymatgen.io.vaspio.Poscar.from_file(os.path.join(testdir, "structure","POSCAR_perfect"))
        myvc = VaspChecker(name="childdir",structure=compare_pos.structure)
        self.assertFalse(os.path.isfile(os.path.join(testdir, "childdir","POSCAR")))
        myvc._vasp_poscar_setup()
        mypos = pymatgen.io.vaspio.Poscar.from_file(os.path.join(testdir, "childdir","POSCAR"))
        self.assertEqual(mypos.structure, compare_pos.structure)

    def test_vasp_poscar_setup_mast_coordinates(self):
        kdict=dict()
        kdict['mast_coordinates'] = ["structure/POSCAR_coordinates"] #note that the input is a list of strings
        myvc = VaspChecker(name="childdir", program_keys=kdict)
        perf_pos = pymatgen.io.vaspio.Poscar.from_file(os.path.join(testdir,"structure","POSCAR_perfect"))
        perf_pos.write_file(os.path.join(testdir,"childdir","POSCAR"))
        coord_pos = pymatgen.io.vaspio.Poscar.from_file(os.path.join(testdir,"structure","POSCAR_coordinates"))
        grafted_pos = pymatgen.io.vaspio.Poscar.from_file(os.path.join(testdir,"structure","POSCAR_grafted"))
        mypos = myvc._vasp_poscar_setup()
        self.assertEqual(mypos.structure, grafted_pos.structure)

    def test_vasp_kpoints_setup(self):
        kdict=dict()
        kdict['mast_kpoints'] = [3,3,3,"M"]
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        mykpt=myvc._vasp_kpoints_setup()
        kpt_compare = pymatgen.io.vaspio.Kpoints.from_file("files/KPOINTS_333M")
        self.assertEqual(kpt_compare.kpts[0][0],mykpt.kpts[0][0])
        self.assertEqual(kpt_compare.kpts[0][1],mykpt.kpts[0][1])
        self.assertEqual(kpt_compare.kpts[0][2],mykpt.kpts[0][2])
        self.assertEqual(kpt_compare.num_kpts, mykpt.num_kpts)
        self.assertEqual(kpt_compare.style, mykpt.style)

        kdict['mast_kpoints'] = [1,2,5,"G"]
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        mykpt=myvc._vasp_kpoints_setup()
        kpt_compare = pymatgen.io.vaspio.Kpoints.from_file("files/KPOINTS_125G")
        self.assertEqual(kpt_compare.kpts[0][0],mykpt.kpts[0][0])
        self.assertEqual(kpt_compare.kpts[0][1],mykpt.kpts[0][1])
        self.assertEqual(kpt_compare.kpts[0][2],mykpt.kpts[0][2])
        self.assertEqual(kpt_compare.num_kpts, mykpt.num_kpts)
        self.assertEqual(kpt_compare.style, mykpt.style)

        kdict['mast_kpoints'] = [2,2,2,"throw error"]
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        self.assertRaises(MASTError, myvc._vasp_kpoints_setup)
