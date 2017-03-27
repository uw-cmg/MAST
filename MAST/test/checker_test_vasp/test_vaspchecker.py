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
from pymatgen.io.vasp import Poscar, Outcar, Potcar, Incar, Kpoints
import numpy as np
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import MASTFile
testname ="checker_test_vasp"
#oldcontrol = os.getenv("MAST_CONTROL")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldscratch
testdir = dirutil.get_test_dir(testname)


class TestVaspChecker(unittest.TestCase):
    """Test Vasp Checker
    """
    def setUp(self):
        os.chdir(testdir)
        if not os.path.exists("childdir"):
            os.mkdir("childdir")
        shutil.copy("files/metadata.txt","childdir")

    def tearDown(self):
        for fname in ["POSCAR","XDATCAR","DYNMAT","OSZICAR","DYNMAT_combined","KPOINTS","POTCAR","INCAR","WAVECAR","CHGCAR","POSCAR_no_sd","XDATCAR_combined","CONTCAR","metadata.txt"]:
            try:
                os.remove("childdir/%s" % fname)
            except OSError:
                pass
        os.rmdir("childdir")

    def test_get_structure_from_file(self):
        myvc = VaspChecker(name = "structure")
        getstr = myvc.get_structure_from_file("structure/POSCAR_perfect")
        mystr = Poscar.from_file("structure/POSCAR_perfect").structure
        self.assertEqual(mystr, getstr)

    def test_get_initial_structure(self):
        myvc = VaspChecker(name = "structure")
        nodir = myvc.get_initial_structure_from_directory()
        compare_nodir = Poscar.from_file("structure/POSCAR").structure
        self.assertEqual(nodir, compare_nodir)
        withdir = myvc.get_initial_structure_from_directory("withdir")
        compare_withdir = Poscar.from_file("withdir/POSCAR").structure
        self.assertEqual(withdir, compare_withdir)

    def test_get_final_structure(self):
        myvc = VaspChecker(name = "structure")
        nodir = myvc.get_final_structure_from_directory()
        compare_nodir = Poscar.from_file("structure/CONTCAR").structure
        self.assertEqual(nodir, compare_nodir)
        withdir = myvc.get_final_structure_from_directory("withdir")
        compare_withdir = Poscar.from_file("withdir/CONTCAR").structure
        self.assertEqual(withdir, compare_withdir)
    def test_forward_final_structure(self):
        myvc = VaspChecker(name = "structure")
        myvc.forward_final_structure_file(os.path.join(testdir,"childdir"))
        fstrp = Poscar.from_file("structure/CONTCAR").structure
        fstrc = Poscar.from_file("childdir/POSCAR").structure
        self.assertEqual(fstrp, fstrc)

    def test_forward_initial_structure(self):
        myvc = VaspChecker(name = "structure")
        myvc.forward_initial_structure_file(os.path.join(testdir,"childdir"))
        istrp = Poscar.from_file("structure/POSCAR").structure
        istrc = Poscar.from_file("childdir/POSCAR").structure
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


    def test_forward_energy_file(self):
        myvc = VaspChecker(name="energy")
        myvc.forward_energy_file(os.path.join(testdir, "childdir"))
        op = Outcar("energy/OSZICAR")
        oc = Outcar("childdir/OSZICAR")
        self.assertEqual(op.run_stats,oc.run_stats)
    def test_is_complete(self):
        vcc = VaspChecker(name="done")
        vcs = VaspChecker(name="started")
        self.assertFalse(vcs.is_complete())
        self.assertTrue(vcc.is_complete())

    def test_is_complete_additional_cases(self):
        kdict=dict()
        kdict['ibrion'] = "-1"
        kdict['nsw'] = "0"
        vcs = VaspChecker(name="done_static",program_keys=kdict)
        self.assertTrue(vcs.is_complete())
        vcs = VaspChecker(name="done_static_butnotconverged",program_keys=kdict)
        self.assertFalse(vcs.is_complete())
        vcs = VaspChecker(name="notdone_static_notconverged",program_keys=kdict)
        self.assertFalse(vcs.is_complete())
        kdict=dict()
        kdict['ibrion'] = "2"
        kdict['nsw'] = "191"
        vcr = VaspChecker(name="notdone_notconverged",program_keys=kdict)
        self.assertFalse(vcr.is_complete())
        vcr = VaspChecker(name="done_butnotconverged",program_keys=kdict)
        self.assertFalse(vcr.is_complete())
        kdict=dict()
        kdict['ibrion'] = "0"
        kdict['nsw'] = "191"
        vcr = VaspChecker(name="done_butnotconverged",program_keys=kdict)
        self.assertTrue(vcr.is_complete())
        vcr = VaspChecker(name="notdone_notconverged",program_keys=kdict)
        self.assertFalse(vcr.is_complete())
        kdict=dict()
        kdict['ibrion'] = "5"
        kdict['nsw'] = "191"
        vcr = VaspChecker(name="done_butnotconverged",program_keys=kdict)
        self.assertTrue(vcr.is_complete())
        vcr = VaspChecker(name="notdone_notconverged",program_keys=kdict)
        self.assertFalse(vcr.is_complete())

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
        myvc.combine_dynamical_matrix_files()
        shutil.move(os.path.join(testdir, "dynamics_split/DYNMAT"),os.path.join(testdir,"childdir"))
        shutil.move(os.path.join(testdir, "dynamics_split/DYNMAT_combined"),os.path.join(testdir,"childdir"))
        dynmat_compare = myvc.read_my_dynamical_matrix_file(myvc.keywords['name'],"DYNMAT_compare")
        dynmat_combined = myvc.read_my_dynamical_matrix_file("childdir","DYNMAT_combined")
        self.assertEqual(dynmat_compare, dynmat_combined)

    def test_vasp_poscar_setup_has_poscar(self):
        compare_pos = Poscar.from_file(os.path.join(testdir, "structure","POSCAR_perfect"))
        compare_pos.write_file(os.path.join(testdir,"childdir","POSCAR"))
        myvc = VaspChecker(name="childdir")
        mypos = myvc._vasp_poscar_setup()
        self.assertEqual(mypos.structure, compare_pos.structure)

    def test_vasp_poscar_setup_no_poscar(self):
        compare_pos = Poscar.from_file(os.path.join(testdir, "structure","POSCAR_perfect"))
        myvc = VaspChecker(name="childdir",structure=compare_pos.structure)
        self.assertFalse(os.path.isfile(os.path.join(testdir, "childdir","POSCAR")))
        myvc._vasp_poscar_setup()
        mypos = Poscar.from_file(os.path.join(testdir, "childdir","POSCAR"))
        self.assertEqual(mypos.structure, compare_pos.structure)

    def test_vasp_poscar_setup_mast_coordinates(self):
        kdict=dict()
        kdict['mast_coordinates'] = ["structure/POSCAR_coordinates"] #note that the input is a list of strings
        myvc = VaspChecker(name="childdir", program_keys=kdict)
        perf_pos = Poscar.from_file(os.path.join(testdir,"structure","POSCAR_perfect"))
        perf_pos.write_file(os.path.join(testdir,"childdir","POSCAR"))
        coord_pos = Poscar.from_file(os.path.join(testdir,"structure","POSCAR_coordinates"))
        grafted_pos = Poscar.from_file(os.path.join(testdir,"structure","POSCAR_grafted"))
        mypos = myvc._vasp_poscar_setup()
        self.assertEqual(mypos.structure, grafted_pos.structure)

    def test_vasp_kpoints_setup_from_metafile(self):
        kdict=dict()
        mymeta = MASTFile("childdir/metadata.txt")
        mymeta.data.append("kpoints = 3x1x7 G 0.5 0.2 .1\n")
        mymeta.to_file("childdir/metadata.txt")
        mymeta2=MASTFile("childdir/metadata.txt")
        print mymeta2.data
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        mykpt = myvc._vasp_kpoints_setup()
        kpt_compare = Kpoints.from_file("files/KPOINTS_317G")
        self.assertEqual(kpt_compare.kpts[0][0],mykpt.kpts[0][0])
        self.assertEqual(kpt_compare.kpts[0][1],mykpt.kpts[0][1])
        self.assertEqual(kpt_compare.kpts[0][2],mykpt.kpts[0][2])
        self.assertEqual(kpt_compare.num_kpts, mykpt.num_kpts)
        self.assertEqual(kpt_compare.style, mykpt.style)
        #self.assertEqual(kpt_compare.kpts_shift, mykpt.kpts_shift)
        self.assertEqual((0.5,0.2,0.1),mykpt.kpts_shift)


    def test_vasp_kpoints_setup_from_keyword(self):
        kdict=dict()
        kdict['mast_kpoints'] = [3,3,3,"M"]
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        mykpt=myvc._vasp_kpoints_setup()
        kpt_compare = Kpoints.from_file("files/KPOINTS_333M")
        self.assertEqual(kpt_compare.kpts[0][0],mykpt.kpts[0][0])
        self.assertEqual(kpt_compare.kpts[0][1],mykpt.kpts[0][1])
        self.assertEqual(kpt_compare.kpts[0][2],mykpt.kpts[0][2])
        self.assertEqual(kpt_compare.num_kpts, mykpt.num_kpts)
        self.assertEqual(kpt_compare.style, mykpt.style)
        
        os.remove("childdir/KPOINTS")
        kdict['mast_kpoints'] = [1,2,5,"G"]
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        mykpt=myvc._vasp_kpoints_setup()
        kpt_compare = Kpoints.from_file("files/KPOINTS_125G")
        self.assertEqual(kpt_compare.kpts[0][0],mykpt.kpts[0][0])
        self.assertEqual(kpt_compare.kpts[0][1],mykpt.kpts[0][1])
        self.assertEqual(kpt_compare.kpts[0][2],mykpt.kpts[0][2])
        self.assertEqual(kpt_compare.num_kpts, mykpt.num_kpts)
        self.assertEqual(kpt_compare.style, mykpt.style)

        os.remove("childdir/KPOINTS")
        kdict['mast_kpoints'] = [2,2,2,"throw error"]
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        self.assertRaises(MASTError, myvc._vasp_kpoints_setup)

    def test_vasp_potcar_setup(self):
        kdict=dict()
        kdict['mast_xc'] = "pbe"
        my_poscar = Poscar.from_file("structure/POSCAR_Al")
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        mypotcar=myvc._vasp_potcar_setup(my_poscar)
        mypotcar = Potcar.from_file("childdir/POTCAR")
        potcar_compare = Potcar.from_file("files/POTCAR_Al_PBE")
        self.assertEqual(potcar_compare[0],mypotcar[0])
        
        kdict['mast_xc'] = "pw91"
        kdict['mast_pp_setup']={'La':'La','Ni':'Ni_pv','O':'O_s'}
        my_poscar = Poscar.from_file("structure/POSCAR_LNO")
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        mypotcar=myvc._vasp_potcar_setup(my_poscar)
        mypotcar = Potcar.from_file("childdir/POTCAR")
        potcar_compare = Potcar.from_file("files/POTCAR_LNO_PW91")
        self.assertEqual(potcar_compare[0],mypotcar[0])

    def test_vasp_incar_setup(self):
        my_structure = Poscar.from_file("structure/POSCAR_LNO").structure
        kdict=dict()
        kdict['mast_xc'] = "pw91"
        kdict['mast_pp_setup']={'La':'La','Ni':'Ni_pv','O':'O_s'}
        kdict['IBRION'] = 3
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        mypos = myvc._vasp_poscar_setup()
        mypot = myvc._vasp_potcar_setup(mypos)
        myvc._vasp_incar_setup(mypot, mypos)
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR_notags"))
        self.assertEqual(myincar, incar_compare)
        #
        kdict['ENCUT'] = 100
        kdict['mast_multiplyencut'] = 1.1
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        myvc._vasp_incar_setup(mypot, mypos)
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR_encut"))
        self.assertEqual(myincar, incar_compare)
        #
        kdict.pop("ENCUT")
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        myvc._vasp_incar_setup(mypot, mypos)
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR_multiply_110"))
        self.assertEqual(myincar, incar_compare)
        #
        kdict.pop("mast_multiplyencut")
        kdict["mast_charge"] = -2
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        myvc._vasp_incar_setup(mypot, mypos)
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR_charge_neg"))
        self.assertEqual(myincar, incar_compare)
        #
        kdict["mast_charge"] = 2
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        myvc._vasp_incar_setup(mypot, mypos)
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR_charge_pos"))
        self.assertEqual(myincar, incar_compare)
        #
        kdict.pop("mast_charge")
        kdict["mast_setmagmom"]="3 5 2"
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        myvc._vasp_incar_setup(mypot, mypos)
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR_setmagmom"))
        self.assertEqual(myincar, incar_compare)
        #
        kdict["mast_setmagmom"]="-1 -2 3 1 2 3 1 2 7 8 7 8 7 8 7 8 5 4 3 2 1 1 2 3 4 5 -5 -4 -3 -2 -1 1 2 3 4 5 4 4 4 -4"
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        myvc._vasp_incar_setup(mypot, mypos)
        myincar = Incar.from_file(os.path.join(testdir,"childdir","INCAR"))
        incar_compare = Incar.from_file(os.path.join(testdir,"files","INCAR_setmagmom_indiv"))
        self.assertEqual(myincar, incar_compare)

    def test_set_up_program_input(self):
        my_structure = Poscar.from_file("structure/POSCAR_Al").structure
        kdict=dict()
        kdict['IBRION']="2"
        kdict['ISIF']="3"
        kdict['ISMEAR']="1"
        kdict['LCHARG']="FALSE"
        kdict['LWAVE']="FALSE"
        kdict['NSW']="191"
        kdict['PREC']="Accurate"
        kdict['SIGMA']="0.2"
        kdict['mast_kpoints']=[1,1,1,"G"]
        kdict['mast_xc']="pw91"
        myvc = VaspChecker(name="childdir",program_keys=kdict,structure=my_structure)
        myvc.set_up_program_input()
        myincar = Incar.from_file("childdir/INCAR")
        myposcar = Poscar.from_file("childdir/POSCAR")
        mypotcar = Potcar.from_file("childdir/POTCAR")
        mykpts = Kpoints.from_file("childdir/KPOINTS")
        compare_incar = Incar.from_file("ready/INCAR")
        compare_poscar = Poscar.from_file("ready/POSCAR")
        compare_potcar = Potcar.from_file("ready/POTCAR")
        compare_kpoints = Kpoints.from_file("ready/KPOINTS")
        self.assertEqual(myincar, compare_incar)
        self.assertEqual(myposcar.structure, compare_poscar.structure)
        self.assertEqual(mykpts.kpts[0][0],compare_kpoints.kpts[0][0])
        self.assertEqual(mykpts.kpts[0][1],compare_kpoints.kpts[0][1])
        self.assertEqual(mykpts.kpts[0][2],compare_kpoints.kpts[0][2])
        self.assertEqual(mykpts.num_kpts, compare_kpoints.num_kpts)
        self.assertEqual(mykpts.style, compare_kpoints.style)


        self.assertEqual(mypotcar, compare_potcar)

    #def test_forward_extra_restart_files(self):
    #    myvc = VaspChecker(name="files")
    #    myvc.forward_extra_restart_files("childdir")
    #    myfiles=dirutil.walkfiles("childdir")
    #    self.assertTrue("childdir/WAVECAR" in myfiles)
    #    self.assertTrue("childdir/CHGCAR" in myfiles)

    def test_softlink_charge_density_file(self):
        import subprocess
        myvc = VaspChecker(name="files")
        myvc.softlink_charge_density_file("childdir")
        myfiles=dirutil.walkfiles("childdir")
        print myfiles
        #self.assertTrue("childdir/CHGCAR" in myfiles)
        listme=subprocess.Popen("ls -l childdir",stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        listres=listme.communicate()[0]
        listme.wait()
        print listres
        self.assertTrue("files/CHGCAR" in listres)
    def test_softlink_wavefunction_file(self):
        import subprocess
        myvc = VaspChecker(name="files")
        myvc.softlink_wavefunction_file("childdir")
        myfiles=dirutil.walkfiles("childdir")
        print myfiles
        #self.assertTrue("childdir/CHGCAR" in myfiles)
        listme=subprocess.Popen("ls -l childdir",stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        listres=listme.communicate()[0]
        listme.wait()
        print listres
        self.assertTrue("files/WAVECAR" in listres)

    def test_forward_charge_density_file(self):
        myvc = VaspChecker(name="files")
        myvc.forward_charge_density_file("childdir")
        self.assertTrue(os.path.isfile("childdir/CHGCAR"))
    def test_forward_wavefunction_file(self):
        myvc = VaspChecker(name="files")
        myvc.forward_wavefunction_file("childdir")
        self.assertTrue(os.path.isfile("childdir/WAVECAR"))
    def test_add_selective_dynamics(self):
        mystruc = Poscar.from_file(os.path.join(testdir,"structure","POSCAR_Al")).structure
        myvc = VaspChecker(name="childdir",structure=mystruc)
        myvc._vasp_poscar_setup()
        sdarr = np.zeros([4,3],bool)
        sdarr[0][0]=True
        sdarr[1][1]=True
        sdarr[2][2]=True
        sdarr[3][0]=True
        sdarr[3][1]=True
        myvc.add_selective_dynamics_to_structure_file(sdarr)
        after_pos = Poscar.from_file(os.path.join(testdir,"childdir","POSCAR"))
        compare_pos = Poscar.from_file(os.path.join(testdir,"structure","POSCAR_Al_sd"))
        self.assertEqual(after_pos.structure, compare_pos.structure)
        self.assertEqual(after_pos.selective_dynamics, compare_pos.selective_dynamics)

    def test_get_non_mast_keywords(self):
        kdict=dict()
        kdict['mast_nokeyword']="hello"
        kdict['mast_charge']="1"
        kdict['IBRION']="3"
        kdict['skip_keyword']="skip"
        myvc = VaspChecker(name="childdir",program_keys=kdict)
        keydict=myvc._vasp_incar_get_non_mast_keywords()
        self.assertEqual(keydict,{"IBRION":"3"})

    def test_get_energy(self):
        myvc = VaspChecker(name="done")
        myenergy=myvc.get_energy()
        self.assertAlmostEqual(myenergy, float(-.83312666E+01),places=7)

    def test_get_max_enmax_from_potcar(self):
        mypotcar = Potcar.from_file(os.path.join(testdir,"files","POTCAR_LNO_PW91"))
        myvc = VaspChecker(name="childdir")
        maxenmax = myvc._get_max_enmax_from_potcar(mypotcar)
        self.assertEqual(maxenmax, 367.921)
        mypotcar = Potcar.from_file(os.path.join(testdir,"files","POTCAR_LTO_PW91"))
        maxenmax = myvc._get_max_enmax_from_potcar(mypotcar)
        self.assertEqual(maxenmax, 400)

    def test_read_my_dynamical_matrix_file(self):
        myvc = VaspChecker(name="dynamics")
        mydm=myvc.read_my_dynamical_matrix_file()
        self.assertEqual(mydm['numspec'],1)
        self.assertEqual(mydm['numatoms'],4)
        self.assertEqual(mydm['numdisp'],3)
        self.assertEqual(mydm['massline']," 55.847\n")
        self.assertEqual(mydm['atoms'][1][1]['displine'],"0.0100 0.0000 0.0000")
        self.assertEqual(mydm['atoms'][1][3]['displine'],"0.0000 0.0000 0.0100")
        self.assertEqual(mydm['atoms'][1][1]['dynmat'][0]," -0.101036   0.000000   0.000000\n")
        self.assertEqual(mydm['atoms'][1][3]['dynmat'][2],"  0.000000   0.000000   0.000000\n")

    def test_write_my_dynamical_matrix_file(self):
        myvc = VaspChecker(name="childdir")
        mydm=myvc.read_my_dynamical_matrix_file("dynamics")
        myvc.write_my_dynamical_matrix_file(mydm)
        mydm2 = myvc.read_my_dynamical_matrix_file()
        self.assertEqual(mydm,mydm2)
    def test_write_my_dynamical_matrix_file_without_disp_mass(self):
        myvc = VaspChecker(name="childdir")
        mydm=myvc.read_my_dynamical_matrix_file("dynamics")
        myvc.write_my_dynmat_without_disp_or_mass(mydm)
        myread = MAST.utility.MASTFile("dynamics/DYNMAT_for_PHON")
        myread2 = MAST.utility.MASTFile("childdir/DYNMAT")
        self.assertEqual(myread2.data,myread.data)

    def test_read_my_displacement_file(self):
        myvc = VaspChecker(name="dynamics")
        mydisp = myvc.read_my_displacement_file(myvc.keywords['name'],"XDATCAR_vasp522")
        self.assertEqual(mydisp['numatoms'],4)
        self.assertEqual(mydisp['configs'][6][0],"   0.00000000  0.00000000  0.00286262\n")
        mydisp = myvc.read_my_displacement_file(myvc.keywords['name'],"XDATCAR_vasp5211")
        self.assertEqual(mydisp['numatoms'],4)
        self.assertEqual(mydisp['configs'][2][1],"   0.50244722  0.50000000  0.00000000\n")

    def test_write_my_displacement_file(self):
        myvc = VaspChecker(name="childdir")
        mydisp = myvc.read_my_displacement_file("dynamics","XDATCAR_vasp5211")
        myvc.write_my_displacement_file(mydisp)
        mydisp2 = myvc.read_my_displacement_file()
        self.assertEqual(mydisp,mydisp2)
        mydisp = myvc.read_my_displacement_file("dynamics","XDATCAR_vasp522")
        myvc.write_my_displacement_file(mydisp)
        mydisp2 = myvc.read_my_displacement_file()
        self.assertEqual(mydisp,mydisp2)

    def test_combine_displacement_files(self):
        myvc = VaspChecker(name="dynamics_split")
        myvc.combine_displacement_files()
        shutil.move(os.path.join(testdir, "dynamics_split/XDATCAR"),os.path.join(testdir,"childdir"))
        shutil.move(os.path.join(testdir, "dynamics_split/XDATCAR_combined"),os.path.join(testdir,"childdir"))
        disp_compare = myvc.read_my_displacement_file(myvc.keywords['name'],"XDATCAR_compare")
        disp_combined = myvc.read_my_displacement_file("childdir","XDATCAR_combined")
        print "COMPARE"
        for key,value in disp_compare.iteritems():
            self.assertEqual(value, disp_combined[key])
    def test_get_total_electrons(self):
        mypos = Poscar.from_file("structure/POSCAR_LTO")
        mypot = Potcar.from_file("files/POTCAR_LTO_PW91")
        myvc = VaspChecker()
        numelec = myvc.get_total_electrons(mypos, mypot)
        self.assertEqual(numelec,328)
    def test_get_valence_list(self):
        mypot = Potcar.from_file("files/POTCAR_LTO_PW91")
        myvc = VaspChecker()
        vlist = myvc.get_valence_list(mypot)
        self.assertEqual(vlist,[11,12,6])
    def test_get_energy_from_energy_file(self):
        myvc = VaspChecker(name="done")
        myenergy=myvc.get_energy_from_energy_file()
        self.assertAlmostEqual(myenergy,float(-.83312666E+01),places=7)
    def test_is_started(self):
        myvc = VaspChecker(name="notready1")
        self.assertFalse(myvc.is_started())
        myvc = VaspChecker(name="ready")
        self.assertFalse(myvc.is_started())
        myvc = VaspChecker(name="started")
        self.assertTrue(myvc.is_started())
        myvc = VaspChecker(name="done")
        self.assertTrue(myvc.is_started())
    def test_write_final_structure_file(self):
        myvc = VaspChecker(name="childdir")
        finalstruc = Poscar.from_file("structure/CONTCAR").structure
        myvc.write_final_structure_file(finalstruc)
        mystruc = Poscar.from_file("childdir/CONTCAR").structure
        self.assertEqual(mystruc,finalstruc)
    def test_has_starting_structure_file(self):
        myvc = VaspChecker(name="notready3")
        self.assertFalse(myvc.has_starting_structure_file())
        myvc = VaspChecker(name="ready")
        self.assertTrue(myvc.has_starting_structure_file())
        myvc = VaspChecker(name="notready1")
        self.assertTrue(myvc.has_starting_structure_file())
        myvc = VaspChecker(name="notready2")
        self.assertTrue(myvc.has_starting_structure_file())
        myvc = VaspChecker(name="notready4")
        self.assertTrue(myvc.has_starting_structure_file())
        myvc = VaspChecker(name="notready5")
        self.assertTrue(myvc.has_starting_structure_file())
    def test_has_ending_structure_file(self):
        myvc = VaspChecker(name="ready")
        self.assertFalse(myvc.has_ending_structure_file())
        myvc = VaspChecker(name="done")
        self.assertTrue(myvc.has_ending_structure_file())
    def test_scale_mesh(self):
        myvc = VaspChecker()
        mystr = Poscar.from_file("ready/POSCAR").structure
        myk=myvc.scale_mesh(mystr, 2000)
        self.assertEqual(myvc.keywords['program_keys']['mast_kpoints'],[8,8,8,'M'])
