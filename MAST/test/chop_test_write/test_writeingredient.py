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
testname="chop_test_write"
testdir = dirutil.get_test_dir(testname)

class TestWriteIngredient(unittest.TestCase):
    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir("writedir"):
            os.mkdir("writedir")
        if not os.path.isdir("writedir/neb_labelinit-labelfin"):
            os.mkdir("writedir/neb_labelinit-labelfin")
        if not os.path.isdir("writedir/single_label1"):
            os.mkdir("writedir/single_label1")
        if not os.path.isdir("writedir/single_phonon_label1"):
            os.mkdir("writedir/single_phonon_label1")

    def tearDown(self):
        for foldername in ["writedir/neb_labelinit-labelfin","writedir/single_label1", "writedir/single_phonon_label1"]:
            try:
                shutil.rmtree(foldername)
            except OSError:
                pass

    def test___init__(self):
        self.assertTrue(True)
        #self.testclass.__init__(**kwargs)

    def test_no_setup(self):
        self.assertTrue(True)
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
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp_neb'
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['mast_neb_settings']['lines']=neblines
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        mywi.write_neb()
        self.assertTrue(mywi.checker.is_ready_to_run())
        pos_compare_00 = MASTFile("files/POSCAR_00")
        pos_compare_01 = MASTFile("files/POSCAR_01")
        pos_compare_02 = MASTFile("files/POSCAR_02")
        pos_compare_03 = MASTFile("files/POSCAR_03")
        pos_compare_04 = MASTFile("files/POSCAR_04")
        mypos00 = MASTFile("writedir/neb_labelinit-labelfin/00/POSCAR")
        mypos01 = MASTFile("writedir/neb_labelinit-labelfin/01/POSCAR")
        mypos02 = MASTFile("writedir/neb_labelinit-labelfin/02/POSCAR")
        mypos03 = MASTFile("writedir/neb_labelinit-labelfin/03/POSCAR")
        mypos04 = MASTFile("writedir/neb_labelinit-labelfin/04/POSCAR")
        self.assertEqual(pos_compare_00.data, mypos00.data)
        self.assertEqual(pos_compare_01.data, mypos01.data)
        self.assertEqual(pos_compare_02.data, mypos02.data)
        self.assertEqual(pos_compare_03.data, mypos03.data)
        self.assertEqual(pos_compare_04.data, mypos04.data)
        #self.testclass.write_neb()
    def test_write_neb_with_parent_image_structures(self):
        ingdir="writedir/neb_labelinit-labelfin"
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
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
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp_neb'
        kdict['mast_coordinates']=list()
        kdict['mast_coordinates'].append('%s/files/POSCAR_coords_01' % testdir)
        kdict['mast_coordinates'].append('%s/files/POSCAR_coords_02' % testdir)
        kdict['mast_coordinates'].append('%s/files/POSCAR_coords_03' % testdir)
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['mast_neb_settings']['lines']=neblines
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        mywi.write_neb()
        self.assertTrue(mywi.checker.is_ready_to_run())
        pos_compare_00 = MASTFile("files/POSCAR_grafted_00")
        pos_compare_01 = MASTFile("files/POSCAR_grafted_01")
        pos_compare_02 = MASTFile("files/POSCAR_grafted_02")
        pos_compare_03 = MASTFile("files/POSCAR_grafted_03")
        pos_compare_04 = MASTFile("files/POSCAR_grafted_04")
        mypos00 = MASTFile("writedir/neb_labelinit-labelfin/00/POSCAR")
        mypos01 = MASTFile("writedir/neb_labelinit-labelfin/01/POSCAR")
        mypos02 = MASTFile("writedir/neb_labelinit-labelfin/02/POSCAR")
        mypos03 = MASTFile("writedir/neb_labelinit-labelfin/03/POSCAR")
        mypos04 = MASTFile("writedir/neb_labelinit-labelfin/04/POSCAR")
        self.assertEqual(pos_compare_00.data, mypos00.data)
        self.assertEqual(pos_compare_01.data, mypos01.data)
        self.assertEqual(pos_compare_02.data, mypos02.data)
        self.assertEqual(pos_compare_03.data, mypos03.data)
        self.assertEqual(pos_compare_04.data, mypos04.data)

    def test_get_parent_structures(self):
        kdict=dict()
        kdict['mast_program'] = 'vasp_neb'
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['mast_neb_settings']['lines']=neblines
        ingdir = "writedir/neb_labelinit-labelfin"
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        unsorted_init = MASTFile("unsorted/parent_structure_labelinit")
        #unsorted_init = MASTFile("unsorted/parent_structure_labelinit_scrambled")
        unsorted_init.to_file("%s/parent_structure_labelinit" % ingdir)
        unsorted_fin = MASTFile("unsorted/parent_structure_labelfin")
        unsorted_fin.to_file("%s/parent_structure_labelfin" % ingdir)
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        [sinit, sfin] = mywi.get_parent_structures()
        init_compare = pymatgen.io.vaspio.Poscar.from_file("files/parent_structure_labelinit").structure
        fin_compare = pymatgen.io.vaspio.Poscar.from_file("files/parent_structure_labelfin").structure
        #print sinit
        #print init_compare
        #print sfin
        #print fin_compare
        self.assertEqual(sinit.sites, init_compare.sites)
        self.assertEqual(sinit.lattice, init_compare.lattice)
        self.assertEqual(sfin.sites, fin_compare.sites)
        self.assertEqual(sfin.lattice, fin_compare.lattice)
        #self.testclass.get_parent_structures()

    def test_get_parent_image_structures(self):
        kdict=dict()
        kdict['mast_program'] = 'vasp_neb'
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['mast_neb_settings']['lines']=neblines
        ingdir = "writedir/neb_labelinit-labelfin"
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        unsorted_01 = MASTFile("unsorted/parent_structure_labelinit-labelfin_01")
        unsorted_01.to_file("%s/parent_structure_labelinit-labelfin_01" % ingdir)
        unsorted_02 = MASTFile("unsorted/parent_structure_labelinit-labelfin_02")
        unsorted_02.to_file("%s/parent_structure_labelinit-labelfin_02" % ingdir)
        unsorted_03 = MASTFile("unsorted/parent_structure_labelinit-labelfin_03")
        unsorted_03.to_file("%s/parent_structure_labelinit-labelfin_03" % ingdir)
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        imstrs = mywi.get_parent_image_structures()
        compare_01 = pymatgen.io.vaspio.Poscar.from_file("files/parent_structure_labelinit-labelfin_01").structure
        compare_02 = pymatgen.io.vaspio.Poscar.from_file("files/parent_structure_labelinit-labelfin_02").structure
        compare_03 = pymatgen.io.vaspio.Poscar.from_file("files/parent_structure_labelinit-labelfin_03").structure
        self.assertEqual(imstrs[0].sites, compare_01.sites)
        self.assertEqual(imstrs[0].lattice, compare_01.lattice)
        self.assertEqual(imstrs[1].sites, compare_02.sites)
        self.assertEqual(imstrs[1].lattice, compare_02.lattice)
        self.assertEqual(imstrs[2].sites, compare_03.sites)
        self.assertEqual(imstrs[2].lattice, compare_03.lattice)
        #self.testclass.get_parent_image_structures()

    def test_place_parent_energy_files(self):
        ingdir="writedir/neb_labelinit-labelfin"
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        peinit = MASTFile("files/parent_energy_labelinit")
        peinit.to_file("%s/parent_energy_labelinit" % ingdir)
        pefin = MASTFile("files/parent_energy_labelfin")
        pefin.to_file("%s/parent_energy_labelfin" % ingdir)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp_neb'
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        os.mkdir(ingdir + '/00')
        os.mkdir(ingdir + '/04')
        mywi.place_parent_energy_files()
        self.assertTrue(os.path.isfile(ingdir + '/00/OSZICAR'))
        self.assertTrue(os.path.isfile(ingdir + '/04/OSZICAR'))
        oszinit = MASTFile(ingdir + "/00/OSZICAR")
        oszfinal = MASTFile(ingdir + "/04/OSZICAR")
        self.assertEqual(oszinit.data, peinit.data)
        self.assertEqual(oszfinal.data, pefin.data)
        #self.testclass.place_parent_energy_files()

    def test_write_neb_subfolders(self):
        ingdir="writedir/neb_labelinit-labelfin"
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        parent_01 = MASTFile("statfiles/parent_structure_labelinit-labelfin_01")
        parent_01.to_file(ingdir + "/parent_structure_labelinit-labelfin_01")
        parent_02 = MASTFile("statfiles/parent_structure_labelinit-labelfin_02")
        parent_02.to_file(ingdir + "/parent_structure_labelinit-labelfin_02")
        parent_03 = MASTFile("statfiles/parent_structure_labelinit-labelfin_03")
        parent_03.to_file(ingdir + "/parent_structure_labelinit-labelfin_03")
        parent_init = MASTFile("statfiles/parent_structure_labelinit")
        parent_init.to_file(ingdir + "/parent_structure_labelinit")
        parent_fin = MASTFile("statfiles/parent_structure_labelfin")
        parent_fin.to_file(ingdir + "/parent_structure_labelfin")
        peinit = MASTFile("statfiles/parent_energy_labelinit")
        peinit.to_file("%s/parent_energy_labelinit" % ingdir)
        pefin = MASTFile("statfiles/parent_energy_labelfin")
        pefin.to_file("%s/parent_energy_labelfin" % ingdir)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp'
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images']=3
        kdict['mast_neb_settings']['lines']=neblines
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        mywi.write_neb_subfolders()
        self.assertFalse(os.path.isfile("%s/00/submit.sh" % ingdir))
        self.assertTrue(os.path.isfile("%s/01/submit.sh" % ingdir))
        self.assertTrue(os.path.isfile("%s/02/submit.sh" % ingdir))
        self.assertTrue(os.path.isfile("%s/03/submit.sh" % ingdir))
        self.assertFalse(os.path.isfile("%s/04/submit.sh" % ingdir))
        oszinit = MASTFile("%s/00/OSZICAR" % ingdir)
        oszfin = MASTFile("%s/04/OSZICAR" % ingdir)
        pos00 = MASTFile("%s/00/POSCAR" % ingdir)
        pos01 = MASTFile("%s/01/POSCAR" % ingdir)
        pos02 = MASTFile("%s/02/POSCAR" % ingdir)
        pos03 = MASTFile("%s/03/POSCAR" % ingdir)
        pos04 = MASTFile("%s/04/POSCAR" % ingdir)
        self.assertEqual(oszinit.data, peinit.data)
        self.assertEqual(oszfin.data, pefin.data)
        self.assertEqual(pos01.data, parent_01.data)
        self.assertEqual(pos02.data, parent_02.data)
        self.assertEqual(pos03.data, parent_03.data)
        self.assertEqual(pos00.data, parent_init.data)
        self.assertEqual(pos04.data, parent_fin.data)


        #self.testclass.write_neb_subfolders()

    def test_write_singlerun(self):
        ingdir="writedir/single_label1"
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp'
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        mywi.write_singlerun()
        self.assertTrue(mywi.checker.is_ready_to_run())
        #self.testclass.write_singlerun()

    def test_write_phonon_multiple(self):
        ingdir="writedir/single_phonon_label1"
        topmetad = MASTFile("files/top_metadata_single_phonon")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single_phonon")
        metad.to_file("%s/metadata.txt" % ingdir)
        shutil.copy("files/CHGCAR", ingdir)
        shutil.copy("files/WAVECAR", ingdir)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp'
        kdict['mast_phonon_settings']=dict()
        kdict['mast_phonon_settings']['phonon_center_site']="0.33 0.25 0.0"
        kdict['mast_phonon_settings']['phonon_center_radius']="1"
        kdict['mast_phonon_settings']['threshold']="0.075"
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        mywi.write_phonon_multiple()
        self.assertTrue(os.path.isdir("%s/phon_01" % ingdir))
        self.assertTrue(os.path.isdir("%s/phon_02" % ingdir))
        self.assertTrue(os.path.isdir("%s/phon_03" % ingdir))
        mywi.checker.keywords['name'] = "%s/phon_01" % ingdir
        self.assertTrue(mywi.checker.is_ready_to_run())
        mywi.checker.keywords['name'] = "%s/phon_02" % ingdir
        self.assertTrue(mywi.checker.is_ready_to_run())
        mywi.checker.keywords['name'] = "%s/phon_03" % ingdir
        self.assertTrue(mywi.checker.is_ready_to_run())
        compare_phon01 = MASTFile("files/phon_01_POSCAR")
        compare_phon02 = MASTFile("files/phon_02_POSCAR")
        compare_phon03 = MASTFile("files/phon_03_POSCAR")
        phon01 = MASTFile("%s/phon_01/POSCAR" % ingdir)
        phon02 = MASTFile("%s/phon_02/POSCAR" % ingdir)
        phon03 = MASTFile("%s/phon_03/POSCAR" % ingdir)
        self.assertEqual(phon01.data, compare_phon01.data)
        self.assertEqual(phon02.data, compare_phon02.data)
        self.assertEqual(phon03.data, compare_phon03.data)
        #self.testclass.write_phonon_multiple()

    def test_write_phonon_single(self):
        ingdir="writedir/single_phonon_label1"
        topmetad = MASTFile("files/top_metadata_single_phonon")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single_phonon")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp'
        kdict['mast_phonon_settings']=dict()
        kdict['mast_phonon_settings']['phonon_center_site']="0.33 0.25 0.0"
        kdict['mast_phonon_settings']['phonon_center_radius']="1"
        kdict['mast_phonon_settings']['threshold']="0.075"
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        mywi.write_phonon_single()
        self.assertTrue(mywi.checker.is_ready_to_run())
        compare_phon_single = MASTFile("files/phon_single_POSCAR")
        phon_singlepos = MASTFile("%s/POSCAR" % ingdir)
        self.assertEqual(phon_singlepos.data, compare_phon_single.data)
        #self.testclass.write_phonon_single()

    def test_get_my_phonon_params(self):
        ingdir="writedir/single_phonon_label1"
        topmetad = MASTFile("files/top_metadata_single_phonon")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single_phonon")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp'
        kdict['mast_phonon_settings']=dict()
        kdict['mast_phonon_settings']['phonon_center_site']="0.33 0.25 0.0"
        kdict['mast_phonon_settings']['phonon_center_radius']="1"
        kdict['mast_phonon_settings']['threshold']="0.075"
        my_structure=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure").structure
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        [mysite, myrad, mythresh] = mywi.get_my_phonon_params()
        self.assertEqual(mysite,"0.33 0.25 0.0")
        self.assertEqual(myrad,"1")
        self.assertEqual(mythresh,"0.075")
        #self.testclass.get_my_phonon_params()

    def test_write_singlerun_automesh(self):
        ingdir="writedir/single_label1"
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_kpoints']=[3,3,3,"G"]
        kdict['mast_xc']='pbe'
        kdict['mast_program']='vasp'
        kdict['mast_kpoint_density']='1000'
        mypos=pymatgen.io.vaspio.Poscar.from_file("files/perfect_structure")
        mypos.write_file("writedir/single_label1/POSCAR")
        mywi = ChopIngredient(name=ingdir,program_keys=kdict,structure=mypos.structure)
        mywi.write_singlerun_automesh()
        mykpts = pymatgen.io.vaspio.Kpoints.from_file("writedir/single_label1/KPOINTS")
        self.assertEqual(mykpts.to_dict['kpoints'][0],[8,8,8])
        #print mykpts
        self.assertTrue(mywi.checker.is_ready_to_run())
        #self.testclass.write_singlerun()
