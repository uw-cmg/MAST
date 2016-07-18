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

testname="chop_test_update"
testdir = dirutil.get_test_dir(testname)
oldscratch = os.getenv("MAST_SCRATCH")

class TestUpdateChildrenIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir("writedir"):
            os.mkdir("writedir")
        if not os.path.isdir("writedir/next_ingred"):
            os.mkdir("writedir/next_ingred")
        if not os.path.isdir("writedir/single_label1"):
            os.mkdir("writedir/single_label1")
        if not os.path.isdir("writedir/neb_labelinit-labelfin"):
            os.mkdir("writedir/neb_labelinit-labelfin")
        if not os.path.isdir("writedir/single_phonon_label1"):
            os.mkdir("writedir/single_phonon_label1")
        os.environ['MAST_SCRATCH'] = testdir

    def tearDown(self):
        tearlist = list()
        tearlist.append("writedir")
        os.environ['MAST_SCRATCH'] = oldscratch
        #tearlist.append("writedir/single_label1")
        #tearlist.append("writedir/next_ingred")
        #tearlist.append("writedir/neb_labelinit-labelfin")
        #tearlist.append("writedir/single_phonon_label1")
        for foldername in tearlist:
            try:
                shutil.rmtree(foldername)
            except OSError:
                pass

    def test___init__(self):
        self.assertTrue(True)
        #self.testclass.__init__(**kwargs)

    def test__fullpath_childname(self):
        ingdir="%s/writedir/single_label1" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myuci = ChopIngredient(name=ingdir,program_keys=kdict,structure=my_structure)
        fullpath = myuci._fullpath_childname("next_ingred")
        self.assertEqual(fullpath, "%s/writedir/next_ingred" % testdir)
        #self.testclass._fullpath_childname(childname)

    def test_give_structure(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        #self.testclass.give_structure(childname)

    def test_give_neb_structures_to_neb(self):
        ingdir="%s/writedir/neb_labelinit-labelfin" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp_neb'
        kdict['mast_neb_settings']=dict()
        kdict['mast_neb_settings']['images'] = 3
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed=dict()
        for subdir in ['00','01','02','03','04']:
            os.mkdir("writedir/neb_labelinit-labelfin/%s" % subdir)
            myrelaxed[subdir] = MASTFile("files/POSCAR_%s" % subdir)
            myrelaxed[subdir].to_file("writedir/neb_labelinit-labelfin/%s/CONTCAR" % subdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_neb_structures_to_neb("next_ingred")
        for subdir in ['01','02','03']:
            givenstr = MASTFile("%s/writedir/next_ingred/parent_structure_labelinit-labelfin_%s" % (testdir,subdir))
            self.assertEqual(myrelaxed[subdir].data, givenstr.data)
        #self.testclass.give_neb_structures_to_neb(childname)

    def test_give_saddle_structure(self):
        ingdir="%s/writedir/neb_labelinit-labelfin" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_neb")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_neb")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp_neb'
        kdict['images'] = 3
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed=dict()
        myosz=dict()
        mywav=dict()
        mychg=dict()
        for subdir in ['00','01','02','03','04']:
            os.mkdir("writedir/neb_labelinit-labelfin/%s" % subdir)
            myrelaxed[subdir] = MASTFile("files/POSCAR_%s" % subdir)
            myrelaxed[subdir].to_file("writedir/neb_labelinit-labelfin/%s/CONTCAR" % subdir)
            myosz[subdir] = MASTFile("files/OSZICAR_%s" % subdir)
            myosz[subdir].to_file("writedir/neb_labelinit-labelfin/%s/OSZICAR" % subdir)
            mychg[subdir] = MASTFile("files/CHGCAR")
            mychg[subdir].to_file("writedir/neb_labelinit-labelfin/%s/CHGCAR" % subdir)
            mywav[subdir] = MASTFile("files/WAVECAR")
            mywav[subdir].to_file("writedir/neb_labelinit-labelfin/%s/WAVECAR" % subdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_saddle_structure("next_ingred") #should be OSZ3
        saddle = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed['03'].data, saddle.data)
        saddledir = myuci.get_saddle_dir()
        self.assertEqual(saddledir, "03")
        #self.testclass.give_saddle_structure(childname)

    def test_give_phonon_multiple_forces_and_displacements(self):
        ingdir="%s/writedir/single_phonon_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single_phonon")
        metad.to_file("%s/metadata.txt" % ingdir)
        mypos = MASTFile("files/phonon_initial_POSCAR")
        mypos.to_file("%s/POSCAR" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myxdat=dict()
        mydynmat=dict()
        for subdir in ['phon_01','phon_02','phon_03']:
            os.mkdir("%s/%s" % (ingdir,subdir))
            myxdat[subdir] = MASTFile("files/XDATCAR_%s" % subdir)
            myxdat[subdir].to_file("%s/%s/XDATCAR" % (ingdir,subdir))
            mydynmat[subdir] = MASTFile("files/DYNMAT_%s" % subdir)
            mydynmat[subdir].to_file("%s/%s/DYNMAT" % (ingdir, subdir))
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_phonon_multiple_forces_and_displacements("next_ingred") 
        newpos = MASTFile("%s/writedir/next_ingred/POSCAR_prePHON" % testdir)
        newdyn = MASTFile("%s/writedir/next_ingred/DYNMAT" % testdir)
        newxdat = MASTFile("%s/writedir/next_ingred/XDATCAR" % testdir)
        comparepos = MASTFile("files/phonon_initial_POSCAR")
        comparedyn = MASTFile("files/DYNMAT_compare")
        comparexdat = MASTFile("%s/XDATCAR" % ingdir)
        self.assertEqual(newpos.data, comparepos.data)
        self.assertEqual(newdyn.data, comparedyn.data)
        self.assertEqual(newxdat.data, comparexdat.data)
        #self.testclass.give_phonon_multiple_forces_and_displacements(childname)

    def test_give_phonon_single_forces_and_displacements(self):
        ingdir="%s/writedir/single_phonon_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single_phonon")
        metad.to_file("%s/metadata.txt" % ingdir)
        mypos = MASTFile("files/phonon_initial_POSCAR")
        mypos.to_file("%s/POSCAR" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myxdat = MASTFile("files/XDATCAR_compare")
        myxdat.to_file("%s/XDATCAR" % ingdir)
        mydynmat = MASTFile("files/DYNMAT_compare")
        mydynmat.to_file("%s/DYNMAT" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_phonon_single_forces_and_displacements("next_ingred") 
        newpos = MASTFile("%s/writedir/next_ingred/POSCAR_prePHON" % testdir)
        newdyn = MASTFile("%s/writedir/next_ingred/DYNMAT" % testdir)
        newxdat = MASTFile("%s/writedir/next_ingred/XDATCAR" % testdir)
        comparepos = MASTFile("files/phonon_initial_POSCAR")
        comparedyn = MASTFile("files/DYNMAT_compare")
        comparexdat = MASTFile("%s/XDATCAR" % ingdir)
        self.assertEqual(newpos.data, comparepos.data)
        self.assertEqual(newdyn.data, comparedyn.data)
        self.assertEqual(newxdat.data, comparexdat.data)
        #self.testclass.give_phonon_single_forces_and_displacements(childname)

    def test_give_structure_and_energy_to_neb(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.data.append("defect_label = label1\n")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        myenergy = MASTFile("files/OSZICAR_relaxed")
        myenergy.to_file("%s/OSZICAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_energy_to_neb("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/parent_structure_label1" % testdir)
        givenenergy = MASTFile("%s/writedir/next_ingred/parent_energy_label1" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertEqual(myenergy.data, givenenergy.data)
        #self.testclass.give_structure_and_energy_to_neb(childname)

    def test_give_structure_and_restart_files(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        mychg = MASTFile("files/CHGCAR")
        mychg.to_file("%s/CHGCAR" % ingdir)
        mywave = MASTFile("files/WAVECAR")
        mywave.to_file("%s/WAVECAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_restart_files("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/WAVECAR" % testdir))
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/CHGCAR" % testdir))
        #self.testclass.give_structure_and_restart_files(childname)

    def test_give_structure_and_restart_files_full_copies(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        mychg = MASTFile("files/CHGCAR")
        mychg.to_file("%s/CHGCAR" % ingdir)
        mywave = MASTFile("files/WAVECAR")
        mywave.to_file("%s/WAVECAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_restart_files_full_copies("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/WAVECAR" % testdir))
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/CHGCAR" % testdir))
        #self.testclass.give_structure_and_restart_files(childname)

    def test_give_structure_and_restart_files_softlinks(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        mychg = MASTFile("files/CHGCAR")
        mychg.to_file("%s/CHGCAR" % ingdir)
        mywave = MASTFile("files/WAVECAR")
        mywave.to_file("%s/WAVECAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_restart_files_softlinks("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/WAVECAR" % testdir))
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/CHGCAR" % testdir))
        #self.testclass.give_structure_and_restart_files(childname)

    def test_give_structure_and_charge_density_softlink(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        mychg = MASTFile("files/CHGCAR")
        mychg.to_file("%s/CHGCAR" % ingdir)
        mywave = MASTFile("files/WAVECAR")
        mywave.to_file("%s/WAVECAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_charge_density_softlink("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/CHGCAR" % testdir))
        #self.testclass.give_structure_and_restart_files(childname)

    def test_give_structure_and_wavefunction_softlink(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        mychg = MASTFile("files/CHGCAR")
        mychg.to_file("%s/CHGCAR" % ingdir)
        mywave = MASTFile("files/WAVECAR")
        mywave.to_file("%s/WAVECAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_wavefunction_softlink("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/WAVECAR" % testdir))
        #self.testclass.give_structure_and_restart_files(childname)
    def test_give_structure_and_wavefunction_full_copy(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        mychg = MASTFile("files/CHGCAR")
        mychg.to_file("%s/CHGCAR" % ingdir)
        mywave = MASTFile("files/WAVECAR")
        mywave.to_file("%s/WAVECAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_wavefunction_full_copy("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/WAVECAR" % testdir))
        #self.testclass.give_structure_and_restart_files(childname)
    def test_give_structure_and_charge_density_full_copy(self):
        ingdir="%s/writedir/single_label1" % testdir
        recipedir="%s/writedir" % testdir
        topmetad = MASTFile("files/top_metadata_single")
        topmetad.data.append("origin_dir = %s/files\n" % testdir) #give origin directory
        topmetad.to_file("writedir/metadata.txt")
        metad = MASTFile("files/metadata_single")
        metad.to_file("%s/metadata.txt" % ingdir)
        kdict=dict()
        kdict['mast_program'] = 'vasp'
        my_structure = Poscar.from_file("files/perfect_structure").structure
        myrelaxed = MASTFile("files/relaxed_structure")
        myrelaxed.to_file("%s/CONTCAR" % ingdir)
        mychg = MASTFile("files/CHGCAR")
        mychg.to_file("%s/CHGCAR" % ingdir)
        mywave = MASTFile("files/WAVECAR")
        mywave.to_file("%s/WAVECAR" % ingdir)
        myuci = ChopIngredient(name=ingdir,program_keys=kdict, structure=my_structure)
        myuci.give_structure_and_charge_density_full_copy("next_ingred")
        givenstr = MASTFile("%s/writedir/next_ingred/POSCAR" % testdir)
        self.assertEqual(myrelaxed.data, givenstr.data)
        self.assertTrue(os.path.exists("%s/writedir/next_ingred/CHGCAR" % testdir))
        #self.testclass.give_structure_and_restart_files(childname)
