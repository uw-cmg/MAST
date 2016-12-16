"""Tests for Atom Indexing"""

from MAST.parsers.inputparser import InputParser
from MAST.recipe.recipesetup import RecipeSetup
from MAST.ingredients.pmgextend.atom_index import AtomIndex
from MAST.ingredients.chopingredient import ChopIngredient
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
from MAST.utility import MASTError
from MAST.controllers.mastmon import MASTMon
import unittest
import numpy as np
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import InputOptions
from MAST.controllers.mastinput import MASTInput
from MAST.utility import MASTFile
from MAST.utility import Metadata
from pymatgen.io.vasp import Poscar
import shutil
testname="atom_index_test"
testdir = dirutil.get_test_dir(testname)

class TestAtomIndexing(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir(os.path.join(testdir,'workdir')):
            os.mkdir('workdir')

    def tearDown(self):
        #pass
        #return
        for dirname in ['workdir','structure_index_files']:
            dfull = os.path.join(testdir, dirname)
            if os.path.isdir(dfull):
                shutil.rmtree(dfull)

    def test_interstitial_neb_setup(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        tdir=os.path.join(testdir,'interstitial_neb_files')
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        #print rwdir
        fnames = os.listdir(os.path.join(rwdir, "structure_index_files"))
        for fname in fnames:
            f1=MASTFile(os.path.join(tdir,"structure_index_files",fname))
            f2=MASTFile(os.path.join(rwdir,"structure_index_files",fname))
            self.assertEqual(f1.data,f2.data)
        #iopt=InputParser(inputfile='neb_pathfinder.inp').parse()
        #rfile = iopt.get_item('personal_recipe', 'personal_recipe_list')
        #print rfile
        #struc = iopt.get_item('structure','structure')
        #myrs=RecipeSetup(recipeFile=rfile, workingDirectory=wdir, inputOptions=iopt, structure=struc)
        #myrs.start()
        #self.assertEqual(type(myrs.metafile),MAST.utility.metadata.Metadata)
        return

    def test__init__(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=mysid)
        self.assertEqual(myai.atomcount, 1)
        return
    
    def test__init__2(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='multidefect.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=mysid)
        self.assertEqual(myai.atomcount, 1)
        return
    
    def test_write_defected(self):
        #raise SkipTest
        tdir=os.path.join(testdir,'multidefect_files')
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='multidefect.inp')
        os.environ['MAST_SCRATCH']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        fnames = os.listdir(os.path.join(rwdir, "structure_index_files"))
        for fname in fnames:
            f1=MASTFile(os.path.join(tdir,"structure_index_files",fname))
            f2=MASTFile(os.path.join(rwdir,"structure_index_files",fname))
            self.assertEqual(f1.data,f2.data)
        return

    def test_find_orig_frac_coord_in_atom_indices(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        os.environ['MAST_SCRATCH']=wdir
        myipparser=InputParser(inputfile=os.path.join(testdir,"neb_pathfinder.inp"))
        myio = myipparser.parse()
        test_sid = os.path.join(testdir, "find_coord_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        orig_coord=np.array([0.0,0.0,0.5],'float')
        self.assertRaises(MASTError, myai.find_frac_coord_in_atom_indices,
                orig_coord,include_orig="only",element="")
        print "subtest1 ok"
        findtest2 = myai.find_frac_coord_in_atom_indices(orig_coord,
                include_orig="only",element="",
                scaling_label="",find_multiple=True)
        self.assertItemsEqual(findtest2, ["0000000000000x12","0000000000000xE2"])
        print "subtest2 ok"
        findtest3 = myai.find_frac_coord_in_atom_indices(orig_coord,
                include_orig="only",element="He")
        self.assertEqual(findtest3, "0000000000000xE2")
        print "subtest3 ok"
        findtest4 = myai.find_frac_coord_in_atom_indices(orig_coord,
                include_orig="only",element="Al",
            find_multiple=True, tol=0.01)
        self.assertItemsEqual(findtest4, ["0000000000000x12","0000000000000xTOL"])
        print "subtest4 ok"
        return
    
    def test_find_frac_coord_in_atom_indices(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        os.environ['MAST_SCRATCH']=wdir
        myipparser=InputParser(inputfile=os.path.join(testdir,"neb_pathfinder.inp"))
        myio = myipparser.parse()
        test_sid = os.path.join(testdir, "find_coord_files")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        relaxed_coord=np.array([-0.00501174,-0.00501174, 0.50707951],'float')
        self.assertRaises(MASTError, myai.find_frac_coord_in_atom_indices,
                relaxed_coord, include_orig="no", element="")
        print "subtest1 ok"
        findtest2 = myai.find_frac_coord_in_atom_indices(relaxed_coord,
                include_orig="no", element="",
                scaling_label="",find_multiple=True)
        self.assertItemsEqual(findtest2, ["0000000000000x12","0000000000000xE2",
                "0000000000000xTOL"])
        print "subtest2 ok"
        findtest3 = myai.find_frac_coord_in_atom_indices(relaxed_coord,
                include_orig="no", element="He")
        self.assertEqual(findtest3, "0000000000000xE2")
        print "subtest3 ok"
        findtest4 = myai.find_frac_coord_in_atom_indices(relaxed_coord,
                include_orig="no",element="Al",
                find_multiple=True, tol=0.001)
        self.assertItemsEqual(findtest4, ["0000000000000x12","0000000000000xTOL"])
        print "subtest4 ok"
        return
    def test_write_defected_phonon_sd_manifests(self):
        raise SkipTest
        return

    def test_write_neb_endpoint_manifests(self):
        raise SkipTest
        return

    def test_write_neb_phonon_sd_manifests(self):
        raise SkipTest
        return

    def test_set_up_initial_index(self):
        raise SkipTest
        return

    def test_update_atom_indices_from_structure(self):
        raise SkipTest
        return

    def test_make_coordinate_and_element_list_from_manifest(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        myip = MASTInput(inputfile='multidefect.inp')
        os.environ['MAST_SCRATCH']=wdir
        os.environ['MAST_CONTROL']=wdir
        myip.set_up_recipe()
        rwdir = myip.working_directory
        myio = myip.input_options
        mysid = os.path.join(rwdir, "structure_index_files")
        test_sid = os.path.join(testdir,"list_files","structure_index_files_orig_only")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        [coordlist, elemlist]=myai.make_coordinate_and_element_list_from_manifest("manifest__group1_")
        mystr = myai.startstr
        newstr = myai.graft_new_coordinates_from_manifest(mystr, "manifest__group1_","")
        comparestr = Poscar.from_file(os.path.join(testdir,"list_files","POSCAR_test1")).structure
        self.assertEqual(comparestr, newstr)
        print "subtest1 ok"
        test_sid = os.path.join(testdir,"list_files","structure_index_files_updated")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        [coordlist, elemlist]=myai.make_coordinate_and_element_list_from_manifest("manifest__group1_", "test2")
        mystr = myai.startstr
        newstr = myai.graft_new_coordinates_from_manifest(mystr, "manifest__group1_","test2")
        comparestr = Poscar.from_file(os.path.join(testdir,"list_files","POSCAR_test2")).structure
        self.assertEqual(comparestr, newstr)
        print "subtest2 ok"
        test_sid = os.path.join(testdir,"list_files","structure_index_files_updated")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        [coordlist, elemlist]=myai.make_coordinate_and_element_list_from_manifest("manifest__group1_", "test3")
        mystr = myai.startstr
        newstr = myai.graft_new_coordinates_from_manifest(mystr, "manifest__group1_","test3")
        comparestr = Poscar.from_file(os.path.join(testdir,"list_files","POSCAR_test3")).structure
        self.assertEqual(comparestr, newstr)
        print "subtest3 ok"
        test_sid = os.path.join(testdir,"list_files","structure_index_files_updated")
        myai = AtomIndex(input_options=myio, structure_index_directory=test_sid)
        [coordlist, elemlist]=myai.make_coordinate_and_element_list_from_manifest("manifest__group1_", "test4")
        mystr = myai.startstr
        newstr = myai.graft_new_coordinates_from_manifest(mystr, "manifest__group1_","test4")
        comparestr = Poscar.from_file(os.path.join(testdir,"list_files","POSCAR_test4")).structure
        self.assertEqual(comparestr, newstr)
        print "subtest4 ok"
        return

    def test_guess_manifest_from_ingredient_metadata(self):
        raise SkipTest
        return

    def test_make_temp_manifest_from_scrambled_structure(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        os.environ['MAST_SCRATCH']=wdir
        os.environ['MAST_ARCHIVE']=wdir
        os.environ['MAST_CONTROL']=wdir
        myip = MASTInput(inputfile='neb_pathfinder.inp')
        myip.set_up_recipe()
        skipfile = MASTFile()
        skipfile.data.append("Skip")
        rwdir = myip.working_directory
        skipfile.to_file(os.path.join(rwdir,"MAST_SKIP"))
        dirutil.lock_directory(wdir)
        mymon = MASTMon()
        #dirutil.unlock_directory(wdir)
        #The MAST_SKIP file inside the folder should prevent any actions
        #from the previous step
        myrp = mymon.set_up_recipe_plan(rwdir, verbose=1)
        os.chdir(rwdir)
        ingname="neb_int1-int2_q=p0_opt1"
        mysid = os.path.join(rwdir, "structure_index_files")
        test_sid = os.path.join(testdir, "nebpathtest_condensed","structure_index_files")
        fnames = os.listdir(test_sid)
        ingdir = os.path.join(rwdir, ingname)
        for fname in fnames: #seed structure index with relaxed files
            shutil.copy(os.path.join(test_sid,fname), mysid) 
        shutil.copy(os.path.join(testdir,"nebpathtest_files","parent_structure_int1"),ingdir)
        shutil.copy(os.path.join(testdir,"nebpathtest_files","parent_structure_int2"),ingdir)
        #try direct match; works
        myai = AtomIndex(input_options=myip.input_options, structure_index_directory=mysid)
        mystr = Poscar.from_file(os.path.join(testdir,"nebpathtest_files",
                    "parent_structure_int2")).structure
        myai.make_temp_manifest_from_scrambled_structure(ingdir,
            mystr, os.path.join(wdir, "manifest_int2_test1"))
        man1=MASTFile(os.path.join(wdir,"manifest_int2_test1"))
        mancomp=MASTFile(os.path.join(testdir,"nebpathtest_files","manifest_int2_test1"))
        self.assertItemsEqual(man1.data, mancomp.data)
        print "subtest1 ok"
        ###chop
        my_ing = ChopIngredient(name = myrp.ingred_input_options[ingname]['name'], 
            program_keys = myrp.ingred_input_options[ingname]['program_keys'], 
            structure = myrp.ingred_input_options[ingname]['structure'])
        parentstructures = my_ing.get_parent_structures()
        sxtend = StructureExtensions(struc_work1=parentstructures[0], struc_work2=parentstructures[1],name=my_ing.keywords['name'])
        pmg_image_structures = sxtend.do_interpolation(my_ing.keywords['program_keys']['mast_neb_settings']['images'])
        #print "IMAGE STRUCTURE -1:"
        #print pmg_image_structures[-1]
        my_ing.checker.write_poscar_with_zero_velocities(Poscar(pmg_image_structures[-1]),os.path.join(my_ing.keywords['name'],"POSCAR_final_image_structure"))
        myai.make_temp_manifest_from_scrambled_structure(my_ing.keywords['name'],pmg_image_structures[-1],os.path.join(my_ing.keywords['name'],'pmg_scrambledep2'))
        my_ing.checker.write_poscar_with_zero_velocities(Poscar(pmg_image_structures[0]),os.path.join(my_ing.keywords['name'],"POSCAR_initial_image_structure"))
        myai.make_temp_manifest_from_scrambled_structure(my_ing.keywords['name'],pmg_image_structures[0],os.path.join(my_ing.keywords['name'],'pmg_scrambledep1'))
        man1_test = MASTFile(os.path.join(my_ing.keywords['name'],'pmg_scrambledep1'))
        man1_comp = MASTFile(os.path.join(testdir,"nebpathtest_files","pmg_scrambledep1"))
        self.assertEqual(man1_test.data, man1_comp.data)
        print "subtest 2.1 ok"
        man2_test = MASTFile(os.path.join(my_ing.keywords['name'],'pmg_scrambledep2'))
        man2_comp = MASTFile(os.path.join(testdir,"nebpathtest_files","pmg_scrambledep2"))
        self.assertEqual(man2_test.data, man2_comp.data)
        print "subtest 2.2 ok"
        return

    def test_unscramble_a_scrambled_structure(self):
        raise SkipTest
        return

    def test_graft_new_coordinates_from_manifest(self):
        raise SkipTest
        return

    def test_get_sd_array(self):
        raise SkipTest
        return

    def test_add_atom_specific_keywords_to_structure_dictionary(self):
        raise SkipTest
        return

    def test_add_element_specific_keywords_to_structure_dictionary(self):
        raise SkipTest
        return

    def test_defect_static(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        tdir=os.path.join(testdir,'nebpathtest_condensed')
        rwdir=os.path.join(wdir,"nebpathtest_condensed")
        shutil.copytree(tdir,rwdir)
        os.environ['MAST_SCRATCH']=wdir
        os.environ['MAST_ARCHIVE']=wdir
        os.environ['MAST_CONTROL']=wdir
        dirutil.lock_directory(wdir)
        mymon = MASTMon()
        #The MAST_SKIP file inside the folder should prevent any actions
        #from the previous step
        myrp = mymon.set_up_recipe_plan(rwdir, verbose=1)
        os.chdir(rwdir)
        #myrp.check_recipe_status(1)
        ingp = "defect_int2_q=p0_opt2"
        ingc = "defect_int2_q=p0_stat"
        iscomplete=myrp.complete_ingredient(ingp)
        self.assertTrue(iscomplete)
        myrp.ingredients[ingp]="C"
        #myrp.update_children(ingp)
        #myrp.do_ingredient_methods(ingp,"mast_update_children_method",ingc)
        #print myrp.update_methods[ingp][ingc]
        #mresult = myrp.run_a_method(ingp, "give_structure", [ingc])
        my_ing = ChopIngredient(name = myrp.ingred_input_options[ingp]['name'], 
            program_keys = myrp.ingred_input_options[ingp]['program_keys'], 
            structure = myrp.ingred_input_options[ingp]['structure'])
        #my_ing.give_structure(ingc)
        childname = my_ing._fullpath_childname(ingc)
        #my_ing.checker.forward_final_structure_file(childname)
        checker = my_ing.checker
        ###from vaspchecker
        childpath=childname
        newname="POSCAR"
        workdir=os.path.dirname(checker.keywords['name'])
        sdir = os.path.join(workdir, "structure_index_files")
        childmeta=Metadata(metafile="%s/metadata.txt" % childpath)
        child_program=childmeta.read_data("program")
        if not "vasp" in child_program: #madelung utility or another folder
            return checker.copy_a_file(childpath, "CONTCAR", newname)
        child_scaling_label=childmeta.read_data("scaling_label")
        child_defect_label=childmeta.read_data("defect_label")
        child_neb_label=childmeta.read_data("neb_label")
        child_phonon_label=childmeta.read_data("phonon_label")
        if child_scaling_label == None:
            child_scaling_label = ""
        if child_defect_label == None:
            child_defect_label = ""
        if child_neb_label == None:
            child_neb_label = ""
        if child_phonon_label == None:
            child_phonon_label = ""
        parentmeta=Metadata(metafile="%s/metadata.txt" % checker.keywords['name'])
        parent_defect_label=parentmeta.read_data("defect_label")
        parent_neb_label=parentmeta.read_data("neb_label")
        if parent_defect_label == None:
            parent_defect_label = ""
        if parent_neb_label == None:
            parent_neb_label = ""
        if (not (child_neb_label == "")) and (not (parent_defect_label == "")):
            child_defect_label = parent_defect_label
        if (not (child_phonon_label == "")):
            if (not (parent_defect_label == "")):
                child_defect_label = parent_defect_label
            if (not (parent_neb_label == "")):
                child_neb_label = parent_neb_label
                child_defect_label = parent_neb_label.split('-')[0].strip() # always define phonons from first endpoint
        #get child manifest
        childmanifest="manifest_%s_%s_%s" % (child_scaling_label, child_defect_label, child_neb_label)
        #build structure from atom indices using parent name_frac_coords
        ing_label=os.path.basename(checker.keywords['name'])
        childmeta.write_data("parent",ing_label)
        mystr=Poscar.from_file("%s/CONTCAR" % checker.keywords['name']).structure
        myatomindex=AtomIndex(structure_index_directory=sdir)
        if "inducescaling" in childpath: #initial scaled coords have no parent
            ing_label="original"
        self.assertEqual(childmanifest,"manifest__int2_")
        self.assertEqual(ing_label,ingp)
        [coordlist, elemlist]=myatomindex.make_coordinate_and_element_list_from_manifest(childmanifest, ing_label)
        for cidx in range(0, len(coordlist)):
            print "%s, %s" % (elemlist[cidx], coordlist[cidx])
        #newstr = mystr.copy()
        #lenoldsites = len(newstr.sites)
        #newstr.remove_sites(range(0, lenoldsites))
        #for cct in range(0, len(coordlist)):
        #    newstr.append(elemlist[cct], 
        #        coordlist[cct],
        #        coords_are_cartesian=False,validate_proximity=True) 
        newstr=myatomindex.graft_new_coordinates_from_manifest(mystr, childmanifest,ing_label)
        newposcar=Poscar(newstr)
        checker.write_poscar_with_zero_velocities(newposcar, os.path.join(childpath, newname))
        pos_test = MASTFile(os.path.join(childpath, "POSCAR"))
        pos_comp = MASTFile(os.path.join(testdir, "nebpathtest_files",
            "POSCAR_defect_int2_stat"))
        self.assertEqual(pos_test.data, pos_comp.data)
        ########
        #from MAST monitor
        #myipparser=InputParser(inputfile=os.path.join(rwdir,"input.inp"))
        #myinputoptions = myipparser.parse()
        #personal_recipe_contents = myinputoptions.get_item('personal_recipe',
        #        'personal_recipe_list')
        #rsetup = RecipeSetup(recipeFile=personal_recipe_contents,
        #        inputOptions=myinputoptions,
        #        structure=myinputoptions.get_item('structure','structure'),
        #        workingDirectory=rwdir)
        #recipe_plan_obj = rsetup.start()
        #recipe_plan_obj.get_statuses_from_file()
        #os.chdir(rwdir)
        #recipe_plan_obj.check_recipe_status(1,"defect_int2_q=p0_opt2")
        return
    
    def test_nebpathtest(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        tdir=os.path.join(testdir,'nebpathtest_condensed')
        rwdir=os.path.join(wdir,"nebpathtest_condensed")
        shutil.copytree(tdir,rwdir)
        os.environ['MAST_SCRATCH']=wdir
        os.environ['MAST_ARCHIVE']=wdir
        os.environ['MAST_CONTROL']=wdir
        dirutil.lock_directory(wdir)
        mymon = MASTMon()
        #The MAST_SKIP file inside the folder should prevent any actions
        #from the previous step
        myrp = mymon.set_up_recipe_plan(rwdir, verbose=1)
        os.chdir(rwdir)
        myrp.check_recipe_status()
        #check that statics set up correctly
        str1_1 = Poscar.from_file(os.path.join(rwdir,"defect_int1_q=p0_stat","POSCAR")).structure
        str1_comp = Poscar.from_file(os.path.join(testdir,"nebpathtest_files","POSCAR_defect_int1_stat")).structure
        self.assertEqual(str1_1,str1_comp)
        str2_1 = Poscar.from_file(os.path.join(rwdir,"defect_int2_q=p0_stat","POSCAR")).structure
        str2_comp = Poscar.from_file(os.path.join(testdir,"nebpathtest_files","POSCAR_defect_int2_stat")).structure
        self.assertEqual(str2_1,str2_comp)
        #mimic completion of statics
        for dnum in [1,2]:
            shutil.copy(os.path.join(rwdir,"defect_int%i_q=p0_opt2" % dnum,
                "OUTCAR"),os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum))
            shutil.copy(os.path.join(rwdir,"defect_int%i_q=p0_opt2" % dnum,
                "OSZICAR"),os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum))
            shutil.copy(os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum,
                "POSCAR"),os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum,
                    "CONTCAR"))
        #check again; try to set up NEB
        myrp = mymon.set_up_recipe_plan(rwdir, verbose=1)
        os.chdir(rwdir)
        myrp.check_recipe_status()
        return
    
    def test_nebtest(self):
        #raise SkipTest
        wdir=os.path.join(testdir,'workdir')
        tdir=os.path.join(testdir,'nebtest_condensed')
        rwdir=os.path.join(wdir,"nebtest_condensed")
        shutil.copytree(tdir,rwdir)
        os.environ['MAST_SCRATCH']=wdir
        os.environ['MAST_ARCHIVE']=wdir
        os.environ['MAST_CONTROL']=wdir
        dirutil.lock_directory(wdir)
        mymon = MASTMon()
        #The MAST_SKIP file inside the folder should prevent any actions
        #from the previous step
        myrp = mymon.set_up_recipe_plan(rwdir, verbose=1)
        os.chdir(rwdir)
        myrp.check_recipe_status()
        #check that statics set up correctly
        str1_1 = Poscar.from_file(os.path.join(rwdir,"defect_int1_q=p0_stat","POSCAR")).structure
        str1_comp = Poscar.from_file(os.path.join(testdir,"nebpathtest_files","POSCAR_defect_int1_stat")).structure
        self.assertEqual(str1_1,str1_comp)
        str2_1 = Poscar.from_file(os.path.join(rwdir,"defect_int2_q=p0_stat","POSCAR")).structure
        str2_comp = Poscar.from_file(os.path.join(testdir,"nebpathtest_files","POSCAR_defect_int2_stat")).structure
        self.assertEqual(str2_1,str2_comp)
        #mimic completion of statics
        for dnum in [1,2]:
            shutil.copy(os.path.join(rwdir,"defect_int%i_q=p0_opt2" % dnum,
                "OUTCAR"),os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum))
            shutil.copy(os.path.join(rwdir,"defect_int%i_q=p0_opt2" % dnum,
                "OSZICAR"),os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum))
            shutil.copy(os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum,
                "POSCAR"),os.path.join(rwdir,"defect_int%i_q=p0_stat" % dnum,
                    "CONTCAR"))
        #check again; try to set up NEB
        myrp = mymon.set_up_recipe_plan(rwdir, verbose=1)
        os.chdir(rwdir)
        myrp.check_recipe_status()
        return
