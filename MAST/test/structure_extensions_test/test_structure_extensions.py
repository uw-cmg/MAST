import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
import shutil
import pymatgen
import numpy as np
from MAST.utility import dirutil
from pymatgen.io.vasp import Poscar
testname ="structure_extensions_test"
testdir = dirutil.get_test_dir(testname)
#oldcontrol = os.getenv("MAST_CONTROL")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldscratch


class TestSE(unittest.TestCase):
    """Test StructureExtensions
    """
    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass
    def test_induce_defect_frac(self):
        perfect = Poscar.from_file("POSCAR_perfect").structure
        compare_vac1 = Poscar.from_file("POSCAR_vac1").structure
        compare_int1 = Poscar.from_file("POSCAR_int1").structure
        compare_sub1 = Poscar.from_file("POSCAR_sub1").structure
        coord_type='fractional'
        threshold=1.e-4
        vac1={'symbol':'O', 'type': 'vacancy', 'coordinates':  np.array([0.25, 0.75, 0.25])}
        int1={'symbol':'Ni', 'type': 'interstitial', 'coordinates': np.array([0.3, 0.3, 0.3])}
        sub1={'symbol':'Fe', 'type': 'substitution','coordinates':np.array([0.25, 0.25,0.75])}
        sxtend = StructureExtensions(struc_work1=perfect)
        struc_vac1 = sxtend.induce_defect(vac1, coord_type, threshold)
        struc_int1 = sxtend.induce_defect(int1, coord_type, threshold)
        struc_sub1 = sxtend.induce_defect(sub1, coord_type, threshold)
        self.assertEqual(struc_vac1,compare_vac1)
        self.assertEqual(struc_int1,compare_int1)
        self.assertEqual(struc_sub1,compare_sub1)
        self.assertEqual(struc_vac1.lattice,compare_vac1.lattice)
        self.assertEqual(struc_int1.lattice,compare_int1.lattice)
        self.assertEqual(struc_sub1.lattice,compare_sub1.lattice)
        self.assertEqual(struc_vac1.sites,compare_vac1.sites)
        self.assertEqual(struc_int1.sites,compare_int1.sites)
        self.assertEqual(struc_sub1.sites,compare_sub1.sites)
    
    def test_sort_structure_and_neb_lines(self):
        perfect1 = Poscar.from_file("POSCAR_defectgroup1").structure
        compare_sorted1 = Poscar.from_file("POSCAR_sorted1").structure
        perfect2 = Poscar.from_file("POSCAR_defectgroup2").structure
        compare_sorted2 = Poscar.from_file("POSCAR_sorted2").structure
        neblines = list()
        neblines.append(["Cr","0.3 0 0","0 0 0"])
        neblines.append(["Ni","0.6 0 0","0.3 0 0"])
        sxtend1 = StructureExtensions(struc_work1=perfect1, name=testdir)
        sorted1 = sxtend1.sort_structure_and_neb_lines(neblines,"00",3)
        sxtend2 = StructureExtensions(struc_work1=perfect2, name=testdir)
        sorted2 = sxtend2.sort_structure_and_neb_lines(neblines,"04",3)
        self.assertEqual(sorted1, compare_sorted1)
        self.assertEqual(sorted2, compare_sorted2)
        self.assertEqual(sorted1.lattice, compare_sorted1.lattice)
        self.assertEqual(sorted2.lattice, compare_sorted2.lattice)
        self.assertEqual(sorted1.sites, compare_sorted1.sites)
        self.assertEqual(sorted2.sites, compare_sorted2.sites)
        neblines = list()
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        sxtend1 = StructureExtensions(struc_work1=perfect1, name=testdir)
        sorted1 = sxtend1.sort_structure_and_neb_lines(neblines,"00",3)
        sxtend2 = StructureExtensions(struc_work1=perfect2, name=testdir)
        sorted2 = sxtend2.sort_structure_and_neb_lines(neblines,"04",3)
        self.assertEqual(sorted1, compare_sorted1)
        self.assertEqual(sorted2, compare_sorted2)
        self.assertEqual(sorted1.lattice, compare_sorted1.lattice)
        self.assertEqual(sorted2.lattice, compare_sorted2.lattice)
        self.assertEqual(sorted1.sites, compare_sorted1.sites)
        self.assertEqual(sorted2.sites, compare_sorted2.sites)
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        perfect3 = Poscar.from_file("POSCAR_defectgroup3").structure
        #print perfect3.get_sorted_structure()
        perfect4 = Poscar.from_file("POSCAR_defectgroup4").structure
        #print perfect4.get_sorted_structure()
        sxtend3 = StructureExtensions(struc_work1=perfect3, name=testdir)
        sorted3 = sxtend3.sort_structure_and_neb_lines(neblines,"00",3)
        #print sorted3
        sxtend4 = StructureExtensions(struc_work1=perfect4, name=testdir)
        sorted4 = sxtend4.sort_structure_and_neb_lines(neblines,"04",3)
        #print sorted4
        compare_sorted3 = Poscar.from_file("POSCAR_sorted3").structure
        #print compare_sorted3
        compare_sorted4 = Poscar.from_file("POSCAR_sorted4").structure
        #print compare_sorted4
        self.assertEqual(sorted3, compare_sorted3)
        self.assertEqual(sorted4, compare_sorted4)
        self.assertEqual(sorted3.lattice, compare_sorted3.lattice)
        self.assertEqual(sorted4.lattice, compare_sorted4.lattice)
        self.assertEqual(sorted3.sites, compare_sorted3.sites)
        self.assertEqual(sorted4.sites, compare_sorted4.sites)

    def test_sort_structure_and_neb_lines_really_scrambled(self):
        raise SkipTest
        neblines = list()
        neblines.append(["Cr","0.0 0.9 0.8","0.0 0.8 0.7"])
        neblines.append(["Cr","0.4 0.2 0.1","0.3 0.3 0.2"])
        neblines.append(["Cr","0.29 0.05 0.05","0.01 0.01 0.98"])
        neblines.append(["Ni","0.61 0.99 0.98","0.25 0.01 0.97"])
        perfect3 = Poscar.from_file("POSCAR_defectgroup3_scrambled").structure
        #print perfect3.get_sorted_structure()
        #print perfect4.get_sorted_structure()
        sxtend3 = StructureExtensions(struc_work1=perfect3)
        sorted3 = sxtend3.sort_structure_and_neb_lines(neblines,"00",3)
        compare_sorted3 = Poscar.from_file("POSCAR_sorted3").structure
        self.assertEqual(sorted3, compare_sorted3)
        self.assertEqual(sorted3.lattice, compare_sorted3.lattice)
        self.assertEqual(sorted3.sites, compare_sorted3.sites)
        perfect3 = Poscar.from_file("POSCAR_defectgroup3_really_scrambled").structure
        #print perfect3.get_sorted_structure()
        #print perfect4.get_sorted_structure()
        sxtend3 = StructureExtensions(struc_work1=perfect3)
        sorted3 = sxtend3.sort_structure_and_neb_lines(neblines,"00",3)
        compare_sorted3 = Poscar.from_file("POSCAR_sorted3").structure
        self.assertEqual(sorted3, compare_sorted3)
        self.assertEqual(sorted3.lattice, compare_sorted3.lattice)
        self.assertEqual(sorted3.sites, compare_sorted3.sites)


    def test_interpolation(self):
        ep1 = Poscar.from_file("POSCAR_ep1").structure
        ep2 = Poscar.from_file("POSCAR_ep2").structure
        compare_im1 = Poscar.from_file("POSCAR_im1").structure
        compare_im2 = Poscar.from_file("POSCAR_im2").structure
        compare_im3 = Poscar.from_file("POSCAR_im3").structure
        sxtend = StructureExtensions(struc_work1 = ep1, struc_work2 = ep2)
        slist = sxtend.do_interpolation(3)
        self.assertEqual(slist[0],ep1)
        self.assertEqual(slist[1],compare_im1)
        self.assertEqual(slist[2],compare_im2)
        self.assertEqual(slist[3],compare_im3)
        self.assertEqual(slist[4],ep2)
    def test_get_sd_array(self):
        perfect = Poscar.from_file("POSCAR_perfect").structure
        sxtend = StructureExtensions(struc_work1=perfect, name=testdir)
        mysd = sxtend.get_sd_array("0.5 0.5 0.5", 3)
        #print mysd
        myarr=np.zeros([40,3],bool)
        for idx in [8,20,25,26,28,30,32,33,36,37,38,39,40]:
            myarr[idx-1][0]=True
            myarr[idx-1][1]=True
            myarr[idx-1][2]=True
        self.assertEqual(sum(np.fabs(sum(mysd-myarr))),0)
    def test_get_sd_array_periodic_boundary(self):
        perfect = Poscar.from_file("POSCAR_perfect").structure
        sxtend = StructureExtensions(struc_work1=perfect, name=testdir)
        mysd = sxtend.get_sd_array("0.95 0.95 0.95", 1, 0.07)
        #print mysd
        myarr=np.zeros([40,3],bool)
        for idx in [1]:
            myarr[idx-1][0]=True
            myarr[idx-1][1]=True
            myarr[idx-1][2]=True
        self.assertEqual(sum(np.fabs(sum(mysd-myarr))),0)
    def test_multiple_sd_array(self):
        perfect = Poscar.from_file("POSCAR_perfect").structure
        sxtend = StructureExtensions(struc_work1=perfect, name=testdir)
        mysdlist = sxtend.get_multiple_sd_array("0.5 0.5 0.5", 1)
        mylist=list()
        mylist.append(np.zeros([40,3],bool))
        mylist.append(np.zeros([40,3],bool))
        mylist.append(np.zeros([40,3],bool))
        mylist[0][8-1][0]=True
        mylist[1][8-1][1]=True
        mylist[2][8-1][2]=True
        self.assertEqual(sum(np.fabs(sum(mysdlist[0]-mylist[0]))),0)
        self.assertEqual(sum(np.fabs(sum(mysdlist[1]-mylist[1]))),0)
        self.assertEqual(sum(np.fabs(sum(mysdlist[2]-mylist[2]))),0)
    def test_graft_coordinates(self):
        perfect = Poscar.from_file("POSCAR_perfect").structure
        coordsonly = Poscar.from_file("POSCAR_coordinates").structure
        compare_grafted = Poscar.from_file("POSCAR_grafted").structure
        sxtend = StructureExtensions(struc_work1=perfect)
        grafted = sxtend.graft_coordinates_onto_structure(coordsonly)
        self.assertEqual(grafted, compare_grafted)
        self.assertEqual(grafted.lattice, compare_grafted.lattice)
        self.assertEqual(grafted.sites, compare_grafted.sites)

    def test_strain_lattice(self):
        perfect = Poscar.from_file("POSCAR_unstrained").structure
        sxtend = StructureExtensions(struc_work1=perfect)
        strained = sxtend.strain_lattice(" 0.98 0.92 1.03  \n")
        strain_compare = Poscar.from_file("POSCAR_strained").structure
        self.assertEqual(strained, strain_compare)
        self.assertEqual(strained.lattice, strain_compare.lattice)
        self.assertEqual(strained.sites, strain_compare.sites)
    def test_scale_structure_three(self):
        hcp = Poscar.from_file("POSCAR_HCP").structure
        sxtend = StructureExtensions(struc_work1=hcp, scaling_size="2,2,2", name=testdir)
        scaled = sxtend.scale_structure()
        compare = Poscar.from_file("POSCAR_HCP_222").structure
        #self.assertEqual(scaled, Poscar.from_file("POSCAR_HCP_222").structure)
        self.assertAlmostEqual(scaled.volume, compare.volume, places=4)
        self.assertEqual(scaled.lattice, compare.lattice)
        self.assertEqual(scaled.sites.sort(), compare.sites.sort())
        return
    def test_scale_structure_nine(self):
        hcp = Poscar.from_file("POSCAR_HCP").structure
        sxtend = StructureExtensions(struc_work1=hcp, scaling_size="2 0 0,0 2 0,0 0 2", name=testdir)
        scaled = sxtend.scale_structure()
        compare = Poscar.from_file("POSCAR_HCP_222").structure
        #self.assertEqual(scaled, Poscar.from_file("POSCAR_HCP_222").structure)
        self.assertAlmostEqual(scaled.volume, compare.volume, places=4)
        self.assertEqual(scaled.lattice, compare.lattice)
        self.assertEqual(scaled.sites.sort(), compare.sites.sort())
        return
    def test_scale_defect(self):
        perfect = Poscar.from_file("POSCAR_perfect").structure
        scalingsize = "2,2,2"
        sxtend = StructureExtensions(struc_work1=perfect, scaling_size=scalingsize, name=testdir)
        scaled = sxtend.scale_structure()
        sxtend2 = StructureExtensions(struc_work1=scaled, struc_work2=perfect, scaling_size=scalingsize, name=testdir)
        vac1={'symbol':'O', 'type': 'vacancy', 'coordinates':  np.array([0.25, 0.75, 0.25])}
        defected =  sxtend2.scale_defect(vac1,'fractional',0.0001)
        int1={'symbol':'Ni', 'type': 'interstitial', 'coordinates': np.array([0.3, 0.3, 0.3])}
        sxtend3 = StructureExtensions(struc_work1=defected, struc_work2=perfect, scaling_size=scalingsize, name=testdir)
        defected2 = sxtend3.scale_defect(int1,'fractional',0.0001)
        sub1={'symbol':'Fe', 'type': 'substitution','coordinates':np.array([0.25, 0.25,0.75])}
        sxtend4 = StructureExtensions(struc_work1=defected2, struc_work2=perfect, scaling_size=scalingsize, name=testdir)
        defected3 = sxtend4.scale_defect(sub1,'fractional',0.0001)
        self.assertEqual(Poscar.from_file("POSCAR_scaled_defected").structure, defected3)
