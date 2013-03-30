#!/usr/bin/env python
# Purpose: tests for classes
# Authors: Tam Mayeshiba
# 10/12/12 created
import os
import sys
import classes_dev
import queue_commands
import unittest
import random
import shutil
import time

testpath=os.path.expanduser("~/ttests")
cgtestname='testperov'
cgtestpath=os.path.join(testpath, cgtestname)
cgtestfiles=os.path.join(testpath, cgtestname + '_files')
cgteststem=os.path.join(testpath, cgtestname, cgtestname)

class TestSequenceFunctions(unittest.TestCase): #example tests
    def setUp(self):
        self.seq=range(10) #for examples
    def test_shuffle(self):
        # make sure the shuffled sequence does not lose any elements
        random.shuffle(self.seq)
        self.seq.sort()
        self.assertEqual(self.seq, range(10))

        # should raise an exception for an immutable sequence
        self.assertRaises(TypeError, random.shuffle, (1,2,3))

    def test_choice(self):
        element = random.choice(self.seq)
        self.assertTrue(element in self.seq)

    def test_sample(self):
        with self.assertRaises(ValueError):
            random.sample(self.seq, 20)
        for element in random.sample(self.seq, 5):
            self.assertTrue(element in self.seq)




class TestCalcgroupFunctions(unittest.TestCase):
    
    def setUp(self):
        self.tcgin = classes_dev.whole_input()
        self.tcgin.set_values_from_file(cgtestfiles + '/inputfile')
        self.tcg = classes_dev.calcgroup(self.tcgin)
        try:
            os.mkdir(cgtestpath)
        except OSError:
            shutil.rmtree(cgtestpath)
            time.sleep(2)
            os.mkdir(cgtestpath)

    def test_complete_bulk_incomplete(self):
        shutil.copytree(cgtestfiles + '/fake_bulk', cgteststem + '_bulk')
        os.mkdir(cgteststem + 'bulk_static')
        os.mkdir(cgteststem + '_ep14')
        os.mkdir(cgteststem + '_ep14_static')
        os.mkdir(cgteststem + '_ep11')
        os.mkdir(cgteststem + '_ep11_static')
        os.mkdir(cgteststem + '_neb14-11')
        os.mkdir(cgteststem + '_neb14-11_static')
        self.tcg.__init__(self.tcgin)
        self.tcg.start_calcgroup()
        self.assertEqual(self.tcg.status, 'B')
        #self.assertEqual(self.tcg.bulk_complete(), False)

    def test_complete_bulk_complete(self):
        shutil.copytree(cgtestfiles + '/fake_bulk', cgteststem + '_bulk')
        shutil.copytree(cgtestfiles + '/fake_bulk_r1', cgteststem + '_bulk_r1')
        shutil.copytree(cgtestfiles + '/fake_bulk_static', cgteststem + '_bulk_static')
        os.mkdir(cgteststem + '_ep14')
        os.mkdir(cgteststem + '_ep14_static')
        os.mkdir(cgteststem + '_ep11')
        os.mkdir(cgteststem + '_ep11_static')
        os.mkdir(cgteststem + '_neb14-11')
        os.mkdir(cgteststem + '_neb14-11_static')
        self.tcg.__init__(self.tcgin)
        self.tcg.start_calcgroup()
        self.assertEqual(self.tcg.status, 'P')
        #self.assertEqual(self.tcg.bulk_complete(), True)
    
    def test_complete_ep11_incomplete(self):
        shutil.copytree(cgtestfiles + '/fake_bulk', cgteststem + '_bulk')
        shutil.copytree(cgtestfiles + '/fake_bulk_r1', cgteststem + '_bulk_r1')
        shutil.copytree(cgtestfiles + '/fake_bulk_static', cgteststem + '_bulk_static')
        shutil.copytree(cgtestfiles + '/fake_ep11', cgteststem + '_ep11')
        os.mkdir(cgteststem + '_ep14')
        os.mkdir(cgteststem + '_ep14_static')
        os.mkdir(cgteststem + '_ep11_static')
        os.mkdir(cgteststem + '_neb14-11')
        os.mkdir(cgteststem + '_neb14-11_static')
        self.tcg.__init__(self.tcgin)
        self.tcg.start_calcgroup()
        self.assertEqual(self.tcg.status, 'P')
        #self.assertEqual(self.tcg.bulk_complete(), True)
        #self.assertEqual(self.tcg.ep_complete(), False)
    
    def test_complete_ep11_ep14_complete(self):
        shutil.copytree(cgtestfiles + '/fake_bulk', cgteststem + '_bulk')
        shutil.copytree(cgtestfiles + '/fake_bulk_r1', cgteststem + '_bulk_r1')
        shutil.copytree(cgtestfiles + '/fake_bulk_static', cgteststem + '_bulk_static')
        shutil.copytree(cgtestfiles + '/fake_ep11', cgteststem + '_ep11')
        shutil.copytree(cgtestfiles + '/fake_ep11_r1', cgteststem + '_ep11_r1')
        shutil.copytree(cgtestfiles + '/fake_ep11_static', cgteststem + '_ep11_static')
        shutil.copytree(cgtestfiles + '/fake_ep14', cgteststem + '_ep14')
        shutil.copytree(cgtestfiles + '/fake_ep14_r1', cgteststem + '_ep14_r1')
        shutil.copytree(cgtestfiles + '/fake_ep14_static', cgteststem + '_ep14_static')
        os.mkdir(cgteststem + '_neb14-11')
        os.mkdir(cgteststem + '_neb14-11_static')
        self.tcg.__init__(self.tcgin)
        self.tcg.start_calcgroup()
        self.assertEqual(self.tcg.status, 'N')
        #self.assertEqual(self.tcg.bulk_complete(), True)
        #self.assertEqual(self.tcg.ep_complete(), True)
    
    def test_complete_neb_incomplete(self):
        shutil.copytree(cgtestfiles + '/fake_bulk', cgteststem + '_bulk')
        shutil.copytree(cgtestfiles + '/fake_bulk_r1', cgteststem + '_bulk_r1')
        shutil.copytree(cgtestfiles + '/fake_bulk_static', cgteststem + '_bulk_static')
        shutil.copytree(cgtestfiles + '/fake_ep11', cgteststem + '_ep11')
        shutil.copytree(cgtestfiles + '/fake_ep11_r1', cgteststem + '_ep11_r1')
        shutil.copytree(cgtestfiles + '/fake_ep11_static', cgteststem + '_ep11_static')
        shutil.copytree(cgtestfiles + '/fake_ep14', cgteststem + '_ep14')
        shutil.copytree(cgtestfiles + '/fake_ep14_r1', cgteststem + '_ep14_r1')
        shutil.copytree(cgtestfiles + '/fake_ep14_static', cgteststem + '_ep14_static')
        shutil.copytree(cgtestfiles + '/fake_neb14-11', cgteststem + '_neb14-11')
        os.mkdir(cgteststem + '_neb14-11_static')
        os.mkdir(cgteststem + '_neb14-11_static/static00')
        os.mkdir(cgteststem + '_neb14-11_static/static01')
        os.mkdir(cgteststem + '_neb14-11_static/static02')
        self.tcg.__init__(self.tcgin)
        self.tcg.start_calcgroup()
        self.assertEqual(self.tcg.status, 'N')
        #self.assertEqual(self.tcg.bulk_complete(), True)
        #self.assertEqual(self.tcg.ep_complete(), True)
        #self.assertEqual(self.tcg.neb_complete(), False)
    
    def test_complete_neb_complete(self):
        shutil.copytree(cgtestfiles + '/fake_bulk', cgteststem + '_bulk')
        shutil.copytree(cgtestfiles + '/fake_bulk_r1', cgteststem + '_bulk_r1')
        shutil.copytree(cgtestfiles + '/fake_bulk_static', cgteststem + '_bulk_static')
        shutil.copytree(cgtestfiles + '/fake_ep11', cgteststem + '_ep11')
        shutil.copytree(cgtestfiles + '/fake_ep11_r1', cgteststem + '_ep11_r1')
        shutil.copytree(cgtestfiles + '/fake_ep11_static', cgteststem + '_ep11_static')
        shutil.copytree(cgtestfiles + '/fake_ep14', cgteststem + '_ep14')
        shutil.copytree(cgtestfiles + '/fake_ep14_r1', cgteststem + '_ep14_r1')
        shutil.copytree(cgtestfiles + '/fake_ep14_static', cgteststem + '_ep14_static')
        shutil.copytree(cgtestfiles + '/fake_neb14-11', cgteststem + '_neb14-11')
    
        shutil.copytree(cgtestfiles + '/fake_neb14-11_static', cgteststem + '_neb14-11_static')
        self.tcg.__init__(self.tcgin)
        self.tcg.start_calcgroup()
        self.assertEqual(self.tcg.status, 'C')
        #self.assertEqual(self.tcg.bulk_complete(), True)
        #self.assertEqual(self.tcg.ep_complete(), True)
        #self.assertEqual(self.tcg.neb_complete(), True)

    def test_get_vacancy_from_endpoint_path(self):
        vac=self.tcg.get_vacancy_from_endpoint_path(cgteststem + '_ep14')
        self.assertEqual(vac, 14)
        vac=self.tcg.get_vacancy_from_endpoint_path(cgteststem + '_ep4_r12')
        self.assertEqual(vac, 4)
        vac=self.tcg.get_vacancy_from_endpoint_path(cgteststem + '_ep4_static')
        self.assertEqual(vac, 4)

    def test_get_vacancies_from_neb_path(self):
        vac=self.tcg.get_vacancies_from_neb_path(cgteststem + '_neb14-11')
        self.assertEqual(vac, ['14','11'])

    def test_set_up_calc_list(self):
       
        self.tcg.__init__(self.tcgin)
        self.tcg.set_up_calc_list(1)
        self.assertEqual(self.tcg.bulk.initialpath, cgteststem + '_bulk')
        self.assertEqual(self.tcg.ep_dict, {})
        self.assertEqual(self.tcg.neb_dict, {})
        
        self.tcg.__init__(self.tcgin)
        self.tcg.set_up_calc_list(2)
        self.assertEqual(self.tcg.bulk.initialpath, cgteststem + '_bulk')
        self.assertEqual(self.tcg.ep_dict['14'].initialpath, cgteststem + '_ep14')
        self.assertEqual(self.tcg.ep_dict['11'].initialpath, cgteststem + '_ep11')
        self.assertEqual(self.tcg.neb_dict, {})

        self.tcg.__init__(self.tcgin)
        self.tcg.set_up_calc_list(3)
        self.assertEqual(self.tcg.bulk.initialpath, cgteststem + '_bulk')
        self.assertEqual(self.tcg.ep_dict['14'].initialpath, cgteststem + '_ep14')
        self.assertEqual(self.tcg.ep_dict['11'].initialpath, cgteststem + '_ep11')
        self.assertEqual(self.tcg.neb_dict['14-11'].relax_path, cgteststem + '_neb14-11')


    def test_make_vacancy_map_from_bulk(self):
        self.tcg.__init__(self.tcgin)
        map = self.tcg.make_vacancy_map_from_bulk(classes_dev.poscar(cgtestfiles + '/CONTCAR_for_bulk'), 14)
        compare = classes_dev.whole_input()
        compare.set_values_from_file(cgtestfiles + '/bulk_to_ep_map_14')
        self.assertEqual(map, compare)
    
    def test_make_vacancy_poscar_from_bulk(self):
        eposcar = self.tcg.make_vacancy_poscar_from_bulk(classes_dev.poscar(cgtestfiles + '/CONTCAR_for_bulk'), 14)
        compare = classes_dev.poscar(cgtestfiles + '/POSCAR_ep14')
        self.assertEqual(eposcar.data, compare.data)

    def test_get_endpoint_path_from_vacancy(self):
        epath = self.tcg.get_endpoint_path_from_vacancy(14)
        self.assertEqual(epath, cgteststem + '_ep14')

    def tearDown(self):
        if os.path.isfile(cgtestpath + '/output.txt'):
            shutil.copy(cgtestpath + '/output.txt',cgtestfiles + '/output' + str(time.time()))
        shutil.rmtree(cgtestpath)
        
class TestNebCalculationFunctions(unittest.TestCase):
    def setUp(self):
        self.testpath=os.path.expanduser("~/ttests")
        self.test1name='testperov'
        self.set_up_neb_calculation()
    def set_up_neb_calculation(self):
        self.tcgin = classes_dev.whole_input()
        epstem = os.path.join(self.testpath, self.test1name, self.test1name + '_ep')
        self.tcgin.set_values_from_file(os.path.join(self.testpath, self.test1name, 'inputfile'))
        self.neb=classes_dev.neb_calculation(self.tcgin, epstem + "_ep14_static",epstem + "_ep11_static")

class TestQueueFunctions(unittest.TestCase):
    def setUp(self):
        self.cgin = classes_dev.whole_input()
        self.cgin.set_values_from_file(os.path.expanduser("~/ttests/submittest/inputfile"))
        self.myrun=classes_dev.run(self.cgin)
        self.myrun.runpath=os.path.expanduser("~/ttests/submittest/")
    
    def test_submit_to_queue(self):
        myjobid=self.myrun.submit_to_queue()
        if not(myjobid == None):
            print myjobid
            self.assertTrue(type(myjobid) == int)
        else:
            self.assertTrue(False)
    
    def test_queue_commands_queue_status_from_text(self):
        jobid=93953
        qtext="93953.bardeen             checkstart       tam                    0 Q morganshort"
        status=queue_commands.queue_status_from_text(jobid, qtext)
        self.assertEqual(status, 'Q')
        qtext="93953.bardeen             checkstart       tam                    0 R morganshort"
        status=queue_commands.queue_status_from_text(jobid, qtext)
        self.assertEqual(status, 'R')
        qtext=""
        status=queue_commands.queue_status_from_text(jobid, qtext)
        self.assertEqual(status, 'X')
        qtext=None
        status=queue_commands.queue_status_from_text(jobid, qtext)
        self.assertEqual(status, 'X')


    def tearDown(self):
        pass





print ("******************************************************")
print ('Unit tests for Automatam.')
print ("******************************************************")


print ("******************************************************")
print ("Calcgroup Suite")
CalcgroupSuite = unittest.TestLoader().loadTestsFromTestCase(TestCalcgroupFunctions)
unittest.TextTestRunner(verbosity=2).run(CalcgroupSuite)
print ("******************************************************")
print ("NEB Calculation Suite")
NebCalculationSuite = unittest.TestLoader().loadTestsFromTestCase(TestNebCalculationFunctions)
unittest.TextTestRunner(verbosity=2).run(NebCalculationSuite)
print ("******************************************************")
print ("Queue Suite")
QueueSuite = unittest.TestLoader().loadTestsFromTestCase(TestQueueFunctions)
#unittest.TextTestRunner(verbosity=2).run(QueueSuite)
print ("******************************************************")
print ("Sequence Suite")
SequenceSuite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
unittest.TextTestRunner(verbosity=2).run(SequenceSuite)
sys.exit()

