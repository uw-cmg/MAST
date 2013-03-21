"""Tests for Nebcalculation"""

from MAST..nebcalculation import Nebcalculation

from nose.plugins.skip import SkipTest

class TestNebcalculation:

    def setUp(self):
        self.testclass = Nebcalculation()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(i_calc_input, epi_path="", epf_path="")

    def test_start_calc(self):
        raise SkipTest
        #self.testclass.start_calc()

    def test_rapid_update(self):
        raise SkipTest
        #self.testclass.rapid_update()

    def test_set_up(self):
        raise SkipTest
        #self.testclass.set_up(epi_path, epf_path)

    def test_reorganize_contcars(self):
        raise SkipTest
        #self.testclass.reorganize_contcars(epinitp, epfinp, nebp)

    def test_set_up_static_runs(self):
        raise SkipTest
        #self.testclass.set_up_static_runs()

    def test_update(self):
        raise SkipTest
        #self.testclass.update()

