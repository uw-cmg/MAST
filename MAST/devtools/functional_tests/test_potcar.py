"""Tests for Potcar"""

from MAST..potcar import Potcar

from nose.plugins.skip import SkipTest

class TestPotcar:

    def setUp(self):
        self.testclass = Potcar()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(file_path="")

    def test_potcar_exists(self):
        raise SkipTest
        #self.testclass.potcar_exists()

    def test_get_max_enmax(self):
        raise SkipTest
        #self.testclass.get_max_enmax()

    def test_get_zval_list(self):
        raise SkipTest
        #self.testclass.get_zval_list()

    def test_get_element_string(self):
        raise SkipTest
        #self.testclass.get_element_string()

