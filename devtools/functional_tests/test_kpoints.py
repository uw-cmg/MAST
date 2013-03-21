"""Tests for Kpoints"""

from MAST..kpoints import Kpoints

from nose.plugins.skip import SkipTest

class TestKpoints:

    def setUp(self):
        self.testclass = Kpoints()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(file_path="")

    def test_create_data(self):
        raise SkipTest
        #self.testclass.create_data()

    def test_multiplied_kpoints(self):
        raise SkipTest
        #self.testclass.multiplied_kpoints()

