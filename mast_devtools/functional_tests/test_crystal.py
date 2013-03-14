"""Tests for Crystal"""

from MAST..crystal import Crystal

from nose.plugins.skip import SkipTest

class TestCrystal:

    def setUp(self):
        self.testclass = Crystal()

    def tearDown(self):
        pass

    def test_make_fcc(self):
        raise SkipTest
        #self.testclass.make_fcc(a,x,y,z,symbol)

    def test_make_hcp(self):
        raise SkipTest
        #self.testclass.make_hcp(a,x,y,z,symbol)

    def test_make_poscar_from_file(self):
        raise SkipTest
        #self.testclass.make_poscar_from_file(x,y,z,symbol,filepath)

    def test_resize_poscar(self):
        raise SkipTest
        #self.testclass.resize_poscar(x,y,z)

