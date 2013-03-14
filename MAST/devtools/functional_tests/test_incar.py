"""Tests for Incar"""

from MAST..incar import Incar

from nose.plugins.skip import SkipTest

class TestIncar:

    def setUp(self):
        self.testclass = Incar()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(file_path="")

    def test_create_data(self):
        raise SkipTest
        #self.testclass.create_data()

    def test_copy_to(self):
        raise SkipTest
        #self.testclass.copy_to(dest)

    def test_create_static(self):
        raise SkipTest
        #self.testclass.create_static(dest, skipencut=0)

    def test_create_static_tetrahedral(self):
        raise SkipTest
        #self.testclass.create_static_tetrahedral(dest, skipencut=0)

    def test_use_vacancy_magmom_string(self):
        raise SkipTest
        #self.testclass.use_vacancy_magmom_string(whichatom)

    def test_adjust_nelect(self):
        raise SkipTest
        #self.testclass.adjust_nelect(adjustment, num_atoms_list, zval_list)

