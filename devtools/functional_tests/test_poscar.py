"""Tests for Poscar"""

from MAST..poscar import Poscar

from nose.plugins.skip import SkipTest

class TestPoscar:

    def setUp(self):
        self.testclass = Poscar()

    def tearDown(self):
        pass

    def test_get_lattice_constant(self):
        raise SkipTest
        #self.testclass.get_lattice_constant()

    def test_get_num_atoms_line_number(self):
        raise SkipTest
        #self.testclass.get_num_atoms_line_number()

    def test_get_species_list(self):
        raise SkipTest
        #self.testclass.get_species_list()

    def test_get_num_atoms_line(self):
        raise SkipTest
        #self.testclass.get_num_atoms_line()

    def test_get_num_atoms_list(self):
        raise SkipTest
        #self.testclass.get_num_atoms_list()

    def test_get_num_atoms(self):
        raise SkipTest
        #self.testclass.get_num_atoms()

    def test_substitute_vacancy_species_string(self):
        raise SkipTest
        #self.testclass.substitute_vacancy_species_string(vac1atom)

    def test_create_new_species_string(self):
        raise SkipTest
        #self.testclass.create_new_species_string(speclist)

    def test_nearest_match(self):
        raise SkipTest
        #self.testclass.nearest_match(coords)

    def test_get_species_map(self):
        raise SkipTest
        #self.testclass.get_species_map()

