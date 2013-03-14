"""Tests for Calcgroup"""

from MAST..calcgroup import Calcgroup

from nose.plugins.skip import SkipTest

class TestCalcgroup:

    def setUp(self):
        self.testclass = Calcgroup()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(i_calcgroup_input)

    def test_vet_calc_type(self):
        raise SkipTest
        #self.testclass.vet_calc_type(cid)

    def test_start_calcgroup(self):
        raise SkipTest
        #self.testclass.start_calcgroup()

    def test_set_up_calc_list(self):
        raise SkipTest
        #self.testclass.set_up_calc_list(cid)

    def test_rapid_update(self):
        raise SkipTest
        #self.testclass.rapid_update()

    def test_update(self):
        raise SkipTest
        #self.testclass.update()

    def test_bulk_complete(self):
        raise SkipTest
        #self.testclass.bulk_complete()

    def test_ep_complete(self):
        raise SkipTest
        #self.testclass.ep_complete()

    def test_make_vacancy_poscar_from_bulk(self):
        raise SkipTest
        #self.testclass.make_vacancy_poscar_from_bulk(initposcar, vacint)

    def test_get_endpoint_path_from_vacancy(self):
        raise SkipTest
        #self.testclass.get_endpoint_path_from_vacancy(vacint)

    def test_get_vacancy_from_endpoint_path(self):
        raise SkipTest
        #self.testclass.get_vacancy_from_endpoint_path(eppath)

    def test_make_vacancy_from_bulk(self):
        raise SkipTest
        #self.testclass.make_vacancy_from_bulk(vacint)

    def test_get_vacancies_from_neb_path(self):
        raise SkipTest
        #self.testclass.get_vacancies_from_neb_path(nebpath)

    def test_neb_complete(self):
        raise SkipTest
        #self.testclass.neb_complete()

