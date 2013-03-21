"""Tests for Creator"""

from MAST..creator import Creator

from nose.plugins.skip import SkipTest

class TestCreator:

    def setUp(self):
        self.testclass = Creator()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(i_inputfile)

    def test_start(self):
        raise SkipTest
        #self.testclass.start()

    def test_create_kpoints(self):
        raise SkipTest
        #self.testclass.create_kpoints()

    def test_create_potcar(self):
        raise SkipTest
        #self.testclass.create_potcar()

    def test_make_bulk_start_poscar(self):
        raise SkipTest
        #self.testclass.make_bulk_start_poscar(getfromprelim="")

    def test_make_poscar_with_vac(self):
        raise SkipTest
        #self.testclass.make_poscar_with_vac()

    def test_make_standard_relax_incar(self):
        raise SkipTest
        #self.testclass.make_standard_relax_incar()

    def test_make_static_incar(self):
        raise SkipTest
        #self.testclass.make_static_incar(baseincar, writepath)

    def test_make_defected_incar(self):
        raise SkipTest
        #self.testclass.make_defected_incar(baseincar, defectpos, writepath)

    def test_make_neb_relax_incar(self):
        raise SkipTest
        #self.testclass.make_neb_relax_incar(defincar, writepath)

    def test_put_entry_in_control_folder(self):
        raise SkipTest
        #self.testclass.put_entry_in_control_folder()

    def test_write_output(self):
        raise SkipTest
        #self.testclass.write_output()

    def test_record(self):
        raise SkipTest
        #self.testclass.record(fileadded, rdir)

    def test_make_paths(self):
        raise SkipTest
        #self.testclass.make_paths()

    def test_populate_paths(self):
        raise SkipTest
        #self.testclass.populate_paths(dirsetuplist)

    def test_make_vacnumber_list(self):
        raise SkipTest
        #self.testclass.make_vacnumber_list(initposcar)

