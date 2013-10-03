"""Tests for Phonchecker"""

from MAST..phonchecker import Phonchecker

from nose.plugins.skip import SkipTest

class TestPhonchecker:

    def setUp(self):
        self.testclass = Phonchecker()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_is_complete(self):
        raise SkipTest
        #self.testclass.is_complete()

    def test_is_ready_to_run(self):
        raise SkipTest
        #self.testclass.is_ready_to_run()

    def test__phon_poscar_setup(self):
        raise SkipTest
        #self.testclass._phon_poscar_setup()

    def test__phon_inphon_get_non_mast_keywords(self):
        raise SkipTest
        #self.testclass._phon_inphon_get_non_mast_keywords()

    def test__phon_inphon_get_allowed_keywords(self):
        raise SkipTest
        #self.testclass._phon_inphon_get_allowed_keywords(allowedpath)

    def test__phon_inphon_setup(self):
        raise SkipTest
        #self.testclass._phon_inphon_setup()

    def test__phon_inphon_get_masses(self):
        raise SkipTest
        #self.testclass._phon_inphon_get_masses()

    def test__phon_forces_setup(self):
        raise SkipTest
        #self.testclass._phon_forces_setup()

    def test__nosd_my_dynmat(self):
        raise SkipTest
        #self.testclass._nosd_my_dynmat()

    def test__replace_my_displacements(self):
        raise SkipTest
        #self.testclass._replace_my_displacements()

    def test_set_up_program_input(self):
        raise SkipTest
        #self.testclass.set_up_program_input()

    def test_is_started(self):
        raise SkipTest
        #self.testclass.is_started()

