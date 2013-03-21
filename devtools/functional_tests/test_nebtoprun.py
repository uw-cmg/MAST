"""Tests for Nebtoprun"""

from MAST..nebtoprun import Nebtoprun

from nose.plugins.skip import SkipTest

class TestNebtoprun:

    def setUp(self):
        self.testclass = Nebtoprun()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(i_input_file="")

    def test_set_up_image_runs(self):
        raise SkipTest
        #self.testclass.set_up_image_runs()

    def test_update_status(self):
        raise SkipTest
        #self.testclass.update_status()

    def test_is_complete(self):
        raise SkipTest
        #self.testclass.is_complete()

    def test_missing_essential_files(self):
        raise SkipTest
        #self.testclass.missing_essential_files()

    def test_missing_restart_files(self):
        raise SkipTest
        #self.testclass.missing_restart_files()

    def test_archive_output(self):
        raise SkipTest
        #self.testclass.archive_output()

    def test_recycle_contcar(self):
        raise SkipTest
        #self.testclass.recycle_contcar()

