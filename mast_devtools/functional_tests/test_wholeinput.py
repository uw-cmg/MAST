"""Tests for Wholeinput"""

from MAST..wholeinput import Wholeinput

from nose.plugins.skip import SkipTest

class TestWholeinput:

    def setUp(self):
        self.testclass = Wholeinput()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__()

    def test___setitem__(self):
        raise SkipTest
        #self.testclass.__setitem__(key, value)

    def test_print_to_file(self):
        raise SkipTest
        #self.testclass.print_to_file(filepath)

    def test_overwrite_file(self):
        raise SkipTest
        #self.testclass.overwrite_file(filepath)

    def test_print_to_screen(self):
        raise SkipTest
        #self.testclass.print_to_screen()

    def test_set_values_from_file(self):
        raise SkipTest
        #self.testclass.set_values_from_file(filepath)

    def test_reconcile_path(self):
        raise SkipTest
        #self.testclass.reconcile_path()

    def test_clone_to_new_input(self):
        raise SkipTest
        #self.testclass.clone_to_new_input(new_input)

