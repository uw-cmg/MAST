"""Tests for Param"""

from MAST..param import Param

from nose.plugins.skip import SkipTest

class TestParam:

    def setUp(self):
        self.testclass = Param()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(str,in_calc_list=None)

    def test_print_selection_string(self):
        raise SkipTest
        #self.testclass.print_selection_string()

    def test_add_option(self):
        raise SkipTest
        #self.testclass.add_option(str)

    def test_print_options(self):
        raise SkipTest
        #self.testclass.print_options()

    def test_output(self):
        raise SkipTest
        #self.testclass.output()

    def test_prompt_user(self):
        raise SkipTest
        #self.testclass.prompt_user()

    def test_get_option(self):
        raise SkipTest
        #self.testclass.get_option(index)

    def test_get_value(self):
        raise SkipTest
        #self.testclass.get_value()

    def test_is_int(self):
        raise SkipTest
        #self.testclass.is_int(i)

