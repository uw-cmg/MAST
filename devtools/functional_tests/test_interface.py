"""Tests for Interface"""

from MAST..interface import Interface

from nose.plugins.skip import SkipTest

class TestInterface:

    def setUp(self):
        self.testclass = Interface()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_prompt_potcar(self):
        raise SkipTest
        #self.testclass.prompt_potcar()

    def test_get_potcar_list(self):
        raise SkipTest
        #self.testclass.get_potcar_list(element,search_dir)

    def test_prompt_script_unknowns(self):
        raise SkipTest
        #self.testclass.prompt_script_unknowns(path)

    def test_retrieve_parameters(self):
        raise SkipTest
        #self.testclass.retrieve_parameters()

    def test_save_parameters(self):
        raise SkipTest
        #self.testclass.save_parameters(file)

    def test_load_parameters(self):
        raise SkipTest
        #self.testclass.load_parameters(path)

    def test_params_to_input(self):
        raise SkipTest
        #self.testclass.params_to_input(start_input)

    def test_start(self):
        raise SkipTest
        #self.testclass.start()

