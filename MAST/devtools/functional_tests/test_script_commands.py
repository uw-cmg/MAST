"""Tests for Script_Commands"""

from MAST..script_commands import Script_Commands

from nose.plugins.skip import SkipTest

class TestScript_Commands:

    def setUp(self):
        self.testclass = Script_Commands()

    def tearDown(self):
        pass

    def test_get_num_nodes(self):
        raise SkipTest
        #self.testclass.get_num_nodes(scriptpath)

    def test_get_num_procs(self):
        raise SkipTest
        #self.testclass.get_num_procs(scriptpath)

    def test_get_memory(self):
        raise SkipTest
        #self.testclass.get_memory(scriptpath)

    def test_modify_jobname(self):
        raise SkipTest
        #self.testclass.modify_jobname(scriptpath, jobname)

    def test_multiply_nodes_and_procs(self):
        raise SkipTest
        #self.testclass.multiply_nodes_and_procs(scriptpath, multiplier)

    def test_modify_walltime(self):
        raise SkipTest
        #self.testclass.modify_walltime(scriptpath, walltime)

    def test_modify_queue(self):
        raise SkipTest
        #self.testclass.modify_queue(scriptpath, queue)

