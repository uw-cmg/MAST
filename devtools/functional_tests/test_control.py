"""Tests for Control"""

from MAST..control import Control

from nose.plugins.skip import SkipTest

class TestControl:

    def setUp(self):
        self.testclass = Control()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__()

    def test_set_up_from_control_path(self):
        raise SkipTest
        #self.testclass.set_up_from_control_path()

    def test_filestring_to_dirstring(self):
        raise SkipTest
        #self.testclass.filestring_to_dirstring(filename)

    def test_dirstring_to_filestring(self):
        raise SkipTest
        #self.testclass.dirstring_to_filestring(dirname)

    def test_send_message(self):
        raise SkipTest
        #self.testclass.send_message(emsubject, emmessage)

    def test_get_queue_snapshot(self):
        raise SkipTest
        #self.testclass.get_queue_snapshot()

    def test_record_output(self):
        raise SkipTest
        #self.testclass.record_output(mypath, outputmsg)

    def test_reset_runlist(self):
        raise SkipTest
        #self.testclass.reset_runlist()

    def test_submit_from_runlist(self):
        raise SkipTest
        #self.testclass.submit_from_runlist()

    def test_control_start(self):
        raise SkipTest
        #self.testclass.control_start(cgpath, cginput)

    def test_control_loop(self):
        raise SkipTest
        #self.testclass.control_loop(this_calcgroup, prelim=0)

