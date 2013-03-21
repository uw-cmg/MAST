"""Tests for Recorder"""

from MAST..recorder import Recorder

from nose.plugins.skip import SkipTest

class TestRecorder:

    def setUp(self):
        self.testclass = Recorder()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__()

    def test_record_output(self):
        raise SkipTest
        #self.testclass.record_output(mypath, outputmsg)

    def test_send_message(self):
        raise SkipTest
        #self.testclass.send_message(emsubject, emmessage)

