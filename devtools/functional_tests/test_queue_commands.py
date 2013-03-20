"""Tests for Queue_Commands"""

from MAST..queue_commands import Queue_Commands

from nose.plugins.skip import SkipTest

class TestQueue_Commands:

    def setUp(self):
        self.testclass = Queue_Commands()

    def tearDown(self):
        pass

    def test_submit_command(self):
        raise SkipTest
        #self.testclass.submit_command(runpath, scriptname)

    def test_queue_status_from_text(self):
        raise SkipTest
        #self.testclass.queue_status_from_text(jobid, queuetext)

    def test_extract_submitted_jobid(self):
        raise SkipTest
        #self.testclass.extract_submitted_jobid(string)

    def test_queue_snap_command(self):
        raise SkipTest
        #self.testclass.queue_snap_command()

