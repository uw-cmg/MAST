"""Tests for Outputfile"""

from MAST..outputfile import Outputfile

from nose.plugins.skip import SkipTest

class TestOutputfile:

    def setUp(self):
        self.testclass = Outputfile()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__()

