"""Tests for Mastobj"""

from MAST..mastobj import Mastobj

from nose.plugins.skip import SkipTest

class TestMastobj:

    def setUp(self):
        self.testclass = Mastobj()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(allowed_keys, **kwargs)

    def test_set_keywords(self):
        raise SkipTest
        #self.testclass.set_keywords(**kwargs)

    def test_help(self):
        raise SkipTest
        #self.testclass.help(keyword)

