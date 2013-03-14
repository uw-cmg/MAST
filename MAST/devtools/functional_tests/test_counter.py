"""Tests for Counter"""

from MAST..counter import Counter

from nose.plugins.skip import SkipTest

class TestCounter:

    def setUp(self):
        self.testclass = Counter()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(i_limit)

    def test_increment(self):
        raise SkipTest
        #self.testclass.increment()

    def test_reset(self):
        raise SkipTest
        #self.testclass.reset()

    def test_count(self):
        raise SkipTest
        #self.testclass.count()

    def test_output(self):
        raise SkipTest
        #self.testclass.output()

    def test_expired(self):
        raise SkipTest
        #self.testclass.expired()

