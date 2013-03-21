"""Tests for Mastfile"""

from MAST..mastfile import Mastfile

from nose.plugins.skip import SkipTest

class TestMastfile:

    def setUp(self):
        self.testclass = Mastfile()

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(file_path=str())

    def test_from_file(self):
        raise SkipTest
        #self.testclass.from_file(file_path)

    def test_from_file_append(self):
        raise SkipTest
        #self.testclass.from_file_append(file_path=str())

    def test_to_file(self):
        raise SkipTest
        #self.testclass.to_file(file_path)

    def test_to_unique_file(self):
        raise SkipTest
        #self.testclass.to_unique_file(parent_path="", try_name="", suffix="", max=10)

    def test_to_stdout(self):
        raise SkipTest
        #self.testclass.to_stdout()

    def test_get_line_number(self):
        raise SkipTest
        #self.testclass.get_line_number(line_number)

    def test_get_line_match(self):
        raise SkipTest
        #self.testclass.get_line_match(string_to_match)

    def test_get_line_match_number(self):
        raise SkipTest
        #self.testclass.get_line_match_number(string_to_match)

    def test_get_last_line_match(self):
        raise SkipTest
        #self.testclass.get_last_line_match(string_to_match)

    def test_get_segment(self):
        raise SkipTest
        #self.testclass.get_segment(string_to_chop, start_string="", end_string="")

    def test_get_segment_from_last_line_match(self):
        raise SkipTest
        #self.testclass.get_segment_from_last_line_match(string_to_match, start_string="", end_string="", match_last_line=True)

    def test_modify_file_by_line_number(self):
        raise SkipTest
        #self.testclass.modify_file_by_line_number(lineno="", mode="", param="")

    def test_file_to_dictionary(self):
        raise SkipTest
        #self.testclass.file_to_dictionary(from_path)

    def test_copy_data_to(self):
        raise SkipTest
        #self.testclass.copy_data_to(other_file)

