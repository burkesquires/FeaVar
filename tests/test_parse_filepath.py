from unittest import TestCase

import logging
import os

class TestParse_filepath(TestCase):

    def test_parse_filepath(self):

        from nvariant.nvariant import parse_filepath

        cwd = os.getcwd()
        test_file_path = "%s/nvariant.py" % cwd
        logging.debug(test_file_path)

        assert test_file_path != ""

        test_absolute_path, test_directory_name, test_base_name, test_file_name, \
        test_file_extension = parse_filepath(test_file_path)

        assert test_absolute_path == test_file_path
        assert test_base_name == "nvariant.py"
        assert test_file_name == "nvariant"
        assert test_file_extension == "py"
        assert test_directory_name == cwd
