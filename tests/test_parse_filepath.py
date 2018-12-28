import logging
import os
from unittest import TestCase


class TestParseFilepath(TestCase):

    def test_parse_filepath(self):

        from nvariant.nvariant import parse_filepath

        cwd = os.getcwd()
        test_file_path = os.path.join(cwd, "nvariant.py")
        logging.debug(test_file_path)

        assert test_file_path != ""

        absolute_path, directory_name, base_name, file_name, file_extension = parse_filepath(test_file_path)

        assert absolute_path == test_file_path
        assert base_name == "nvariant.py"
        assert file_name == "nvariant"
        assert file_extension == "py"
        assert directory_name == cwd
