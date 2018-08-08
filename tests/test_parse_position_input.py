from unittest import TestCase

# "-a ../data/flu_test_small.clw -r CY010795 -p "124,142,143,268,331,332" -r CY021716"

class TestParsePositionInput(TestCase):

    def test_parse_position_input(self):

        from nvariant.nvariant import parse_position_input

        def test_parse_position_input_1():
            positions1 = '10, 21, 32, 43'
            assert parse_position_input(positions1) == [10, 21, 32, 43]

        def test_parse_position_input_2():
            positions2 = '10,21,32,43'
            assert parse_position_input(positions2) == [10, 21, 32, 43]

        def test_parse_position_input_3():
            positions3 = "10-16"
            assert parse_position_input(positions3) == [10, 11, 12, 13, 14, 15, 16]

        def test_parse_position_input_test_4():
            positions4 = "10 - 21, 32, 43"
            assert parse_position_input(positions4) == [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 32, 43]

        def test_parse_position_input_test_5():
            positions5 = "10-21"
            assert parse_position_input(positions5) == [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]

        def test_parse_position_sorted_output():
            unsorted_positions = '10, 32, 21, 43'
            assert parse_position_input(unsorted_positions) == [10, 21, 32, 43]

        # TODO: Test if a period was actually typed instead of a comma

        def test_parse_position_input_fail_empty_string():
            assert parse_position_input("") == None

        def test_parse_position_input_fail_missing_comma():
            assert parse_position_input("10,21,32 43") == None

        def test_parse_position_input_fail_missing_dash():
            assert parse_position_input("10 21,32,43") == None

        def test_parse_position_input_fail_missing_range_end():
            assert parse_position_input("10 -") == None