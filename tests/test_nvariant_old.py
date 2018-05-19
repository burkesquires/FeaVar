from nvariant.nvariant import parse_position_input
import pytest


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

def test_parse_position_input_fail_empty_string():
     assert parse_position_input("") == None

def test_parse_position_input_fail_missing_comma():
    assert parse_position_input("10,21,32 43") == None

def test_parse_position_input_fail_missing_dash():
    assert parse_position_input("10 21,32,43") == None

def test_parse_position_input_fail_missing_range_end():
    assert parse_position_input("10 -") == None
