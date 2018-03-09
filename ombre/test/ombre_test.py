from ombre import ombre

# test entire script







# Test various position strings to make sure that they all work

def parse_position_input_test():
    positions1 = "10, 21, 32, 43"
    positions2 = "10,21,32,43"
    positions3 = "10-21,32,43"
    positions4 = "10 - 21, 32, 43"
    positions5 = "10-21"

    assert ombre.parse_position_input(positions1) == "10, 21, 32, 43"
    assert ombre.parse_position_input(positions2) == "10, 21, 32, 43"
    assert ombre.parse_position_input(positions3) == "10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 32, 43"
    assert ombre.parse_position_input(positions4) == "10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 32, 43"
    assert ombre.parse_position_input(positions5) == "10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21"


def parse_position_input_fail_test():
    fail_positions1 = ""
    fail_positions2 = "10,21,32 43"
    fail_positions3 = "10 21,32,43"
    fail_positions4 = "10 21, 32, 43"
    fail_positions5 = "10 -"

    try:
        ombre.parse_position_input(fail_positions1)
        assert False
    except ValueError:
        assert True

    try:
        ombre.parse_position_input(fail_positions2)
        assert False
    except ValueError:
        assert True

    try:
        ombre.parse_position_input(fail_positions3)
        assert False
    except ValueError:
        assert True

    try:
        ombre.parse_position_input(fail_positions4)
        assert False
    except ValueError:
        assert True

    try:
        ombre.parse_position_input(fail_positions5)
        assert False
    except ValueError:
        assert True

# Test confirm_seq_feature_in_ref






# check_reference_positions




# import_metadata




