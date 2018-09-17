from unittest import TestCase


class TestCorrectIndexDict(TestCase):
    def test_correct_index_dict(self):

        from nvariant.nvariant import correct_index_dict

        # Do I need to test for an empy sequence?

        test_ref_seq = "-M---SP---QTE----TKAS-"
        corrected_positions = correct_index_dict(test_ref_seq)

        assert corrected_positions == {1: 2, 2: 6, 3: 7, 4: 11, 5: 12, 6: 13, 7: 18, 8: 19, 9: 20, 10: 21}

        test_ref_seq = "---MSPQTETKAS--VPPEE---WTTVWTDGL--RCYHI---VAYVKTFQGP"
        corrected_positions = correct_index_dict(test_ref_seq)

        assert corrected_positions == {1: 4, 2: 5, 3: 6, 4: 7, 5: 8, 6: 9, 7: 10, 8: 11, 9: 12, 10: 13, 11: 16,
                                       12: 17, 13: 18, 14: 19, 15: 20, 16: 24, 17: 25, 18: 26, 19: 27, 20: 28,
                                       21: 29, 22: 30, 23: 31, 24: 32, 25: 35, 26: 36, 27: 37, 28: 38, 29: 39,
                                       30: 43, 31: 44, 32: 45, 33: 46, 34: 47, 35: 48, 36: 49, 37: 50, 38: 51, 39: 52}
