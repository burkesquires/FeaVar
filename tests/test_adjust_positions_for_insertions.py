from unittest import TestCase

from FeaVar import FeaVar


class TestAdjustPositionsForInsertions(TestCase):

    def test_check_adjusted_positions1(self):

        test_ref_seq = "-A--BB---CCC----DDDD"

        test_positions1 = [1, 3, 5]

        adjusted_positions1 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions1)

        assert adjusted_positions1 == [2, 6, 11]

    def test_check_adjusted_positions_no_insertions(self):

        test_ref_seq = "ABBCCCDDDD"

        test_positions1 = [1, 3, 5]

        adjusted_positions1 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions1)

        assert adjusted_positions1 == [1, 3, 5]

    def test_check_adjusted_positions_right_insertions_only(self):

        test_ref_seq = "ABBCCCDDDD-----"

        test_positions1 = [1, 3, 5]

        adjusted_positions1 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions1)

        assert adjusted_positions1 == [1, 3, 5]

    def test_check_adjusted_positions_right_no_sequence(self):

        test_ref_seq = ""

        test_positions1 = [1, 3, 5]

        adjusted_positions1 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions1)

        assert adjusted_positions1 is None, "No reference sequence"

    def test_check_adjusted_positions_right_not_enough_sequence(self):

        test_ref_seq = "-A--BB----"

        test_positions1 = [1, 3, 5]

        adjusted_positions1 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions1)

        assert adjusted_positions1 == [2, 6, 11]
        #self.assertTrue('This is broken' in str(context.exception))

        #self.assertRaises(ExpectedException, afunction)
