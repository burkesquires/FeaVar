from unittest import TestCase

from FeaVar import FeaVar


class TestCheckReferencePositions1(TestCase):

    def test_check_adjusted_positions(self):

        test_ref_seq = "-A--BB---CCC----DDDD"

        test_positions1 = [1, 3, 5]

        adjusted_positions1 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions1)

        assert adjusted_positions1 == [2, 6, 11]


    def test_check_reference_positions(self):

        test_ref_seq = "---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFRVTPQPG-----------------" \
                       "VPPEEAGAAVAAESSTGT---------WTTVWTDGLTSLDRYKG-----RCYHIEPVPG-------------------" \
                       "EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPVAYVKTFQGP"

        test_positions1 = [1, 3, 5]

        adjusted_positions1 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions1)

        test_result = FeaVar.check_reference_positions(test_ref_seq, adjusted_positions1)

        assert test_result is True


    def test_check_reference_positions2(self):

        test_ref_seq = "---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFRVTPQPG-----------------" \
                       "VPPEEAGAAVAAESSTGT---------WTTVWTDGLTSLDRYKG-----RCYHIEPVPG-------------------" \
                       "EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPVAYVKTFQGP"

        test_positions2 = [50, 63, 75]

        adjusted_positions2 = FeaVar.adjust_positions_for_insertions(test_ref_seq, test_positions2)

        test_result = FeaVar.check_reference_positions(test_ref_seq, adjusted_positions2)

        assert test_result is True
