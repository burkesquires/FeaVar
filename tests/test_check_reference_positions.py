from unittest import TestCase


class TestCheckReferencePositions(TestCase):

    def test_check_reference_positions(self):

        import py.test

        from nvariant.nvariant import check_reference_positions

        test_ref_seq = "---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFRVTPQPG-----------------" \
                       "VPPEEAGAAVAAESSTGT---------WTTVWTDGLTSLDRYKG-----RCYHIEPVPG-------------------" \
                       "EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPVAYVKTFQGP"

        test_positions = [1, 3, 5]

        test_result = check_reference_positions(test_ref_seq, test_positions)

        assert test_result == True
