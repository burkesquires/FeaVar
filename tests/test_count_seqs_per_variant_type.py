from unittest import TestCase

from nvariant.nvariant import count_seqs_per_variant_type

class TestCountSeqsPerVariantType(TestCase):

    def test_count_seqs_per_variant_type(self):

        import pandas as pd

        input_df = pd.DataFrame.from_dict({
            'id': ['sequence1', 'sequence2', 'sequence3', 'sequence4', 'sequence5'],
            'variant_type': ['VT001', 'VT002', 'VT001', 'VT002', 'VT001']
        })

        input_df = pd.DataFrame.from_dict({
            'id': ['sequence1', 'sequence2', 'sequence3', 'sequence4', 'sequence5'],
            'variant_type': ['VT001', 'VT002', 'VT001', 'VT002', 'VT001']
        })

        expected = {
            'result': [...]
        }

        file_path = ""

        #assert_dict_equal(expected, my_func(input).to_dict(), "oops, there's a bug...")

        df_test = count_seqs_per_variant_type(input_df, file_path)

        TestCase.assertDictEqual(df_test, count_seqs_per_variant_type(input_df).to_dict(), "oops, there's a bug...")

    """
    df_by_variant_type = pd.DataFrame({'count' : dataframe.groupby(
        ["variant_type"]).size()}).reset_index()
    df_by_variant_type.sort_values('count', ascending=False, inplace=True)

    row_length = len(df_by_variant_type)
    df_by_variant_type["VT"] = [VT_count(i) for i in range(1, row_length + 1)]
    df_by_variant_type.reindex(index=["VT"])

    df_by_variant_type.to_csv("sfvt_%s.csv" % file_path)
    report = df_by_variant_type.to_string()
    logging.info("Sequences per variant type:\n%s" % report)

    return df_by_variant_type
    """