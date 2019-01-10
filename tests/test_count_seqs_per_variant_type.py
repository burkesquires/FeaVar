from unittest import TestCase

from nvariant.nvariant import count_seqs_per_variant_type

sequences = [
    "gb:CY021716: TCAATTATATTCAATATGGAAAGAATAAAAGAACTAAGGAATCTAATGTCGCAATCTCGC\n\
     gb:CY020292: TCAATTATATTCAATATGGAAAGAATAAAAGAACTAAGGAATCTAATGTCGCAGTCTCGC\n\
     gb:CY083917: TCAAATATATTCAATATGGAGAGAATAAAAGAACTGAGAGATGTTATGTCGCAGTCCCGC\n\
     gb:CY063613: --------------TATGGAGAGAATAAAAGAACTGAGAGAACTAATGTCGCAGTCCCGC\n\
     gb:CY083782: TCAAATATATTCAATATGGAGAGAATAAAAGAACTGAGAGATCTTATGTCGCAGTCCCGC\n\
     gb:CY073732: --------------TATGGAGAGAATAAAAGAACTGAGAGAACTAATGTCGCAGTCCCGC\n\
     gb:CY062698: --------------TATGGAGAGAATAAAAGAACTGAGAGATCTAATGTCGCAGTCCCGC\n\
     gb:CY062706: --------------TATGGAGAGAATAAAAGAACTGAGAGATCTTATGTCGCAGTCCCGC\n\
     gb:CY066942: --------------TATGGAGAGAATAAAAGAACTGAGAGATCTAATGTCGCAGTCCCGC\n\
     gb:CY066950: --------------TATGGAGAGAATAAAAGAACTGAGAGAACTTATGTCGCAGTCCCGC\n"
]

,variant_type,count,VT
3,AGAACTGAGAGATCTAATGTC,5,VT-001
1,AGAACTAAGGAATCTAATGTC,2,VT-002
0,AGAACTAAGAGATCTAATGTC,1,VT-003
2,AGAACTGAGAAATCTAATGTC,1,VT-004
4,AGAACTGAGGGATCTAATGTC,1,VT-005


region = "25 - 35"
test_ref_seq = 'CY063613'

class TestCountSeqsPerVariantType(TestCase):

    def test_count_seqs_per_variant_type(self):

        import pandas as pd

        input_df = pd.DataFrame.from_dict({
            'id': ['sequence1', 'sequence2', 'sequence3', 'sequence4', 'sequence5'],
            'variant_type': ['VT001', 'VT002', 'VT001', 'VT002', 'VT001']
        })

        d = {'variant_type': ['AT', 2], 'count': [3, 4], 'VT': [3, 4]}
        expected = pd.DataFrame(data=d)

        }


        file_path = ""

        # assert_dict_equal(expected, my_func(input).to_dict(), "oops, there's a bug...")

        df1 = pd.DataFrame({'a': [1, 2, 3, 4, 5]})
        df2 = pd.DataFrame({'a': [6, 7, 8, 9, 10]})

        expected_res = pd.Series([7, 9, 11, 13, 15])
        pd.testing.assert_series_equal((df1['a'] + df2['a']), expected_res, check_names=False)


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