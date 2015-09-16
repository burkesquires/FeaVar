#!/usr/bin/env python
"""

This module computes the variant type in a given sequence feature and creates
plots for each different type of metadata given.

"""
from Bio import AlignIO
from pandas import DataFrame
import pandas as pd

__authors__ = 'R. Burke Squires, Carolyn Komatsoulis'


def parse_position_input(raw_positions):
    """
    Takes a string argument for positions and splits it out into individual
    start and stop positions. The positions can be one or more columns or
    integers and linear or non-linear.

    Parameters
    ----------
    raw_positions : string
        The raw positions of the reference sequence to assemble a sequence
        feature from.

    """

    position_groupings = raw_positions.split(",")

    positions_coordinates = []
    for position_grouping in position_groupings:
        positions = position_grouping.split("-")
        if len(positions) == 2:
            temp_list = list(range(int(positions[0]), int(positions[1]) + 1))
            positions_coordinates += temp_list
        elif len(positions) == 1:
            positions_coordinates.append(int(positions[0]))

    return positions_coordinates


def confirm_ref_seq_in_alignment(reference_identifier, alignment,
                                 msa_format="clustal"):
    """
    Tests to see if the reference identifier (accession, etc) can be found in
    one of the alignment sequence identifiers.

    Parameters
    ----------
    reference_identifier : string
        The identifier (accession, gid) of the reference sequence

    alignment :
        The alignment to be used for the SFVT analysis

    msa_format : string
        The msa_format of the alignment, default is clustal
    """
    reference_sequence = ""
    test = False
    for alignment in AlignIO.parse(alignment, msa_format):
        for record in alignment:
            if reference_identifier in record.id:
                test = True
                reference_sequence = record.seq
    return test, reference_sequence


def confirm_seq_feature_in_ref(reference_sequence,
                               sequence_feature_positions):
    """
    Tests to see if the sequence feature positions can be found in the
    reference sequence.

    Parameters
    ----------
    reference_identifier : string
        The identifier (accession, gid) of the reference sequence

    alignment : string
        The alignment to be used for the SFVT analysis

    sequence_feature_positions : string
        A list of the positions of teh variant type in the sequence feature
    """

    test = True
    length = len(reference_sequence)
    for position in sequence_feature_positions:
        if position > length:
            test = False

    logging.debug("confirm seq feature in reference test result: %s" % test)
    return test


def check_reference_positions(reference_sequence, positions):
    """
    Takes the aligned reference sequence and the list of parsed positions
    and checks to see if there are any dashes at the beginning of the sequence.
    If there are, they are removed and the positions are corrected.

    Parameters
    ----------
    reference_sequence : string
        The reference sequence from the alignment

    positions : string
        The list of parsed positions
    """

    length_raw = len(reference_sequence)
    length_adjusted = len(reference_sequence.lstrip('-'))
    adjustment = length_raw - length_adjusted
    if adjustment != 0:
        positions[:] = [x - 13 for x in positions]
    return positions


def parse_filepath(file_path):
    """
    Parses a file path, file name and extension. The conventions used follow a
    PHP convention (apparently).

        /path/to/file.zip     # path
        /path/to              # dirname
        file.zip              # basename
        file                  # filename
        zip                   # extension

    Parameters
    ----------
    file_path : string
        The input relative or full path of a file
    """
    import os.path

    if os.path.exists(file_path):

        absolute_path = os.path.abspath(file_path)
        base_name = os.path.basename(file_path)
        file_name = base_name.split(".")[0]
        file_extension = base_name.split(".")[1]
        directory_name = os.path.dirname(file_path)

        return absolute_path, directory_name, base_name, file_name, \
               file_extension

    else:

        logging.error("File was not found: %s" % file_path)
        return None, None, None, None, None


def import_metadata(file_path):
    """
    Import delimited file with metadata for each sequence by accession number

    Parameters
    ----------
    file_path : string
        The file path of the metadata file to be imported.
    """
    df_metadata = pd.read_table(file_path)  # , skiprows=[0,1,2], header=0)?
    df_metadata.to_csv("df_metadata.csv")
    return df_metadata


def VT_count(i):
    vt_id = "VT-%03d" % (i,)
    return vt_id


def count_seqs_per_variant_type(dataframe, file_path):
    """
    Counts sequences per variant type

    Parameters
    ----------
    dataframe : dataframe
        The pandas dataframe with all data inti

    file_path : string
        The file path of the output file to be saved.
    """
    df_by_variant_type = DataFrame({'count' : dataframe.groupby(
        ["variant_type"]).size()}).reset_index()
    df_by_variant_type.sort('count', ascending=False, inplace=True)

    row_length = len(df_by_variant_type)
    df_by_variant_type["VT"] = [VT_count(i) for i in range(1, row_length + 1)]
    df_by_variant_type.reindex(index=["VT"])

    df_by_variant_type.to_csv("sfvt_%s.csv" % file_path)
    report = df_by_variant_type.to_string()
    logging.info("Sequences per variant type:\n%s" % report)

    return df_by_variant_type


def plot_variant_type_data(df_all_data, field):
    """

    Parameters
    ----------
    reference_identifier : string
        The identifier (accession, gid) of the reference sequence

    Parameters
    ----------
    df_all_data : df_all_data
        The pandas df_all_data with all data inti

    field : string
        The metadata field to be plotted.
    """
    df_all_data.to_csv("df_by_field.csv")
    df_by_one_field = df_all_data.groupby(['VT', field]).size()
    df_by_one_field.to_csv("df_by_%s.csv" % field)

    unpacked = df_by_one_field.unstack(level=1)
    plot = unpacked.plot(kind='bar', subplots=False)
    fig = plot.get_figure()
    fig.savefig("sfvt_%s.svg" % field)

    plot = unpacked.plot(kind='bar', stacked=True, subplots=False)
    fig = plot.get_figure()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig("sfvt_stacked_%s.svg" % field, dpi=100)

def select_vts_to_plot(df, count):
    """
    Selects top variant top for plotting

    """
    vts_to_select =  ["VT-%03d" % i for i in range(count)]
    df_selected = df[df['VT'].isin(vts_to_select)]
    return df_selected


def main(arguments):
    """
    Main method for the sequence feature variant type python script

    """
    test, reference_sequence = confirm_ref_seq_in_alignment(
        arguments.reference_identifier, arguments.alignment)
    logging.info("Reference seqeunce test result: %s" % test)

    # Parse positions
    positions = parse_position_input(arguments.positions)
    logging.info("The positions (parsed): %s" % positions)

    test = confirm_seq_feature_in_ref(reference_sequence, positions)

    if test:
        # read in multiple sequence alignment
        alignment = AlignIO.read(arguments.alignment,
                                 arguments.alignment_format)

        checked_positions = check_reference_positions(reference_sequence,
                                                      positions)

        variants = []
        for record in alignment:
            sequence = record.seq
            sequence_feature_temp = ''.join([sequence[index] for index in
                                             checked_positions])
            variants.append([record.id, sequence_feature_temp])
            # logging.debug(sequence_feature_temp)

        null, null, null, file_name, null = parse_filepath(arguments.alignment)

        headers = ['accession', 'variant_type']
        df_starter = pd.DataFrame(variants, columns=headers)
        if arguments.loglevel == 'debug':
            df_starter.to_csv('df_accession_index.csv')

        df_by_variant_type = count_seqs_per_variant_type(df_starter, file_name)
        df_by_variant_type.to_csv("variant_types.csv")

        if arguments.metadata_file is not None:

            logging.debug("Metadata file is present at: %s" %
                          arguments.metadata_file)

            df_metadata = import_metadata(arguments.metadata_file)
            df_all_data = pd.merge(df_starter, df_metadata, on='accession',
                                   how='outer')
            if arguments.loglevel == 'debug':
                df_all_data.to_csv("df_all_data.csv")

            df_all_data_with_variant_type = pd.merge(df_all_data,
                                                     df_by_variant_type,
                                                     on='variant_type',
                                                     how='outer')
            df_file_name = "df_all_data_with_variant_type.csv"
            df_all_data_with_variant_type.to_csv(df_file_name)

            # for large dataframes select the top X rows
            df_top = select_vts_to_plot(df_all_data_with_variant_type,
                                        arguments.top)
            if arguments.loglevel == 'debug':
                df_top.to_csv("df_top.csv")

            columns = list(df_all_data.columns)
            logging.debug("The columns in the %s dataframe are: %s" %
                          ("df_all_data", columns))

            variant_types = list(pd.unique(df_all_data.variant_type.ravel()))
            logging.debug("The variant types are: %s" % variant_types)

            for field in columns[2:]:

                logging.info("Now plotting graph for: %s" % field)
                plot_variant_type_data(df_top, field)

    else:
        logging.error("No reference identifier found: %s" %
                      arguments.reference_identifier)


if __name__ == "__main__":
    import sys
    import logging
    import datetime
    import argparse

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-a', "--alignment",
                        type=str,
                        required=True,
                        help="The name (and path) of the alignment file. "
                             "Format: clustal")
    PARSER.add_argument("-f", "--alignment_format",
                        required=False,
                        type=str,
                        default="clustal",
                        help="The alignment file format. Default = clustal")
    PARSER.add_argument("-r", "--reference_identifier",
                        required=True,
                        type=str,
                        help="The reference sequence identifier; "
                             "An accession or gid: AB01223")
    PARSER.add_argument("-p", "--positions",
                        required=True,
                        type=str,
                        help="The position(s) of the sequence feature. No"
                             "spaces are allowed. Example: '100-110' or "
                             "'100-110,120,130'")
    PARSER.add_argument("-m", "--metadata_file",
                        required=False,
                        type=str,
                        help="The metadata file (tab delimited).")
    PARSER.add_argument("-t", "--top",
                        required=False,
                        type=str,
                        default=25,
                        help="The number (top) of variant types to plot (default=25).")
    PARSER.add_argument("-log", "--loglevel",
                        required=False,
                        default="debug",
                        help="What level of logging would you like to see? "
                             "(info, debug, error; default=info).")
    ARGS = PARSER.parse_args()

    D = datetime.date.today()
    logging.FileHandler('./%s.log' % D.isoformat())
    logging.basicConfig(stream=sys.stdout, level=ARGS.loglevel.upper())

    main(ARGS)
