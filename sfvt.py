__author__ = 'komatsouliscy'

import sys
from Bio import AlignIO
import pandas as pd


def parse_position_input(raw_positions):
    """
    This method will take a string argument for positions and plit it out into individual start and stop positions. The
    positions can be one or more columns or integers and linear or non-linear.

    :param raw_positions:
    :return:
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


def confirm_reference_seq_in_alignment(reference_identifier, alignment, format="clustal"):
    """
    This method will test to see if the reference identifier (accession, etc) can be found in one of the alignemnt
    sequence identifiers.

    :param reference_identifier: The identifier (accession, gid) of the reference sequence
    :param alignment: The alignment to be used for the SFVT analysis
    :param format: The format of the alignment, default is clustal
    :return:
    """

    test = False
    for alignment in AlignIO.parse(alignment, format):
        for record in alignment:
            if reference_identifier in record.id:
                test = True
                reference_sequence = record.seq
    return test, reference_sequence


def confirm_sequence_feature_in_reference(reference_sequence, sequence_feature_positions):
    """
    This method will test to see if the sequence feature positions can be found in the reference sequence.

    :param reference_identifier: The identifier (accession, gid) of the reference sequence
    :param alignment: The alignment to be used for the SFVT analysis
    :param sequence_feature_positions: A list of the positions of teh variant type in the sequence feature
    :return:
    """

    test = True
    length = len(reference_sequence)
    for position in sequence_feature_positions:
        if position > length:
            test = False

    logging.info("confirm_sequence_feature_in_reference test result: %s" % test)
    return test


def check_reference_positions(reference_sequence, positions):
    """
    This function takes the aligned reference sequence and the lsit of parsed positions and checks to see if there are
    any dashes at the beginning of the sequence. If there are, they are removed and the positions are corrected.

    :param reference_sequence: The reference sequence from the alignment
    :param positions: The list of parsed positions
    :return: The original positions or adjusted positions in adjusted.
    """

    length_raw = len(reference_sequence)
    length_adjusted = len(reference_sequence.lstrip('-'))
    adjustment = length_raw - length_adjusted
    if adjustment != 0:
        positions[:] = [x - 13 for x in positions]
    return positions


def parse_directory_filename_and_extension(file_path):
    """
    This class will parse out a file path, file name and extension. The conventions used follow a PHP convention
     (apparently).

        /path/to/file.zip     # path
        /path/to              # dirname
        file.zip              # basename
        file                  # filename
        zip                   # extension

    :rtype : tuple
    :param file_path: The input realtive or full path of a file
    :return:
    """
    #TODO: Test if this works for Mac, Linux and windows

    import os.path

    if os.path.exists(file_path):

        absolute_path = os.path.abspath(file_path)
        base_name = os.path.basename(file_path)
        file_name = base_name.split(".")[0]
        file_extension = base_name.split(".")[1]
        directory_name = os.path.dirname(file_path)

        return absolute_path, directory_name, base_name, file_name, file_extension

    else:

        logging.error("File was not found: %s" % file_path)
        return None, None, None


def import_metadata(file_path):
    df_metadata = pd.read_table(file_path)#, skiprows=[0,1,2], header=0)?
    df_metadata.to_csv("df_metadata.csv")
    return df_metadata


def count_sequences_per_variant_type(dataframe, file_name):
    df_by_variant_type = dataframe.groupby('variant_type')
    count_by_variant_type = df_by_variant_type.count()
    count_by_variant_type.sort('accession', ascending=False, inplace=True)
    count_by_variant_type.to_csv("sfvt_%s.csv" % file_name)
    report = count_by_variant_type.to_string()
    logging.info("Sequences per variant type:\n%s" % report)


def main(args):

    from Bio import AlignIO

    test = False

    test, reference_sequence = confirm_reference_seq_in_alignment(args.reference_identifier, args.alignment)
    logging.info("Reference seqeunce test result: %s" % test)

    # Parse positions
    positions = parse_position_input(args.positions)
    logging.info("The positions (parsed): %s" % positions)

    test = confirm_sequence_feature_in_reference(reference_sequence, positions)

    if test:
        # read in multiple sequence alignment
        alignment = AlignIO.read(args.alignment, args.alignment_format)
        #logging.debug("The alignment is: %s" % alignment)

        checked_positions = check_reference_positions(reference_sequence, positions)

        variants = []
        for record in alignment:
            sequence = record.seq
            sequence_feature_temp = ''.join([sequence[index] for index in checked_positions ])
            variants.append([record.id, sequence_feature_temp])
            #logging.debug(sequence_feature_temp)

        null, null, null, file_name, file_extension = parse_directory_filename_and_extension(args.alignment)

        headers = ['accession', 'variant_type']
        df = pd.DataFrame(variants, columns=headers)
        df.to_csv('df_accession_index.csv')

        count_sequences_per_variant_type(df, file_name)

        if args.metadata_file is not None:
            df_metadata = import_metadata(args.metadata_file)
            df_metadata.to_csv("df_metadata.csv")
            df_all_data = pd.merge(df, df_metadata, on='accession', how='outer')
            df_all_data.to_csv("df_all_data.csv")

            columns = list(df_all_data.columns)
            print(columns)


        # for x in range(0, 3):
        #     """This is the error message for if the primary output is not what is desired."""
        #     if len(two_list[0]) == 1:
        #         print("ERROR")
        #
        # for x in range(0, 1):
        #     """This is the helpful error message for the previously mentioned problem as well as different length issues."""
        #     if len(two_list[0]) == 1:
        #         print("SEQUENCE FEATURE '%s' NOT IN ALIGNMENTS" % two_list[0][0])
        #         if len(two_list[0][0]) != len(two_list[9][0]):
        #             print("LENGTH OF SEQUENCE FEATURE DOES NOT MATCH LENGTH OF FEATURE IN ALIGNMENTS")

    else:
        logging.error("No reference identifier found: %s" % args.reference_identifier)


if __name__ == "__main__":
    import os
    import sys
    import logging
    import datetime
    import argparse

    parser = argparse.ArgumentParser(
        #prog='sfvt.py',
        #usage="%(prog)s -a alignment\n",
        #description='Scans a directory of directories for peptides that have been '
        #            'predicted and processes them, uploading results to the '
        #            'DBAASP database.',
        #formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=15),
        #add_help=False
    )
    parser.add_argument('-a', "--alignment",
                        type=str,
                        required=True,
                        help="The name (and path) of the alignment file. Format: clustal")
    parser.add_argument("-f", "--alignment_format",
                        required=False,
                        type=str,
                        default="clustal",
                        help="The alignment file format. Default = clustal")
    parser.add_argument("-r", "--reference_identifier",
                        required=True,
                        type=str,
                        help="The reference sequence identifier; An accession or gid: AB01223")
    parser.add_argument("-p", "--positions",
                        required=True,
                        type=str,
                        help="The position(s) of the sequence feature. Example: '100-110' or '100-110, 120, 130'")
    parser.add_argument("-m", "--metadata_file",
                        required=False,
                        type=str,
                        help="The metadata file (tab delimited).")
    parser.add_argument("-log", "--loglevel",
                        required=False,
                        default="debug",
                        help="What level of logging would you like to see? (info, debug, error; default=debug).")
    args = parser.parse_args()


    d = datetime.date.today()
    fh = logging.FileHandler('./%s.log' % d.isoformat())
    logging.basicConfig(stream=sys.stdout,  level=args.loglevel.upper())

    main(args)


