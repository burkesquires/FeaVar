#!~/anaconda3/bin python

"""
FeaVar

This module computes the variant type(s) in a given sequence feature and
creates plots for each different type of metadata given.
"""

__author__ = 'R. Burke Squires'
__copyright__ = "Copyright 2025"
__credits__ = ["Carolyn Komatsoulis"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "R. Burke Squires"
__email__ = "burkesquires (at) gmail.com"
__status__ = "Beta"

import os
import sys
import logging
import datetime
import argparse
from typing import Optional, Any, Union

import pandas
from Bio import AlignIO
from pandas import DataFrame

# Global variable placeholders (initially defined in main/pre-flight)
# It is better practice to pass these, but keeping globals to minimize logic breakage for now.
df_starter = None
df_by_variant_type = None
checked_positions = None
parsed_positions = None
corrected_positions = None
output_dir = None # Defined in run_cli, used in functions


def parse_position_input(raw_positions: str) -> list:
    """
    Takes a string argument of positions and splits it out into individual start and stop positions.
    """
    logging.debug("parse_position_input raw positions: {0}".format(raw_positions))

    position_coordinates = []
    position_groupings = raw_positions.split(",")

    if len(position_groupings) > 0:
        for position_group in position_groupings:
            positions = position_group.split("-")
            if len(positions) == 2:
                temp_list = list(range(int(positions[0]), int(positions[1]) + 1))
                position_coordinates.extend(temp_list)
            elif len(positions) == 1:
                position_coordinates.append(int(positions[0]))

    sorted_positions = sorted(position_coordinates)
    logging.debug("parse_position_input parsed positions: " + str(sorted_positions))

    return sorted_positions


def create_index_offset_list(ref_seq) -> list:
    """Adjust the indices of the reference sequence if necessary."""
    logging.debug("create_index_offset_list ref_seq: {0}".format(ref_seq))
    seq_length = len(ref_seq)
    dash_indices = [pos for pos, char in enumerate(ref_seq) if char == "-"]
    index_correction_factor = [0] * seq_length

    for dash_idx in dash_indices:
        for idx in range(dash_idx, seq_length):
            index_correction_factor[idx] += 1

    logging.debug("create_index_offset_list index_correction_factor: {0}".format(str(index_correction_factor)))
    return index_correction_factor


def correct_index_dict(ref_seq: str) -> dict:
    """Create a dictionary of corrected position for each position in a sequence"""
    corrected_index = {}
    dash_count = 0
    char_count = 1

    for idx, char in enumerate(ref_seq):
        if char == "-":
            dash_count += 1
            continue
        corrected_index[char_count] = char_count + dash_count
        char_count += 1

    return corrected_index


def adjust_positions_for_insertions(ref_seq: str, positions: list) -> Optional[list[Any]]:
    """Corrects the sequence feature positions for any insertions in the reference sequence."""
    if positions:
        if ref_seq:
            corrected_indices = correct_index_dict(ref_seq)
            return [corrected_indices[position] for position in positions]
        else:
            print("No reference sequence found.")
            return None
    else:
        print("No positions to adjust.")
        return None


def check_for_ref_seq_in_alignment(reference_identifier, alignment_path, msa_format="clustal"):
    """Tests to see if the reference identifier can be found in the alignment."""
    reference_sequence = ""
    test = False
    
    # Note: alignment passed here is the path string based on args
    for alignment in AlignIO.parse(alignment_path, msa_format):
        for record in alignment:
            if reference_identifier in record.id:
                test = True
                reference_sequence = str(record.seq)

    return test, reference_sequence


def confirm_seq_feature_in_ref(reference_sequence: str, sequence_feature_positions: list) -> bool:
    """Tests to see if the sequence feature positions can be found in the reference sequence."""
    test = True
    raw_sequence = reference_sequence.replace("-", "")
    length = len(raw_sequence)

    for position in sequence_feature_positions:
        if position > length:
            test = False

    logging.debug('Confirm seq feature in reference tests result: {}'.format(test))
    return test


def check_reference_positions(reference_sequence: str, positions: list) -> bool:
    """Checks to see if there are any dashes at the beginning of the sequence."""
    seq_length = len(reference_sequence)
    test = True

    if all(i <= seq_length for i in positions):
        for position in positions:
            if reference_sequence[position - 1] == "-":
                test = False
        return test
    else:
        return False


def import_metadata(metadata_file_path: str) -> pandas.DataFrame:
    """Import delimited file with metadata for each sequence by accession number"""
    df_metadata = pandas.read_table(metadata_file_path)
    if output_dir:
        df_metadata.to_csv(os.path.join(output_dir, "df_metadata.csv"))
    return df_metadata


def vt_count(i):
    return "VT-%03d" % (i,)


def count_seqs_per_variant_type(dataframe: pandas.DataFrame, file_path) -> pandas.DataFrame:
    """Counts sequences per variant type"""
    df_data_by_variant_type = pandas.DataFrame({'count': dataframe.groupby(["variant_type"]).size()}).reset_index()
    df_data_by_variant_type.sort_values('count', ascending=False, inplace=True)

    row_length = len(df_data_by_variant_type)
    df_data_by_variant_type["VT"] = [vt_count(i) for i in range(1, row_length + 1)]
    # Fixed reindex usage which was not assigning back
    df_data_by_variant_type = df_data_by_variant_type.set_index("VT", drop=False)

    if output_dir:
        df_data_by_variant_type.to_csv(os.path.join(output_dir, "feavar_{}.csv".format(file_path)))
    
    report = df_data_by_variant_type.to_string()
    logging.info("Sequences per variant type:\n{}".format(report))

    return df_data_by_variant_type


def plot_variant_type_data(df_all_data: pandas.DataFrame, field: str):
    """Plots data based on metadata field."""
    if output_dir:
        df_all_data.to_csv(os.path.join(output_dir, "df_by_field.csv"))
    
    df_by_one_field = df_all_data.groupby(['VT', field]).size()
    # df_by_one_field.to_csv("df_by_{}.csv".format(field)) # Optional: save locally

    unpacked = df_by_one_field.unstack(level=1)
    
    # Simple bar plot
    plot = unpacked.plot(kind='bar', subplots=False)
    fig = plot.get_figure()
    # fig.savefig("feavar_{}.svg".format(field)) # Optional save

    # Stacked bar plot
    plot = unpacked.plot(kind='bar', stacked=True, subplots=False)
    fig = plot.get_figure()
    fig.set_size_inches(18.5, 10.5)
    
    if output_dir:
        fig.savefig(os.path.join(output_dir, "feavar_stacked_{}.svg".format(field)), dpi=100)


def select_var_types_to_plot(df: pandas.DataFrame, count: int) -> pandas.DataFrame:
    """Selects top variant for plotting"""
    # Ensure count is int
    count = int(count)
    vts_to_select = ["VT-%03d" % i for i in range(1, count + 1)] # Adjusted range to match VT-001 start
    df_selected = df[df['VT'].isin(vts_to_select)]
    return df_selected


def pre_flight_check(arguments):
    """Runs checks before main computation."""
    logging.info("Pre-flight starting.")

    global checked_positions, parsed_positions, corrected_positions

    ref_seq_in_alignment, reference_sequence = check_for_ref_seq_in_alignment(arguments.reference_identifier,
                                                                              arguments.alignment)

    if ref_seq_in_alignment:
        logging.info("Reference sequence found in alignment: {}".format(arguments.reference_identifier))
        parsed_positions = parse_position_input(arguments.positions)
        logging.info("Parsed positions: {}".format(parsed_positions))
        corrected_positions = adjust_positions_for_insertions(reference_sequence, parsed_positions)
        logging.info("Corrected positions: {}".format(corrected_positions))
        checked_positions = check_reference_positions(reference_sequence, corrected_positions)
    else:
        logging.error("No reference identifier found: {}".format(arguments.reference_identifier))
        corrected_positions = [] # ensure it exists
        checked_positions = False

    rules = [ref_seq_in_alignment,
             confirm_seq_feature_in_ref(reference_sequence, parsed_positions) if ref_seq_in_alignment else False,
             len(corrected_positions) > 0 if corrected_positions else False,
             checked_positions]

    print("Analysis complete.")
    return corrected_positions, rules


def compute_variant_types(alignment_file_path: str, alignment_format: str, vt_positions: list, log_level: str):
    """Computes the variant types from alignment."""
    global df_starter, df_by_variant_type

    try:
        alignment = AlignIO.read(alignment_file_path, alignment_format)
        variants = []
        for record in alignment:
            sequence = record.seq
            # Note: vt_positions are 1-based from arguments usually, need to ensure 0-based indexing for slicing
            # But the logic in adjust_positions seems to return 1-based indices? 
            # Assuming logic holds from original script:
            sequence_feature_temp = ''.join([sequence[index-1] for index in vt_positions]) 
            variants.append([record.id, sequence_feature_temp])

        _, file_name = os.path.split(alignment_file_path)
        headers = ['accession', 'variant_type']
        df_starter = pandas.DataFrame(variants, columns=headers)

        if log_level == 'debug' and output_dir:
            df_starter.to_csv(os.path.join(output_dir, 'df_accession_index.csv'))

        df_by_variant_type = count_seqs_per_variant_type(df_starter, file_name)
        if output_dir:
            df_by_variant_type.to_csv(os.path.join(output_dir, "variant_types.csv"))

    except OSError as err:
        print("OS error: {0}".format(err))
    except ValueError:
        print("Could not load or read alignment.")
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

    return df_by_variant_type, df_starter


def process_metadata(metadata_file, df_data_by_variant_type, df_starter, no_of_vt_to_plot, log_level):
    """Merges variant data with metadata and plots."""
    try:
        logging.debug("Metadata file is present at: {}".format(metadata_file))

        df_metadata = import_metadata(metadata_file)
        df_all_data = pandas.merge(df_starter, df_metadata, on='accession', how='outer')

        if log_level == 'debug' and output_dir:
            df_all_data.to_csv(os.path.join(output_dir, "df_all_data.csv"))

        df_all_data_with_variant_type = pandas.merge(df_all_data,
                                                     df_data_by_variant_type,
                                                     on='variant_type',
                                                     how='outer')

        if output_dir:
            df_file_name = "df_all_data_with_variant_type.csv"
            df_all_data_with_variant_type.to_csv(os.path.join(output_dir, df_file_name))

        # for large dataframes select the top X rows
        df_top = select_var_types