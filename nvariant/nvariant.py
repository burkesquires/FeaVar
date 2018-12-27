#!~/anaconda3/bin python
"""
nVariant

This module computes the variant type in a given sequence feature and creates
plots for each different type of metadata given.

"""

__author__ = 'R. Burke Squires'
__copyright__ = "Copyright 2018"
__credits__ = ["Carolyn Komatsoulis"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "R. Burke Squires"
__email__ = "burkesquires (at) gmail.com"
__status__ = "Beta"

# TODO Add ability to use native IEDB format (H25, H45, V46, N47, L496, S306, L307, P308, T333;
# B: D363, G364, W365, Q382, T385, Q386, I389, D390, T393, V396, N397, I400

def parse_position_input(raw_positions: str) -> list:
    """
    Takes a string argument of positions and splits it out into individual
    start and stop positions. The positions can be one or more columns or
    integers and linear or non-linear.

    Args:
        raw_positions : The raw positions of the reference sequence to assemble a sequence feature from.

    Returns:
        A list is returned

    """

    logging.debug("parse_position_input raw positions: " + raw_positions)

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
    """

    :param ref_seq:
    :return:
    """

    logging.debug("create_index_offset_list ref_seq: " + ref_seq)

    l = len(ref_seq)

    dash_indices = [pos for pos, char in enumerate(ref_seq) if char == "-"]

    # NOTE: remember that the indices begin at 0 and the sequence feature positions begin at 1

    index_correction_factor = [0] * l

    for dash_idx in dash_indices:

        for idx in range(dash_idx, l):
            index_correction_factor[idx] += 1

    logging.debug("create_index_offset_list index_correction_factor: " + index_correction_factor)

    return index_correction_factor


def correct_index_dict(ref_seq):

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


def adjust_positions_for_insertions(positions, ref_seq) -> list:
    """
    Takes a string argument of positions and the aligned reference sequence and
    corrects the sequence feature positions for any insertions in the reference sequence.

    Parameters
    ----------
    positions : string
        The positions of the reference sequence to assemble a sequence feature from.

    ref_seq : string
        The reference sequence pulled from the alignment (with insertion dashes included (if applicable))
        :rtype: list

    """
    # SFVT      [   23  45   89   ]
    # Positions [--123--456-78901-]
    # Reference [--TTA--GGG-TAGGG-]
    # Aligned   [AATTAAAGGGTTAGGGG]
    # Positions [12345678901234567]
    # Output    [   45  89   13,14]
    # dashes    [0, 1, 5, 6, 10, 17]


    # TODO: check to see if there are any other characters?, like a period or substitute characters

    # get positions of all dashes in the reference sequence [0, 1, 5, 6, 10,16]

    if positions is not None:

        corrected_indices = correct_index_dict(ref_seq)

        corrected_positions = []

        for position in positions:

            corrected_positions.append(corrected_indices[position])

        return corrected_positions

    else:

        print("No positions to adjust.")

def test_for_ref_seq_in_alignment(reference_identifier, alignment, msa_format="clustal"):
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
    from Bio import AlignIO

    reference_sequence = ""

    test = False
    for alignment in AlignIO.parse(alignment, msa_format):
        for record in alignment:
            if reference_identifier in record.id:
                test = True
                reference_sequence = record.seq
    return test, reference_sequence


def confirm_seq_feature_in_ref(reference_sequence, sequence_feature_positions):
    """
    Tests to see if the sequence feature positions can be found in the
    reference sequence.

    Parameters
    ----------
    reference_sequence : string
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

    logging.debug("confirm seq feature in reference tests result: %s" % test)
    return test


def check_reference_positions(reference_sequence: str, positions: str):
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

    # Adjust all positions based on leading gap
    if adjustment != 0:
        positions[:] = [x + adjustment for x in positions]

    # how to discover if adjusted positions are actually deletion (dashes)?
    test = True

    for position in positions:
        if reference_sequence[position - 1] == "-":
            test = False

    return test


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
    df_metadata.to_csv(os.path.join(output_dir, "df_metadata.csv"))
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
    df_by_variant_type = pd.DataFrame({'count' : dataframe.groupby(
        ["variant_type"]).size()}).reset_index()
    df_by_variant_type.sort_values('count', ascending=False, inplace=True)

    row_length = len(df_by_variant_type)
    df_by_variant_type["VT"] = [VT_count(i) for i in range(1, row_length + 1)]
    df_by_variant_type.reindex(index=["VT"])

    df_by_variant_type.to_csv(os.path.join(output_dir, "sfvt_%s.csv" % file_path))
    report = df_by_variant_type.to_string()
    logging.info("Sequences per variant type:\n%s" % report)

    return df_by_variant_type


def plot_variant_type_data(df_all_data, field):
    """

    Parameters
    ----------
    df_all_data : df_all_data
        The pandas df_all_data with all data inti

    field : string
        The metadata field to be plotted.
    """
    df_all_data.to_csv(os.path.join(output_dir, "df_by_field.csv"))
    df_by_one_field = df_all_data.groupby(['VT', field]).size()
    df_by_one_field.to_csv("df_by_%s.csv" % field)

    unpacked = df_by_one_field.unstack(level=1)
    plot = unpacked.plot(kind='bar', subplots=False)
    fig = plot.get_figure()
    fig.savefig("sfvt_%s.svg" % field)

    plot = unpacked.plot(kind='bar', stacked=True, subplots=False)
    fig = plot.get_figure()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(os.path.join(output_dir,"sfvt_stacked_%s.svg" % field, dpi=100))


def select_vts_to_plot(df, count):
    """
    Selects top variant top for plotting

    """
    vts_to_select =  ["VT-%03d" % i for i in range(count)]
    df_selected = df[df['VT'].isin(vts_to_select)]
    return df_selected


def set_output_directory(output_dir, output_file_name):
    """
    Create an OS independent path for output

    Parameters
    ----------
    output_dir : string
        The directory to save output in

    output_file_name : string
        The name of the output file to save
    """
    import os

    return os.path.join(os.path.dirname(output_dir), output_file_name)


def main(arguments):
    """
    Main method for the sequence feature variant type python script

    """

    # TODO add report of algorithm version, results, etc; look at importing pandas-html I think

    ref_seq_in_alignment = False
    reference_sequence = ""

    ref_seq_in_alignment, reference_sequence = test_for_ref_seq_in_alignment(arguments.reference_identifier, arguments.alignment)

    logging.info("Reference sequence tests result: %s" % ref_seq_in_alignment)

    logging.info("Positions : %s" % arguments.positions)
    parsed_positions = parse_position_input(arguments.positions)

    corrected_positions = adjust_positions_for_insertions(parsed_positions, reference_sequence)

    logging.info("Corrected positions: %s" % corrected_positions)

    checked_positions = check_reference_positions(reference_sequence, corrected_positions)

    rules = [len(corrected_positions) > 0,
             confirm_seq_feature_in_ref(reference_sequence, corrected_positions) == True,
             ref_seq_in_alignment == True,
             checked_positions == True]

    if all(rules):

        # read in multiple sequence alignment
        alignment = AlignIO.read(arguments.alignment, arguments.alignment_format)

        variants = []
        for record in alignment:
            sequence = record.seq
            sequence_feature_temp = ''.join([sequence[index] for index in corrected_positions])
            variants.append([record.id, sequence_feature_temp])
            # logging.debug(sequence_feature_temp)

        null, null, null, file_name, null = parse_filepath(arguments.alignment)

        headers = ['accession', 'variant_type']
        df_starter = pd.DataFrame(variants, columns=headers)
        if arguments.loglevel == 'debug': df_starter.to_csv(os.path.join(output_dir, 'df_accession_index.csv'))

        df_by_variant_type = count_seqs_per_variant_type(df_starter, file_name)
        df_by_variant_type.to_csv(os.path.join(output_dir, "variant_types.csv"))

        if arguments.metadata_file is not None:

            logging.debug("Metadata file is present at: %s" %
                          arguments.metadata_file)

            df_metadata = import_metadata(arguments.metadata_file)
            df_all_data = pd.merge(df_starter, df_metadata, on='accession',
                                   how='outer')

            if arguments.loglevel == 'debug':
                df_all_data.to_csv(os.path.join(output_dir, "df_all_data.csv"))

            df_all_data_with_variant_type = pd.merge(df_all_data,
                                                     df_by_variant_type,
                                                     on='variant_type',
                                                     how='outer')
            df_file_name = "df_all_data_with_variant_type.csv"
            df_all_data_with_variant_type.to_csv(os.path.join(output_dir, df_file_name))

            # for large dataframes select the top X rows
            df_top = select_vts_to_plot(df_all_data_with_variant_type,
                                        arguments.top)
            if arguments.loglevel == 'debug':
                df_top.to_csv(os.path.join(output_dir, "df_top.csv"))

            columns = list(df_all_data.columns)
            logging.debug("The columns in the %s dataframe are: %s" %
                          ("df_all_data", columns))

            variant_types = list(pd.unique(df_all_data.variant_type.ravel()))
            logging.debug("The variant types are: %s" % variant_types)

            for field in columns[2:]:

                logging.info("Now plotting graph for: %s" % field)
                plot_variant_type_data(df_top, field)

    else:
        logging.error("No reference identifier found: %s" % arguments.reference_identifier)


if __name__ == "__main__":
    import os
    import os.path
    import logging
    import datetime
    import argparse
    from Bio import AlignIO
    import pandas as pd

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
                        help="The position(s) of the sequence feature, enclosed in quotes, "
                             "comma separated, dashes for ranges. "
                             "Example: '100-110' or '100-110,120,130'")
    PARSER.add_argument("-d", "--project_directory",
                        required=False,
                        type=str,
                        help="The directory of the project.")
    PARSER.add_argument("-m", "--metadata_file",
                        required=False,
                        type=str,
                        help="The metadata file (tab delimited).")
    PARSER.add_argument("-t", "--top",
                        required=False,
                        type=str,
                        default=10,
                        help="The number (top) of variant types to plot (default=10).")
    PARSER.add_argument("-log", "--loglevel",
                        required=False,
                        default="debug",
                        help="What level of logging would you like to see? "
                             "(info, debug, error; default=info).")
    ARGS = PARSER.parse_args()

    # TODO: Is it better to have one log and append to it with date and time or have a log for each day or run?
    import os

    cwd = os.getcwd()
    log_directory = os.path.join(cwd, "logs")
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)

    output_dir = os.path.join(cwd, "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    D = datetime.date.today()
    f_handler = logging.FileHandler('%s/%s.log' % (log_directory, D.isoformat()))
    # logging.basicConfig(stream=sys.stdout, level=ARGS.loglevel.upper())
    logging.basicConfig(level=ARGS.loglevel.upper(),
                        format='%(name)s - %(levelname)s - %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=os.path.join(log_directory, '%s.log' % D.isoformat()),
                        filemode='w')
    logging.info("Logging started")
    logging.debug("Arguments:" + str(ARGS))

    # Create a custom logger
    logger = logging.getLogger(__name__)

    # Create handlers
    c_handler = logging.StreamHandler()
    #f_handler = logging.FileHandler('file.log')
    c_handler.setLevel(logging.WARNING)
    f_handler.setLevel(logging.ERROR)

    # Create formatters and add it to handlers
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

    main(ARGS)
