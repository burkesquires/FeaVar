__author__ = 'komatsouliscy'


import sys

def main(argv):

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("afile", type=argparse.FileType("r"), help="input the name of an alignment file. format: example.example")
    parser.add_argument("rsfile", type=argparse.FileType("r"), help="input reference sequence file name. format: example.example")
    parser.add_argument("positions", type=str, help="input the positions for the ref. seq.")
    parser.add_argument("alignpositions", type=str, help="input the positions for the alignments")
    args = parser.parse_args()

    first_slice1 = []
    first_slice2 = []
    second_slice1 = []
    second_slice2 = []
    iterator = 0
    the_positions = []
    literal_positions = args.positions.split(",")
    while iterator <= (len(literal_positions)-1):
        """Splits the positions input and adds splits to a list. Example: 12-13,14 -> [ [12-13], [14] ]."""
        actual_positions = literal_positions[iterator].split("-")
        the_positions.append(actual_positions)
        iterator += 1

    list_positions = 0
    while list_positions < len(the_positions):
        """
        If the input is linear, start and end positions are added to a separate list. If not, the end positions are created.
        Then, the start and end positions are added to a separate list.
        """
        if len(the_positions[list_positions]) == 2:
            first_slice1.append(int(the_positions[list_positions][0]))
            first_slice2.append(int(the_positions[list_positions][1]))
            list_positions += 1
        elif len(the_positions[list_positions]) == 1:
            second_slice1.append(int(the_positions[list_positions][0]))
            second_slice2.append(int(int((the_positions[list_positions][0])) + 1))
            list_positions += 1
        elif list_positions > len(the_positions):
            break

    a_iterator = 0
    a_the_positions = []
    aliteral_positions = args.alignpositions.split(",")
    while a_iterator <= (len(aliteral_positions)-1):
        """Splits the alignment sequences positions."""
        aactual_positions = aliteral_positions[a_iterator].split("-")
        a_the_positions.append(aactual_positions)
        a_iterator += 1

    a_first_slice1 = []
    a_first_slice2 = []
    a_second_slice1 = []
    a_second_slice2 = []
    list_positions2 = 0
    while list_positions2 < len(a_the_positions):
        """Creates end positions for non linear and separates start and end for linear input. This is for the alignments."""
        if len(a_the_positions[list_positions2]) == 2:
            a_first_slice1.append(int(a_the_positions[list_positions2][0]))
            a_first_slice2.append(int(a_the_positions[list_positions2][1]))
            list_positions2 += 1
        elif len(a_the_positions[list_positions2]) == 1:
            a_second_slice1.append(int(a_the_positions[list_positions2][0]))
            a_second_slice2.append(int(int((a_the_positions[list_positions2][0])) + 1))
            list_positions2 += 1
        elif list_positions2 > len(a_the_positions):
            break

    from Bio.Alphabet import IUPAC
    from Bio import SeqIO

    record = SeqIO.read(args.rsfile, "fasta", IUPAC.extended_protein)

    # Already Aligned Sequences
    from Bio import AlignIO

    my_accession = []

    sequences = []

    align = AlignIO.read(args.afile, "clustal")
    for record in align:
        """Adds all the sequences and their accession numbers to their respective lists."""
        sequences.append(record.seq)
        my_accession.append(record.id)

    original_2D_list = []   # this is the list with ALL of the sequences and accession numbers

    original_2D_list.append([])

    original_2D_list.append([])

    original_2D_list[0].extend(my_accession)  # the accessions are at index 0

    original_2D_list[1].extend(sequences)     # the sequences are at index 1

    two_list = []
    two_list.append([])  # 0
    two_list.append([])  # 1 for counts
    two_list.append([])  # for accession in 2

    record = SeqIO.read("ProteinFastaResults.fasta", "fasta", IUPAC.extended_protein)
    ref_seq = []
    ref_seq.append(record.seq)
    actual_seq = ref_seq[0]

    lp = 0  # list position
    slp = 0  # second list position
    second_part = ''
    third_part = ''
    sfirst_part = ''  # second first part
    x = 0
    for item in ref_seq:
        """Slices the reference sequence based on the processed input."""
        while lp < len(first_slice1):
            first_part = actual_seq[first_slice1[lp]:first_slice2[lp]]
            second_part += first_part
            lp += 1
            if lp >= len(first_slice1):
                x += 1
                if len(second_slice1) == 0:
                    x = 2
            elif len(first_slice1) <= 1 and len(second_slice1) == 0:
                x = 2
        while slp < len(second_slice1):
            sfirst_part = actual_seq[second_slice1[slp]:second_slice2[slp]]
            third_part += sfirst_part
            slp += 1
            if slp >= len(second_slice1):
                x += 1
                if len(second_slice1) <= 1 and len(first_slice1) == 0:
                    x = 2
            elif len(second_slice1) == 0:
                x += 1
        if x == 2:
            two_list[0].append((second_part + third_part))
            two_list[1].append(1)
            two_list[2].append(record.id)
            lp = 0
            slp = 0
            second_part = ''
            third_part = ''
            sfirst_part = ''
            if slp >= len(second_slice1) and lp >= len(first_slice1):
                break

    counter = 0
    while counter <= 101:
        """Adds lists to the 2D list. These will be needed later."""
        two_list.append([])
        counter += 1

    slice_list = []
    slice_list.append([])
    slice_list.append([])

    lp1 = 0
    slp1 = 0
    the_third_part = ''
    t_third_part = ''
    the_second_part = ''
    t_second_part = ''
    y = 0
    item = 0
    while item < 40:
        """Slices the alignments sequences based on the processed input."""
        while lp1 < len(a_first_slice1):
            the_second_part = original_2D_list[1][item][a_first_slice1[lp1]:a_first_slice2[lp1]]
            the_third_part += the_second_part
            lp1 += 1
            if lp1 >= len(a_first_slice1):
                y += 1
                if len(a_first_slice1) <= 1 and len(a_second_slice1) == 0:
                    y = 2
        while slp1 < len(a_second_slice1):
            t_second_part = original_2D_list[1][item][a_second_slice1[slp1]:a_second_slice2[slp1]]
            t_third_part += t_second_part
            slp1 += 1
            if slp1 >= len(a_second_slice1):
                y += 1
                if len(a_second_slice1) <= 1 and len(a_first_slice1) == 0:
                    y = 2
            elif len(a_second_slice1) <= 1:
                y += 1
        if y == 2:
            slice_list[0].append(the_third_part + t_third_part)
            slice_list[1].append(original_2D_list[0][item])
            the_third_part = ''
            the_second_part = ''
            t_third_part = ''
            t_second_part = ''
            lp1 = 0
            slp1 = 0
            item += 1
            y = 0
            if slp1 >= len(a_second_slice1) and lp1 >= len(a_first_slice1):
                break

    count = 0
    if count <= 41:
        for item in sequences:
            """
            If the reference sequences matches an alignment sequence, then the alignment sequence is added to the same
            list as the reference sequence.
            """
            if str(two_list[0][0]).find(str(slice_list[0][count])) == 0:
                two_list[0].append(slice_list[0][count])
                two_list[2].append(slice_list[1][count])
                two_list[1].append(1)
                count += 1
            elif str(two_list[0][0]).find(str(slice_list[0][count])) == -1:
                two_list[3].append(slice_list[0][count])
                # adds all the rest of the sequences to a list to be sorted
                two_list[5].append(slice_list[1][count])           # by later code
                count += 1

    another_item_finder = []
    for x in range(0, len(slice_list[0])+1):
        "Adds numbers over and over again to a list that will be used to run through the sequences."
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(0)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(0)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(0)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(x)
        another_item_finder.append(0)
        another_item_finder.append(x)

    item_finder = 0  # goes through the list of sequences that need to be sorted
    literal_number = 0  # is the same as item_finder but actually is incremented
    literal_number += 1  # also needs to be 1 - cannot be 0!!
    searcher = 0
    searcher += 1  # keeps it at two_list[3][1] in terms of comparing that to two_list[3][#]
    sequence_list = 0    # accesses/appends to the list that has the sequences
    count_list = sequence_list + 1    # accesses/appends to the list with the strain count
    accession_number_list = count_list + 1   # accesses/appends to the list with the accession numbers
    sequence_list += 9
    count_list += 9
    accession_number_list += 9
    while sequence_list <= 90 and count_list <= 91 and accession_number_list <= 92 and len(two_list[3]) > 0 and literal_number > 0:
        """This clusters the rest of the sequences."""
        for item in two_list[3]:
            if str(two_list[3][searcher]).find(str(two_list[3][item_finder])) == 0:
                two_list[sequence_list].append(two_list[3][item_finder])
                two_list[3].remove(two_list[3][item_finder])
                two_list[count_list].append(1)
                two_list[accession_number_list].append(two_list[5][item_finder])
                two_list[5].remove(two_list[5][item_finder])
                if literal_number < len(two_list[3]):
                    if str(two_list[sequence_list][0]).find(str((two_list[3][literal_number]))) == 0:
                        two_list[sequence_list].append(two_list[3][literal_number])
                        two_list[3].remove(two_list[3][literal_number])
                        two_list[count_list].append(1)
                        two_list[accession_number_list].append(two_list[5][literal_number])
                        two_list[5].remove(two_list[5][literal_number])
                        literal_number += 1
                    elif str(two_list[sequence_list][0]).find(str(two_list[3][literal_number])) == -1:
                        literal_number += 1
                    elif literal_number >= len(two_list[3]):
                        literal_number = 0
                        sequence_list += 3
                        count_list += 3
                        accession_number_list += 3
            elif str(two_list[3][searcher]).find(str(two_list[3][item_finder])) == -1:
                two_list[sequence_list].append(two_list[3][item_finder])
                two_list[3].remove(two_list[3][item_finder])
                two_list[count_list].append(1)
                two_list[accession_number_list].append(two_list[5][item_finder])
                two_list[5].remove(two_list[5][item_finder])
                for number in another_item_finder:
                    if number in range(0, len(two_list[3])):
                        if str(two_list[sequence_list][0]).find(str((two_list[3][number]))) == 0:
                            two_list[sequence_list].append(two_list[3][number])
                            two_list[3].remove(two_list[3][number])
                            two_list[count_list].append(1)
                            two_list[accession_number_list].append(two_list[5][number])
                            two_list[5].remove(two_list[5][number])
                            number -= 2

                else:   # this is adding it every time the number in the list is above the length of the list
                    sequence_list += 3  # - leads to irregular iteration (when it is in elif)
                    count_list += 3
                    accession_number_list += 3
        if len(two_list[3]) == 1:
            searcher = 0
        elif len(two_list[3]) == 0:
            break

    print("")
    for x in range(0, 100, 3):
        if len(two_list[0]) > 1:
            """Prints the lengths of the list."""
            if x == 0:
                print("The following numbers are the lengths of the VT sorted lists:")
                print(len(two_list[x]))
            elif x > 0:
                print(len(two_list[x]))

    for x in range(0, 33, 3):
        if len(two_list[0]) > 1:
            """This list comprehension prints the results"""
            print("\nVT-%s." % int(((x / 3) + 1)))
            print("sequences: %s" % two_list[x])
            print("strain count: %s" % len(two_list[x + 1]))
            print("accession numbers: %s" % two_list[x + 2])


    for x in range(0, 3):
        """This is the error message for if the primary output is not what is desired."""
        if len(two_list[0]) == 1:
            print("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR"
                   " ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR "
                   "ERROR ERROR ERROR ERROR")

    for x in range(0, 1):
        """This is the helpful error message for the previously mentioned problem as well as different length issues."""
        if len(two_list[0]) == 1:
            print("SEQUENCE FEATURE '%s' NOT IN ALIGNMENTS" % two_list[0][0])
            if len(two_list[0][0]) != len(two_list[9][0]):
                print("LENGTH OF SEQUENCE FEATURE DOES NOT MATCH LENGTH OF FEATURE IN ALIGNMENTS")


if __name__ == "__main__":
    main(sys.argv)