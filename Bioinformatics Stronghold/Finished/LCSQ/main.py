import sys
sys.path.append('../../../MolBio')
import molbio as mb


def lcs_table_approach(sequence1: str, sequence2: str) -> list[str]:
    # THIS FUNCTION IS TOO SLOW BECAUSE OF HIGH MEMORY USAGE
    # LCS-table approach.
    # Inspired by: https://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Solution_for_two_sequences

    # No. rows (sequence 1 and null column)
    rows = len(sequence1) + 1
    # No. columns (sequence 2 and null column)
    columns = len(sequence2) + 1

    # Creating table; Outer list is rows, inner lists are columns
    lcs_table = []
    for row in range(rows):
        lcs_table.append([])
        for column in range(columns):
            lcs_table[row].append([])

    # Filling null columns
    for row in range(rows):
        for column in range(columns):
            if (row == 0) or (column == 0):
                lcs_table[row][column] = ['']

    # Evaluating cells rowwise then columnwise
    for row in range(1, rows):

        for column in range(1, columns):

            # Checking if there is a sequence match
            s_1 = sequence1[row-1]
            s_2 = sequence2[column-1]

            # If match, concat the sequence element to the upperleft cell
            if s_1 == s_2:
                upperleft = lcs_table[row-1][column-1]
                lcs_table[row][column] = [s + s_1 for s in upperleft]

            # If not match, choose the unique largest (maybe multiple) elements from the left and upper cells
            else:
                upper = lcs_table[row-1][column]
                left = lcs_table[row][column-1]
                largest_length = len(max(upper + left, key=lambda k: len(k)))
                largest = [s for s in upper + left if len(s) == largest_length]
                lcs_table[row][column] = list(set(largest))

        longest_set = len(max(lcs_table[row], key=lambda k: len(k)))

    # The last (lowerright) cell of the table are all LCS sequences
    lowerright = lcs_table[-1][-1]
    return lowerright


def lcs_traceback_table_approach(sequence1: str, sequence2: str) -> str:
    # LCS-table traceback approach.
    # Inspired by: https://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Solution_for_two_sequences
    # Inspired by: https://github.com/Hosseinem/Rosalind-solutions/blob/master/17%20%5Blcsq%5D%20Finding%20a%20Shared%20Spliced%20Motif

    # No. rows (sequence 1 and null column)
    rows = len(sequence1) + 1
    # No. columns (sequence 2 and null column)
    columns = len(sequence2) + 1

    # Creating empty table
    table = []
    for _ in range(rows):
        row = []
        table.append(row)
        for _ in range(columns):
            row.append(0)

    # Filling table with lcs max length
    for row in range(1, rows):
        for column in range(1, columns):

            # If match, take upperleft and add 1
            if sequence1[row-1] == sequence2[column-1]:
                upperleft = table[row-1][column-1]
                table[row][column] = upperleft + 1

            # If not, take largest of upper and left cell
            else:
                left = table[row][column-1]
                upper = table[row-1][column]
                table[row][column] = max(left, upper)

    # Tracing back through lengths to get sequence
    lcs = ''
    row = rows-1
    column = columns-1
    # Looping through table until Ø-border is hit, starting in lowerright corner
    while table[row][column] != 0:

        # Getting cell values
        cell = table[row][column]
        left = table[row-1][column]
        upper = table[row][column-1]
        upperleft = table[row-1][column-1]

        # If cell value matches left, then cell has inherited from this cell (or upper) because of no match; move left
        if cell == left:
            row -= 1

        # If cell value matches upper, then cell has inherited from this cell because of no match; move up
        elif cell == upper:
            column -= 1

        # If cell is equal to upperleft + 1, there was a match; append element to lcs and move to upper left cell
        elif cell - 1 == upperleft:
            lcs += sequence1[row-1]
            row -= 1
            column -= 1

        else:
            raise AssertionError("Algorithm failed!")

    # Check that the maximum length matches the final length of the lcs
    lcs_max_len = table[-1][-1]
    assert len(lcs) == lcs_max_len, "Algorithm failed!"

    # Traced sequence was written in reverse order
    lcs_lst = list(lcs)
    lcs_lst.reverse()
    lcs = ''.join(lcs_lst)

    return lcs


if __name__ == '__main__':
    sequences = mb.parse_fasta('sequences.fasta')
    lcs = lcs_traceback_table_approach(*sequences.values())
    print(lcs)