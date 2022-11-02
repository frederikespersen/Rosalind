########################################################################################################################
#
#
#   `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`                               -:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
#      `=`,'=/     `=`,'=/     `=`,'=/       _  _      _  ___  o          `=`,'=/     `=`,'=/     `=`,'=/
#        y==/        y==/        y==/        )\/,) __  )) ))_) _  __        y==/        y==/        y==/
#      ,=,-<=`.    ,=,-<=`.    ,=,-<=`.     ((`(( ((_)(( ((__)(( ((_)     ,=,-<=`.    ,=,-<=`     ,=,-<=`
#    ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_                               ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
#
#
#   The MolBio module
#     - A module for data and functions related to molecular biology - i.e. information processes in the
#       central dogma and omics analysis.
#
#   Content categories
#     - Biological data       (Variables)
#     - The central dogma     (Functions)
#     - Sequencing            (Functions)
#     - Sequence analysis     (Functions)
#     - File parsing          (Functions)
#
#   ~ Frederik Espersen Knudsen, 2022
#
########################################################################################################################
# ~ TODO: Introduce macromolecule classes

import math



##### BIOLOGICAL DATA

dna_nucleobases = ['A', 'C', 'G', 'T']
rna_nucleobases = ['A', 'C', 'G', 'U']
dna_complimentary_nb = {'A': 'T',
                        'T': 'A',
                        'G': 'C',
                        'C': 'G'}
rna_complimentary_nb = {'A': 'U',
                        'U': 'A',
                        'G': 'C',
                        'C': 'G'}
purines = ['A', 'G']
pyrimidines = ['C', 'T', 'U']


amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
aas_isotropic_mass = {'A':  71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
                      'F': 147.06841, 'G':  57.02146, 'H': 137.05891, 'I': 113.08406,
                      'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
                      'P':  97.05276, 'Q': 128.05858, 'R': 156.10111, 'S':  87.03203,
                      'T': 101.04768, 'V':  99.06841, 'W': 186.07931, 'Y': 163.06333}


genetic_code = {'UUU': 'F',    'CUU': 'L',    'AUU': 'I',    'GUU': 'V',
                'UUC': 'F',    'CUC': 'L',    'AUC': 'I',    'GUC': 'V',
                'UUA': 'L',    'CUA': 'L',    'AUA': 'I',    'GUA': 'V',
                'UUG': 'L',    'CUG': 'L',    'AUG': 'M',    'GUG': 'V',
                'UCU': 'S',    'CCU': 'P',    'ACU': 'T',    'GCU': 'A',
                'UCC': 'S',    'CCC': 'P',    'ACC': 'T',    'GCC': 'A',
                'UCA': 'S',    'CCA': 'P',    'ACA': 'T',    'GCA': 'A',
                'UCG': 'S',    'CCG': 'P',    'ACG': 'T',    'GCG': 'A',
                'UAU': 'Y',    'CAU': 'H',    'AAU': 'N',    'GAU': 'D',
                'UAC': 'Y',    'CAC': 'H',    'AAC': 'N',    'GAC': 'D',
                'UAA': 'Stop', 'CAA': 'Q',    'AAA': 'K',    'GAA': 'E',
                'UAG': 'Stop', 'CAG': 'Q',    'AAG': 'K',    'GAG': 'E',
                'UGU': 'C',    'CGU': 'R',    'AGU': 'S',    'GGU': 'G',
                'UGC': 'C',    'CGC': 'R',    'AGC': 'S',    'GGC': 'G',
                'UGA': 'Stop', 'CGA': 'R',    'AGA': 'R',    'GGA': 'G',
                'UGG': 'W',    'CGG': 'R',    'AGG': 'R',    'GGG': 'G'}




##### CENTRAL DOGMA
def dna_check(strand: str) -> bool:
    """
    Checks whether a string is in a valid DNA strand format.
    * Valid format: Only containing uppercase characters A, C, T, and G

    :param strand: A DNA strand
    :return: Whether the strand is a valid DNA strand
    """
    # Finding unique characters
    chars = [*set(strand)]

    # Checking that only allowed characters are in string
    return all([char in dna_nucleobases for char in chars])


def rna_check(strand: str) -> bool:
    """
    Checks whether a string is in a valid RNA strand format.
    * Valid format: Only containing uppercase characters A, C, U, and G

    :param strand: A RNA strand
    :return: Whether the strand is a valid RNA strand
    """
    # Finding unique characters
    chars = [*set(strand)]

    # Checking that only allowed characters are in string
    return all([char in rna_nucleobases for char in chars])


def complimentary(template: str) -> str:
    """
    Takes a template DNA-strand and returns the complimentary (coding) strand.
    Both strings are in 5'-3' polarity.

    :param template: A 5'-3' DNA strand
    :return: The complimentary 5'-3' DNA strand
    """
    assert dna_check(template), "Invalid DNA strand format! (Can only contain 'A', 'C', 'G', or 'T')"

    coding = ''
    # Mapping to complimentary nucleobases and reversing polarity to 5'-3'
    for nb in template[::-1]:
        coding += dna_complimentary_nb[nb]

    return coding


def transcribe(strand: str, template=True) -> str:
    """
    Transcribes a template (or coding, see params) DNA-strand and returns the RNA strand (pre-mRNA).
    Both strings are in 5'-3' polarity.

    :param strand: A 5'-3' DNA template strand
    :param template: Whether the strand is the template strand (True) or coding (False) strand (Default 'True')
    :return: The transcribed 5'-3' RNA
    """
    # Getting coding strand
    if template:
        strand = complimentary(strand)

    # Turning coding strand to RNA
    pre_mrna = strand.replace('T', 'U')

    return pre_mrna


def orfs(mrna: str) -> list:
    """
    Identifies and returns all open reading frames substrands in a mRNA strand.
    All strands are in 5'-3' polarity
    :param mrna: A 5'-3' mRNA strand to be searched for ORFs
    :return: A list of all possible open reading frame 5'-3' mRNA substrands
    """
    assert rna_check(mrna), "Invalid RNA format! (Can only contain 'A', 'C', 'G', or 'U')"

    # Looping over all possible reading frames
    orfs = []
    for frame in [mrna[offset:] for offset in [0, 1, 2]]:

        # Splitting frame into codons, amino acids
        frame_codons = [frame[i:i + 3] for i in range(0, len(frame), 3) if i + 3 <= len(frame)]
        frame_aas = [genetic_code[codon] for codon in frame_codons]

        # Identifying orfs by start and stop codons; Assembling codon (not AA) substrings
        frame_orfs = []
        translating = False
        for aa, codon in zip(frame_aas, frame_codons):

            # Identifying start codons; Initiating an ORF
            if aa == 'M':  # Start codon; Note that there may be several starting codons in one ORF, creating intra-ORFs
                translating = True
                frame_orfs = [orf + codon for orf in frame_orfs]
                frame_orfs.append(codon)

            # Identifying start codons; Terminating active ORFs
            elif aa == 'Stop':
                translating = False
                if len(frame_orfs) != 0:
                    orfs += frame_orfs
                    frame_orfs = []

            # Adding codons to active ORFs
            else:
                if translating:
                    frame_orfs = [orf + codon for orf in frame_orfs]

    return orfs


def translate(orf: str) -> str:
    """
    Translates an open reading frame (ORF) mRNA strand into a peptide.
    The ORF mRNA must be in 5'-3' polarity.
    The peptide will be in N-C polarity.
    [ORF: As-is reading frame - all codons are translated]

    :param orf: An open reading frame of 5'-3' RNA
    :return: The translated N-C peptide
    """
    assert rna_check(orf), "ORF has invalid RNA format! (Can only contain 'A', 'C', 'T', or 'U')"
    leftover = len(orf) % 3
    assert leftover == 0, f"Sequence cannot include unfinished triplets! (sequence is {leftover} nucleobases too long)"

    # Splitting ORF into codons
    peptide = ''
    for codon, c_start in [(orf[i:i+3], i) for i in range(0, len(orf), 3)]:

        if c_start == 0:
            assert codon == 'AUG', "ORF must start with a Start codon (AUG)!"

        # Translating triplet codons by the 'universal' genetic code
        aa = genetic_code[codon]
        peptide += aa

    # Removing potential Stop codon if included
    peptide.replace('Stop', '')

    return peptide




##### SEQUENCING
def unique_reads(reads: dict) -> dict:
    """
    Takes a dict of sequence reads.
    Returns unique reads; i.e. reads that are not a subset of other reads in terms of sequence.

    :param reads: A dict with keys of read tag and values of read sequence
    :return: A similar dict to 'reads' with unique reads
    """
    # Looping over reads
    unique = {}
    for read, seq in reads.items():

        # Checking whether each read is contained within the sequence of other reads
        other_seqs = filter(lambda s: s != seq, reads.values())
        if all([seq not in s for s in other_seqs]):

            # If not contained in other read sequences, consider it unique
            unique[read] = seq

    return unique


def contig(reads: dict, rmol: float, verbose=True) -> str:
    """
    Takes a dict of substrand sequences from a contig (I.e. an overlap is assumed between strands).
    Returns a contig sequence.
    Sequences are stitched together, starting with those that have the biggest overlap.
    Sequence overlaps are checked from maximum overlap length down to a minimum overlap length.
    The relative minimum overlap length (rmol) is a float between 0 and 1.
    It defines the minimum overlap needed between strings.
    If no contig is formable, an AssertionError will be raised.
    The contig is not necessarily unique.

    ~ TODO: Rewrite loops; Reset to longest mol after each stitching to secure maximum overlap for each stitch
    ~ TODO: Evaluate chance of two strands having exact same left/right overlap given a mol; They could each be stitched

    :param reads: A dict with keys of strand tags and values of strand sequences with continuous overlap
    :param rmol: The relative minimum overlap length between strands ]0:1]
    :return: A contig stitched from the overlapping sequences.
    """

    # Checking params
    assert len(reads) > 0, "No reads provided!"
    assert 0 < rmol <= 1, "Relative minimum overlap length (rmol) must be 0 < rmol ≤ 1"

    # Extracting strand sequences
    unique = unique_reads(reads)
    non_unique = {k_n: s_n for k_n, s_n in reads.items() if k_n not in unique}
    if verbose and len(non_unique) > 0:
        print("contig(): Some reads were non-unique and were not used for contig stitching:")
        print(*[f' - {k}: {s}' for k, s in non_unique.items()], sep='\n')
    strands = list(unique.values())

    # Starting contig from the largest string
    contig = strands.pop(strands.index(max(strands)))

    # Setting temporary placeholder values for number of strands left (General loop)
    no_strands_i = len(strands)
    no_strands_j = 0

    # Stitching while there are strand lefts
    while len(strands) > 0:

        # Check for progress in general loop; Else stitching is not possible
        if no_strands_i == no_strands_j:
            non_stitching_reads = [k for k, y in unique.items() if y in strands]
            err_1 = "Strand sequences do not overlap sufficiently to stitch a contig!\n"
            if verbose:
                err_2 = "The following sequences cannot be stitched into the contig:\n"
                err_3 = "\n".join([" - " + k for k in non_stitching_reads])
                err_4 = f"\nConsider lowering rmol (Current rmol = {rmol}) or consider whether strings are part of contig."
                err_5 = " (See config_rmol_scan())"
            else:
                err_2, err_3, err_4, err_5 = [""]*4
            raise AssertionError(err_1 + err_2 + err_3 + err_4 + err_5)

        # Checking number of strands left before doing full general loop
        no_strands_i = len(strands)

        # Defining minimum overlap length (mol) for the general round
        largest_strand = max(strands)
        mol = math.floor(len(largest_strand) * rmol) + 1

        # Looping over descending overlap lengths
        for l in [len(largest_strand) - i for i in range(len(largest_strand) - mol + 1)]:

            # Setting temporary placeholder values for number of strings (Length loop)
            no_strands_m = len(strands)
            no_strands_n = 0

            # Looping over this overlap length while progress is made
            while no_strands_m > no_strands_n:

                # Checking number of strands left before loop over this length
                no_strands_m = len(strands)

                # Looping over (sufficiently long) strands
                for strand in [s for s in strands if len(s) >= l]:

                    # Defining potential overlap areas
                    str_left = strand[:l]
                    sup_left = contig[:l]
                    str_right = strand[len(strand)-l:]
                    sup_right = contig[len(contig)-l:]

                    # Checking for overlap between strand and contig ends
                    if str_left == sup_right:  # Strand-right = Superstrand-left

                        # Extending contig with non-overlap sequence
                        contig = contig + strand[l:]
                        strands.remove(strand)

                    elif sup_left == str_right: # Strand-left = Superstrand-right

                        # Extending contig with non-overlap sequence
                        contig = strand[:len(strand)-l] + contig
                        strands.remove(strand)

                    # Checking number of strands left after loop over this length
                    no_strands_n = len(strands)

        # Checking number of strands left after doing full general loop
        no_strands_j = len(strands)

    return contig


def contig_rmol_scan(reads: dict, p_stepsize=2) -> (str, float):
    """
    Scans the contig() function for the first succesfull rmol value.
    Scanning occurs with a set stepsize.
    See contig() for requirements for 'reads'.
    Note that though a rmol value may be returned, it may be insignificant;
    A small rmol may correspond to just a 1 bp overlap.
    Compare the final rmol to the length of the shortest read for a worst-case-scenario reference.

    ~TODO: Evaluate chance of two strands having exact same left/right overlap given a rmol respective to genome sizes

    :param reads: See config() for requirements
    :param p_stepsize: The negative log10 stepsize for rmol scan (default 2)
    :return: A tuple of the stitched contig and max rmol value
    """
    # Defining number of steps in rmol scan
    assert type(p_stepsize) == int and p_stepsize > 0, "p_stepsize must be a positive integer!"
    n_scan = math.floor(1 // 10**(-p_stepsize)) + 1

    # Looping over rmol values
    for rmol in [i / n_scan for i in range(n_scan)]:
        # Checking if rmol value is low enough to assemble contig
        try:
            c = contig(reads, rmol)
            return c, rmol

        except AssertionError:
            continue

    raise Exception("Parameter scan of rmol yielded no results! Strands most likely do not have a continuous overlap.")



##### SEQUENCE ANALYSIS
def peptide_mass(peptide: str) -> float:
    """
    Calculates the mass of a peptide using isotropic amino acid masses.

    :param peptide: The amino acid sequence of the peptide
    :return: The mass of the peptide in Da
    """
    # Summing amino acid masses
    mass = 0
    for aa in peptide:
        assert aa in aas_isotropic_mass, f'Peptide must only include natural amino acids! ({aa} given)'
        mass += aas_isotropic_mass[aa]

    # Rounding final mass (Ala as reference)
    uncertainty = len(str(aas_isotropic_mass['A']).split('.')[1])
    return round(mass, uncertainty)


def locate_continuous_motif(string: str, motif: str) -> list:
    """
    Locates and returns the starting positions of a motif in a DNA/RNA/Peptide string.
    String and motif polarity must match.

    :param string: A DNA/RNA/Peptide string to search in
    :param motif: A -||- motif to search for
    :return: 1-indexed locations of the motif in the string
    """
    assert len(string) >= len(motif), f"Motif length longer than string length! ({len(motif)} > {len(string)})"

    # Looping over all possible motif locations
    motif_locations = []
    for i in range(len(string) - len(motif) + 1):

        # Adding absolute motif matches
        if string[i:i + len(motif)] == motif:
            motif_locations.append(i + 1)  # String position is 1-indexed

    return motif_locations


def largest_common_motif(strands: list) -> list:
    """
    Finds the longest common motif among the given DNA/RNA/Peptide strands.
    The longest common motif may not be unique; A list of all largest common motifs of the same length is returned.

    :param strands: A list of DNA/RNA/Peptide strands to search among
    :return: All largests common motifs
    """
    # We need only check the substrings of one strand, since any substrings not contained herein are not common
    # We choose the shortest strand as to check the fewest amount of substrings
    short = min(strands)

    # Generating substrings, starting with the largest
    subs = []
    for length in reversed(range(len(short))):
        for i in range(len(short) - length + 1):
            sub = short[i:i + length + 1]

            # Checking if substring is universally common
            if all([sub in strand for strand in strands]):
                subs.append(sub)

        # If at least one match is found, return all motifs found with that length
        if len(subs) != 0:
            break

    # Returning unique motifs
    motifs = [*set(subs)]
    return motifs


def palindromes(strand: str, min_length=4, max_length=12) -> list:
    """
    Locates palindromes in a DNA strand of length within set boundaries.
    Returns a list of tuples of string 1-indexed DNA-positions and palindrome lengths

    :param strand: A DNA strand to search
    :param min_length: The minimum palindrome length (must be even, default 4)
    :param max_length: The maximum palindrome length (must be even, default 12)
    :return: A list of tuples of 1-indexed positions and lengths
    """
    assert dna_check(strand), "Invalid DNA strand format! (Can only contain 'A', 'C', 'G', or 'T')"
    assert min_length % 2 == 0 and max_length % 2 == 0, "Palindrome lengths must be even!"

    # Looping over all palindrome lengths (must be even)
    palins = []
    for length in range(min_length, max_length + 1, 2):
        half = length // 2

        # Inspecting each substring of given length
        for i in range(len(strand)-length+1):
            substrand = strand[i:i+length]

            # Dividing substring down the middle into 5' part and 3' part
            substrand5 = substrand[:half]
            substrand3 = substrand[half:]

            # A palindrome must have the first half of the strand be complimentary to the second
            if substrand5 == complimentary(substrand3):
                # Returning palindrome position (1-indexed) and length
                palins.append((i + 1, length))

    return palins




##### FILE HANDLING

def parse_fasta(filename: str) -> dict:
    """
    Parses a FASTA file to a dict.

    :param filename: Relative filename of FASTA file
    :return: A dict of with FASTA tags as keys and FASTA strings as values
    """
    # Reading file
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines()]

    # Parsing file tags and strings
    strings = {}
    for line in lines:
        if line[0] == '>':  # FASTA tag
            k = line[1:]
            strings[k] = ''
        else:  # FASTA (Sub)string
            strings[k] += line

    return strings




########################################################################################################################

def credits():
    print("""
                                                 Performed using
    `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`                               -:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
       `=`,'=/     `=`,'=/     `=`,'=/       _  _      _  ___  o          `=`,'=/     `=`,'=/     `=`,'=/  
         y==/        y==/        y==/        )\/,) __  )) ))_) _  __        y==/        y==/        y==/   
       ,=,-<=`.    ,=,-<=`.    ,=,-<=`.     ((`(( ((_)(( ((__)(( ((_)     ,=,-<=`.    ,=,-<=`     ,=,-<=`  
     ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_                               ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
                                             Frederik Espersen, 2022
    """)
