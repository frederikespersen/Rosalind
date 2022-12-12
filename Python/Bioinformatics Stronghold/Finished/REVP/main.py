import sys
sys.path.append('../..')
import molbio as mb


def locate_palindromes(strand: str, min_length=4, max_length=12) -> dict:
    assert mb.dna_check(strand), "Invalid DNA strand format! (Can only contain 'A', 'C', 'T', or 'G')"
    assert min_length % 2 == 0 and max_length % 2 == 0, "Palindrome lengths must be even!"

    # Looping over all palindrome lengths (must be even)
    palindromes = {}
    for length in range(min_length, max_length + 1, 2):
        half = length // 2

        # Inspecting each substring of required length
        for i in range(len(strand)-length+1):
            substrand = strand[i:i+length]

            # Dividing substring down the middle into 5' part and 3' part
            substrand5 = substrand[:half]
            substrand3 = substrand[half:]

            # A palindrome must have the first half of the string be complimentary to the second
            if substrand5 == mb.complimentary(substrand3):
                # Returning palindrome position (1-indexed) with length as value
                palindromes[str(i + 1)] = length

    return palindromes


fasta_filename = 'strand.fasta'
if __name__ == '__main__':
    strand = [*mb.parse_fasta(fasta_filename).values()][0]
    palindromes = locate_palindromes(strand)
    for position, length in palindromes.items():
        print(position, length)

