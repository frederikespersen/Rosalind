import sys
sys.path.append('../..')
import molbio as mb


def largest_common_motifs(strands:list) -> list:

    # We need only check the substrings of one strand, since any substrings not contained herein are not common
    # We choose the shortest strand to check the fewest amount of substrings
    short = min(strands)

    # Generating substrings, starting with the largest
    subs = []
    for length in reversed(range(len(short))):
        for i in range(len(short)-length):
            sub = short[i:i+length+1]

            # Checking if substring is universally common
            if all([sub in strand for strand in strands]):
                subs.append(sub)

        # If at least match is found, don't break loop till all substrings of this length have been checked
        if len(subs) != 0:
            break

    # Returning unique motifs
    motifs = [*set(subs)]
    return motifs


fasta_filename = 'strands.fasta'
if __name__ == '__main__':
    strands = mb.parse_fasta(fasta_filename)
    lcm = largest_common_motifs([*strands.values()])
    print(lcm[0])
    print(lcm)