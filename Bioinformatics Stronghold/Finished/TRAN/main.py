import sys
sys.path.append('../..')
import molbio as mb


def ratio_transition_transversion(strand1: str, strand2: str) -> float:
    assert mb.dna_check(strand1) and mb.dna_check(strand2), "Strands must have correct DNA format!"
    assert len(strand1) == len(strand2), "Strands must be of same length!"

    transitions = 0
    transversions = 0
    for nb1, nb2 in zip(strand1, strand2):
        if nb1 != nb2:
            nb1_puri = nb1 in mb.purines
            nb2_puri = nb2 in mb.purines
            if (nb1_puri and nb2_puri) or ((not nb1_puri) and (not nb2_puri)):
                transitions += 1
            else:
                transversions += 1

    return transitions / transversions


fasta_filename = 'strands.fasta'
if __name__ == '__main__':
    strands = mb.parse_fasta(fasta_filename)
    r = ratio_transition_transversion(*strands.values())
    print(r)