import sys
sys.path.append('../../../MolBio')
import molbio as mb


def relative_hamming_distance(seq1: str, seq2: str) -> float:
    assert len(seq1) == len(seq2)
    distance = 0
    for s1, s2 in zip(seq1, seq2):
        if s1 != s2:
            distance += 1
    relative_distance = distance / len(seq1)
    return relative_distance


def relative_hamming_distance_matrix(reads: dict[str]) -> list[list[float]]:
    distance_matrix = []
    for read_i in reads.values():
        row = []
        for read_j in reads.values():
            row.append(relative_hamming_distance(read_i, read_j))
        distance_matrix.append(row)
    return distance_matrix


if __name__ == '__main__':
    reads = mb.parse_fasta('reads.fasta')
    dm = relative_hamming_distance_matrix(reads)
    for row in dm:
        print(*row)