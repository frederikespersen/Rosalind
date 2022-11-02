def locate_first_spliced_motif(strand: str, motif: str) -> list:
    locations = []
    _motif = motif
    for i, nb in enumerate(strand):
        if nb == _motif[0]:
            locations.append(i+1)
            if len(_motif) == 1:
                break
            _motif = _motif[1:]

    assert len(motif) == len(locations)
    return locations


import sys
sys.path.append('../..')
import molbio as mb

if __name__ == '__main__':
    strands = mb.parse_fasta('strands.fasta')
    locations = locate_first_spliced_motif(*strands.values())
    print(*locations)
