import sys
sys.path.append('../../../MolBio')
import molbio as mb


def error_corrections(reads: dict[str]) -> dict[str]:
    # Collecting reads together
    grouped_reads = {}
    for tag, read in reads.items():
        for read_pair in grouped_reads.keys():
            if read in read_pair:
                grouped_reads[read_pair] += [tag]
                break
        if not any([read in read_pair for read_pair in grouped_reads.keys()]):
            read_pair = f'{read}/{mb.complimentary(read)}'
            grouped_reads[read_pair] = [tag]

    # Defining error reads as cases with only 1 read instance
    grouped_correct_reads = {k: v for k, v in grouped_reads.items() if len(v) > 1}
    error_reads = {v[0]: reads[v[0]] for v in grouped_reads.values() if len(v) == 1}

    corrections = {}
    for tag, read in error_reads.items():
        for i, bp in enumerate(read):
            for mut in ['A', 'C', 'G', 'T']:
                if bp == mut:
                    continue
                mut_read = list(read)
                mut_read[i] = mut
                mut_read = ''.join(mut_read)
                if any([mut_read in read_pair for read_pair in grouped_correct_reads.keys()]):
                    corrections[tag] = f'{read}->{mut_read}'
                    break

    return corrections


if __name__ == '__main__':
    reads = mb.parse_fasta('reads.fasta')
    corrections = error_corrections(reads)
    print(*corrections.values(), sep='\n')