import sys
sys.path.append('../../../MolBio')
import molbio as mb

def failure_array(sequence: str) -> list[int]:
    # Inspired by:
    # http://saradoesbioinformatics.blogspot.com/2016/08/speeding-up-motif-finding.html

    P = [0] * len(sequence)

    # l is the length of the longest substring;
    # Also the index of the last element of the prefix
    l = 0

    # Looping over elements i (except first) over sequence
    for i in range(1, len(sequence)):

        # Using l from previous loop

        # While the longest length is not 0 and
        # while the next element of the previous longest prefix (index (l-1)+1=l) does not match the current character
        # then step back l one time. ?The reason why this works is unclear to me?.
        # In short, decrease l until a sequence match is found, or l is 0
        while l > 0 and sequence[l] != sequence[i]:
            l = P[l - 1]

        # If the next element of the previous longest prefix (index (l-1)+1=l) matches the current character
        # then increase the longest length
        if sequence[l] == sequence[i]:
            l += 1

        # Set P[i] equal to current l
        P[i] = l

    return P


if __name__ == '__main__':
    sequences = mb.parse_fasta('sequence.fasta')
    sequences = [*sequences.values()]
    P = failure_array(sequences[0])
    P = map(str, P)
    P = ' '.join(P)
    with open('output.txt', 'w') as file:
        file.write(P)
        print()
    print(P)

# print(max(open('output.txt', 'r').readlines()[0][0]))