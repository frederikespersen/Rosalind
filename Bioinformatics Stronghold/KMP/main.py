import sys
sys.path.append('../../MolBio')
import molbio as mb


def failure_array(sequence: str) -> list[int]:
    # Initialising failure array
    P = [0] * len(sequence)

    # Looping over sequence elements
    for k in range(len(sequence)):

        print(f"\r[{k}/{len(sequence)-1}]", end='', flush=True)

        # P[0] = 0 by definition
        if k == 0:
            continue

        # If the previous k was successful...
        if P[k-1] > 0:

            # Check if the next character of the prefix is this character
            if sequence[k] == sequence[P[k-1]+1-1]:
                P[k] = P[k-1] + 1
                continue

            # Else start a new search; Note that the subsequence can be no longer than the previous; fewer j to search
            # Shorter-assumption: k - j < P[k-1] --> k - P[k-1] < j
            else:
                # Looping over increasing j; decreasing length
                for j in filter(lambda j: k - P[k-1] < j, [*range(k)]):

                    if sequence[k-j-1] != sequence[k]:
                        continue

                    prefix = sequence[:k - j]
                    subsequence = sequence[j + 1:k + 1]
                    assert len(prefix) == len(subsequence)

                    # If match, longest sequence found. Note length and break loop
                    if prefix == subsequence:
                        P[k] = (k - j)
                        break

        # Else, start over; If previous was 0, then this can have length of max 1
        else:
            if sequence[k] == sequence[0]:
                P[k] = 1

    return P


if False:
    with open('output.txt', 'r') as file:
        lines = file.readlines()
    print(max(lines[0]))
elif __name__ == '__main__':
    sequences = mb.parse_fasta('sequence.fasta')
    sequences = [*sequences.values()]
    P = failure_array(sequences[1])
    with open('output.txt', 'w') as file:
        file.write(' '.join([str(p) for p in P]))
        print()
    print(*P)
