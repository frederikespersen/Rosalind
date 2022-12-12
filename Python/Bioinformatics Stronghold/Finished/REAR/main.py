# Solution inspired by
# https://medium.com/@matthewwestmk/calculating-reversal-distance-using-parks-exact-greedy-algorithm-87c62d690eef

def reversal(sequence: list, start_index: int, end_index: int) -> list:
    prefix = sequence[:start_index]
    reversed = sequence[start_index:end_index][::-1]
    suffix = sequence[end_index:]
    return prefix + reversed + suffix


def breakpoints(sequence: list, target_sequence: list) -> list:
    # Find all breakpoints;
    # Pick a pair of adjacent elements in the sequence, find their corresponding indices in the target sequence.
    # Their distance will be either 1 or -1, if the elements are supposed to be adjacent.
    # If not, there has occured a break - we refer to the index of the second element as the breakpoint
    breakpoints = []

    # Check all adjacent element pairs in sequence
    for i in range(len(sequence)-1):
        current_element = sequence[i]
        adjacent_element = sequence[i+1]

        # Find the difference between the corresponding index of the elements in the target sequence
        target_index_difference = target_sequence.index(current_element) - target_sequence.index(adjacent_element)

        # If the difference is not 1, then a reversal has occcured; Reversals are confined within two breakpoints
        if abs(target_index_difference) != 1:
            breakpoints.append(i+1)

    return breakpoints


def generate_reversals(sequences: list, target_sequence: list) -> list:
    # Generate all logical (defined by two breakpoints) 1-reversal variations for each given sequence
    reversals = []
    for sequence in sequences:

        # Find breakpoints of sequence
        sequence_breakpoints = breakpoints(sequence, target_sequence)

        # Generate reversals with pairs of breakpoints; We don't know how the breakpoints pair up
        break_start_max_index = -1
        for break_start in sequence_breakpoints[:break_start_max_index]:
            break_end_min_index = sequence_breakpoints.index(break_start)+1
            for break_end in sequence_breakpoints[break_end_min_index:]:
                reversals.append(reversal(sequence, break_start, break_end))

    return reversals


def generate_best_reversals(sequences: list, target_sequence: list) -> list:
    # Generate all 1-reversal variations of sequences
    reversals = generate_reversals(sequences, target_sequence)

    # Setting variable minimum number of breakpoints; Starting with max
    min_no_breakpoints = len(target_sequence)

    # Finding reversals with the largest reduction in breakpoint number - i.e. closer to a true reversal
    best_reversals = []
    for reversal in reversals:

        # Checking number of breakpoints after reversal
        no_breakpoints = len(breakpoints(reversal, target_sequence))

        # If this sequence is better than former best sequences, replace the best sequences
        if no_breakpoints < min_no_breakpoints:
            best_reversals = [reversal]
            # Reset best breakpoint number
            min_no_breakpoints = no_breakpoints

        # If this sequence is as good as former best sequences, add it to list
        elif no_breakpoints == min_no_breakpoints:
            best_reversals.append(reversal)

    return best_reversals


def reversal_distance(permutations: tuple) -> int:
    # Setting arbitrary sequence and target sequence
    sequence = permutations[0]
    target_sequence = permutations[1]

    # Adding polarity to sequences with padding; otherwise infinite while loop
    sequence = ['-'] + sequence + ['+']
    target_sequence = ['-'] + target_sequence + ['+']

    reversal_sequences = [sequence]
    reversals = 0
    # Creating the best reversals until we create the target sequence
    while target_sequence not in reversal_sequences:
        reversal_sequences = generate_best_reversals(reversal_sequences, target_sequence)
        # Count the number of reversals
        reversals += 1

    return reversals



if __name__ == '__main__':
    # Parsing input
    with open("permutations.txt", 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    pairs = [lines[i*3:i*3+2] for i in range(len(lines)//3 + 1)]
    pairs = [tuple(s.split(' ') for s in pair) for pair in pairs]
    reversal_distances = [reversal_distance(pair) for pair in pairs]
    print(*reversal_distances)
