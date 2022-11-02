def edges_for_tree(n: int, adjacency_list: list):
    # Min. required edges = number of trees in adjacency list - 1 + no. of nodes not in adjacency list

    # Counting nodes not in existing trees
    nodes = set(sum(adjacency_list, []))
    lonely_nodes = 0
    for k in [m+1 for m in range(n)]:
        if k not in [*nodes]:
            lonely_nodes += 1

    # Sorting to make sure that assembly of trees occurs correctly
    adjacency_list = sorted([sorted(pair) for pair in adjacency_list])

    # Assembling trees
    adjacency_list_trees = [adjacency_list.pop(0)]
    for pair in adjacency_list:
        appended = False
        for tree in adjacency_list_trees:
            if pair[0] in tree:
                tree.append(pair[1])
                appended = True
                break
            elif pair[1] in tree:
                tree.append(pair[0])
                appended = True
                break
        if not appended:
            adjacency_list_trees.append(pair)
    lacking_branches = len(adjacency_list_trees) - 1

    min_edges = lonely_nodes + lacking_branches

    return min_edges


if __name__ == '__main__':
    with open('input.txt', 'r') as file:
        lines = file.readlines()
    n = int(lines.pop(0))
    adjacency_list = [line.strip().split() for line in lines]
    adjacency_list = [[int(i) for i in pair] for pair in adjacency_list]
    print(edges_for_tree(n, adjacency_list))