""" Author: Yuchao He
    Aim of this file: To identify the order of alpha-carbon atoms in the chain test_q1 and data_q1
    How to run the file:  python task_1.py test_q1.txt  for test_q1.txt
                          python task_1.py data_q1.txt  for data_q1.txt
        """
import numpy as np
import sys
from collections import Counter, defaultdict
MAX_DIST, MIN_DIST = 3.9, 3.7

def read_file(text_path):
    # Read x,y,z positions of alpha-carbon atoms
    pos = []
    with open(text_path, 'r') as f:
        for line in f:
            line = line.split()
            x, y, z = map(float, line[1:])
            pos.append(np.array([x, y, z]))  
    pos = np.array(pos).astype(float)
    return pos

def closest_pairs(pos):
    
    distance_map = np.linalg.norm(
        pos[:, np.newaxis, :] - pos[np.newaxis, :, :], axis=2
    )
    dist_mask = (MIN_DIST <= distance_map) & (distance_map <= MAX_DIST)
    i_indices, j_indices = np.where(dist_mask)

    # Filter for i < j
    pairs = [(i,j) for i, j in zip(i_indices, j_indices) if i<j]

    return pairs
    

def gen_chain(pairs):

    ele_ls = [element for tup in pairs for element in tup]
    ele_counts = Counter(ele_ls)
    uni_ele = [elem for elem, count in ele_counts.items() if count == 1]
    start_atom = max(uni_ele)

    neighbors = defaultdict(list)
    for a, b in pairs:
        neighbors[a].append(b)
        neighbors[b].append(a)

    path = [start_atom]
    visited = set(path)

    while True:
        current = path[-1]
        # Find the next unvisited neighbor
        next_nodes = [node for node in neighbors[current] if node not in visited]
        if not next_nodes:  # No unvisited neighbors left
            break
        next_node = next_nodes[0]  # Select the next neighbor (arbitrary choice if multiple)
        path.append(next_node)
        visited.add(next_node)
    path = np.array(path) +1
    print("\nfound Path:")
    for ele in path:
        print(f"{ele}\n")

def main(text_path):
    pos = read_file(text_path)
    pairs = closest_pairs(pos)
    gen_chain(pairs)

if __name__ == '__main__':
    #text_path = "test_q1.txt"
    text_path = sys.argv[1]
    main(text_path)
    