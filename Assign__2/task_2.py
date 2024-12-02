""" Author: Yuchao He
    Aim of this file: To identify the order of alpha-carbon atoms in the main chain from data_q2.txt
    How to run the file:  python task_2.py 
        """
import numpy as np
import copy
from collections import defaultdict
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

def find_atom_relate(pos):
    # Calculate distance matrix
    distance_map = np.linalg.norm(
        pos[:, np.newaxis, :] - pos[np.newaxis, :, :], axis=2
    )
    dist_mask = (MIN_DIST <= distance_map) & (distance_map <= MAX_DIST)
    i_indices, j_indices = np.where(dist_mask)

    # Find pairs and neighbors
    pairs = [(i, j) for i, j in zip(i_indices, j_indices) if i < j]
    neighbors = defaultdict(list)
    for a, b in pairs:
        neighbors[a].append(b)
        neighbors[b].append(a)

    # Separate single and multi-neighbor atoms
    single_neighbor_keys = [key for key, value in neighbors.items() if len(value) == 1]
    multi_neighbor_keys = [key for key, value in neighbors.items() if len(value) > 0]
    single_neighbor_atoms = {key: neighbors[key] for key in single_neighbor_keys}  # only one neighbor
    multi_neighbor_atoms = {key: neighbors[key] for key in multi_neighbor_keys}       # more than one neighbor

    return single_neighbor_atoms, multi_neighbor_atoms


def calculate_mean_distance(pos, multi_neighbor_atoms):
    # Calculate mean distance of all distances between neighbor atoms and selected atom
    mean_distances = {}
    for atom, neighbors in multi_neighbor_atoms.items():
        mean_pos = np.mean([pos[neighbor] for neighbor in neighbors], axis=0)
        distances = [np.linalg.norm(pos[neighbor] - mean_pos) for neighbor in neighbors]
        mean_distances[atom] = distances
    return mean_distances

def find_most_related_atoms(mean_distances, multi_neighbor_atoms):
    # Find the most related atoms based on center distances
    most_related_atoms = {}
    for atom, distances in mean_distances.items():
        min_distance = np.min(distances)
        most_related_atoms[atom] = [
            multi_neighbor_atoms[atom][i] 
            for i, dist in enumerate(distances) 
            if abs(dist - min_distance) < 1.2
        ]
    return most_related_atoms

def find_paths(start_atom, most_related_atoms, chain, chains):
    # Recursive function to find all paths
    chain = copy.deepcopy(chain)
    chain.append(start_atom)

    # Remove atoms already in the chain
    related_atoms = most_related_atoms[start_atom]

    for atom in related_atoms:
        if atom in chain:
            related_atoms.remove(atom)
        
        if atom not in list(most_related_atoms.keys()):
            related_atoms.remove(atom)

    if len(related_atoms) == 0:
        chains.append(chain)
        return

    for atom in related_atoms:
        find_paths(atom, most_related_atoms, chain, chains)


def main(text_path):
    pos = read_file(text_path)
    single_neighbor_atoms, multi_neighbor_atoms = find_atom_relate(pos)
    mean_distances = calculate_mean_distance(pos, multi_neighbor_atoms)
    most_related_atoms = find_most_related_atoms(mean_distances, multi_neighbor_atoms)

    # Find all chains starting from single neighbor atoms
    chains = []
    for start_atom in single_neighbor_atoms.keys():
        find_paths(start_atom, most_related_atoms, [], chains)

    # Find the longest chain
    longest_chain = max(chains, key=len)

    for atom in longest_chain:
        print(atom+1)
    print(f"\nThe total number of alpha-carbon atoms in the main chain: {len(longest_chain)}")
    
if __name__ == '__main__':
    text_path = "data_q2.txt"
    main(text_path)