""" Author: Yuchao He 
    Aim of this file: To calculate distances between consecutive alpha-carbon atoms
    How to run the file:  python cal_dist.py
        """
import numpy as np

def read_pdb(pdb_file):
    # Reads a PDB file to extract index and coordinates of C-alpha atoms.
    ca_atoms = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                index = int(line[6:11].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ca_atoms.append((index, np.array([x, y, z])))
    return ca_atoms

def compute_consecutive_distances(ca_atoms):
    # Computes and analyse distances between consecutive C-alpha atoms.
    distances = []
    for i in range(len(ca_atoms) - 1):
        dist = np.linalg.norm(ca_atoms[i][1] - ca_atoms[i + 1][1])
        distances.append(dist)
    aim_dist= 3.8
    mean_dist= np.mean(distances)
    std_dist = np.std(distances)
    deviations = [dist - aim_dist for dist in distances]
    min_dev, max_dev = min(deviations), max(deviations)

    return mean_dist, std_dist, min_dev, max_dev

def main(pdb_pth):
    #Print out final results.
    ca_atoms = read_pdb(pdb_pth)
    mean_dist, std_dist, min_dev, max_dev = compute_consecutive_distances(ca_atoms)
    print(f"Mean distance: {mean_dist:.3f} Å")
    print(f"Standard deviation: {std_dist:.3f} Å")
    print(f"Min deviation from 3.8 Å: {min_dev:.3f} Å")
    print(f"Max deviation from 3.8 Å: {max_dev:.3f} Å")


if __name__ == "__main__":
    pdb_file = "1crn.pdb"  
    main(pdb_file)
    

    