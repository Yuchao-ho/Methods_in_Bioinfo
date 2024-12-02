import numpy as np
from task_1_3_global_alignment import global_alignment

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_PENALTY = -2

def smith_waterman(X, Y):
    m, n = len(X), len(Y)
    
    # Initialize score and traceback matrices
    F = np.zeros((m + 1, n + 1), dtype=int)
    trace = np.zeros((m + 1, n + 1), dtype=int)

    # Variables to store the highest scoring cell for traceback
    max_score = 0
    max_i, max_j = 0, 0
    
    # Fill the matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = F[i - 1, j - 1] + (MATCH_SCORE if X[i - 1] == Y[j - 1] else MISMATCH_SCORE)
            delete = F[i - 1, j] + GAP_PENALTY
            insert = F[i, j - 1] + GAP_PENALTY
            F[i, j] = max(0, match, delete, insert)
            
            # Traceback information for alignment
            if F[i, j] == match:
                trace[i, j] = 3  # Diagonal
            elif F[i, j] == delete:
                trace[i, j] = 1  # Up
            elif F[i, j] == insert:
                trace[i, j] = 2  # Left

            # Track the cell with the highest score
            if F[i, j] > max_score:
                max_score = F[i, j]
                max_i, max_j = i, j

    alignX, alignY = [], []
    i, j = max_i, max_j

    while F[i, j] != 0:
        if trace[i, j] == 3:  # Diagonal
            alignX.append(X[i - 1])
            alignY.append(Y[j - 1])
            i -= 1
            j -= 1
        elif trace[i, j] == 1:  # Up
            alignX.append(X[i - 1])
            alignY.append('-')
            i -= 1
        elif trace[i, j] == 2:  # Left
            alignX.append('-')
            alignY.append(Y[j - 1])
            j -= 1

    # Reverse the alignments since they were built backwards
    alignX = alignX[::-1]
    alignY = alignY[::-1]
    match_line = ''.join('|' if alignX[k] == alignY[k] else ' ' for k in range(len(alignX)))
    
    # Display the result
    global_alignment(X,Y)
    print("\nMax alignment score:", max_score)
    print("Optimal Local Alignment:")
    print("".join(alignX))
    print(match_line)
    print("".join(alignY))

# Test sequences
if __name__ == "__main__":    

    X = "KQTGKGS"
    Y = "KSAGKGAI"
    smith_waterman(X, Y)
