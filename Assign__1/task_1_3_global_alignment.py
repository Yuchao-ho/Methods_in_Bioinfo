""" Modify the program global_alignment.c (or equivalent) so that an extra line of output is printed between 
the two aligned sequence, indicating exact matches with the character "|", e.g.

    AT-CGAT
    || || |
    ATACG-T """

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_PENALTY = 2
STOP, UP, LEFT, DIAG = 0, 1, 2, 3

def global_alignment(X, Y):
    m, n = len(X), len(Y)
    F = [[0] * (n + 1) for _ in range(m + 1)]  # Score matrix
    trace = [[0] * (n + 1) for _ in range(m + 1)]  # Trace matrix

    # Initialize matrices
    for i in range(1, m + 1):
        F[i][0] = F[i - 1][0] - GAP_PENALTY
        trace[i][0] = UP
    for j in range(1, n + 1):
        F[0][j] = F[0][j - 1] - GAP_PENALTY
        trace[0][j] = LEFT

    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if X[i - 1] == Y[j - 1]:
                score = F[i - 1][j - 1] + MATCH_SCORE
            else:
                score = F[i - 1][j - 1] + MISMATCH_SCORE
            trace[i][j] = DIAG

            tmp = F[i - 1][j] - GAP_PENALTY
            if tmp > score:
                score = tmp
                trace[i][j] = UP

            tmp = F[i][j - 1] - GAP_PENALTY
            if tmp > score:
                score = tmp
                trace[i][j] = LEFT

            F[i][j] = score

    # Traceback
    alignX, alignY = [], []
    i, j = m, n

    while i > 0 or j > 0:
        if trace[i][j] == DIAG:
            alignX.append(X[i - 1])
            alignY.append(Y[j - 1])
            i -= 1
            j -= 1
        elif trace[i][j] == UP:
            alignX.append(X[i - 1])
            alignY.append('-')
            i -= 1
        elif trace[i][j] == LEFT:
            alignX.append('-')
            alignY.append(Y[j - 1])
            j -= 1

    # Reverse alignments
    alignX = alignX[::-1]
    alignY = alignY[::-1]

    match_indicator = []
    match_count = 0
    mismatch_count = 0
    non_gap = len(alignX)

    for x_char, y_char in zip(alignX, alignY):
        if x_char == y_char:
            match_indicator.append('|')
            match_count += 1
        else:
            match_indicator.append(' ')
            mismatch_count += 1
            if x_char == "-" or y_char == "-":
                non_gap -= 1

    # Calculate Percent Identity: 
    # Divided by the length of alignments
    alignment_length = len(alignX)
    percent_identity = (match_count / alignment_length) * 100
    # Divided by the number of non-gap positions:
    percent_identity_gap = (match_count / non_gap) * 100
    # Hamming Distance (count of mismatches in the aligned sequences)
    if m == n:
        hamming_distance = mismatch_count
    else:
        hamming_distance = None

    # Print alignment with match indicator
    print(f"Global Alignment: ")
    print("".join(alignX))
    print("".join(match_indicator))
    print("".join(alignY))
    print(f"Percent Identity: \n")
    print(f"  - Divided by the length of the alignment ({alignment_length}): {percent_identity:.2f}%")
    print(f"  - Divided by the number of non-gap positions ({non_gap}): {percent_identity_gap:.2f}%")
    print(f"Hamming Distance: {hamming_distance}\n")

# Test sequences
if __name__ == "__main__":

    X = "ATCGAT"
    Y = "ATACGT"
    global_alignment(X, Y)

    X = "ACGATAGCGAAACCAAAA"
    Y = "CACGTAGCCGATGTC"
    global_alignment(X, Y)

    X = "ATTA"
    Y = "ATTTTA"
    global_alignment(X, Y)
