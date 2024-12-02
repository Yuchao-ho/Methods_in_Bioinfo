def levenshtein_distance(X: str, Y: str):
        """Calculate the Levenshtein distance"""
        m, n = len(X), len(Y)
        # Initialize the matrix
        lev = [[0] * (n + 1) for _ in range(m + 1)]
        
        # Fill first row and column
        for i in range(m + 1):
            lev[i][0] = i
        for j in range(n + 1):
            lev[0][j] = j
            
        # Fill the rest of the matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if X[i-1] == Y[j-1]:
                    lev[i][j] = lev[i-1][j-1]
                else:
                    lev[i][j] = min(lev[i-1][j] + 1,    # deletion
                                 lev[i][j-1] + 1,       # insertion
                                 lev[i-1][j-1] + 1)     # substitution
        
        print(f"\n")
        print(X)
        print(Y)
        print(f"Levenshtein Distance: {lev[m][n]}\n")
        return lev[m][n]

# Test sequences
if __name__ == "__main__":    

    X = "KQTGKGS"
    Y = "KSAGKGAI"
    levenshtein_distance(X, Y)

    X = "TAAACGTCGT"
    Y = "AAACGTCGTA"
    levenshtein_distance(X, Y)
