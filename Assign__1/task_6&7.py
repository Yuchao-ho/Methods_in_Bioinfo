import numpy as np

class SequenceAlignment:
    def __init__(self,):
        self.MATCH_SCORE = 2
        self.MISMATCH_SCORE = -1
        self.GAP_PENALTY = 2
        self.STOP = 0
        self.UP = 1
        self.LEFT = 2
        self.DIAG = 3

    def count_optimal_alignments(self, X: str, Y: str):
        """Count the total number of optimal alignments"""
        m, n = len(X), len(Y)
        # Score matrix
        F = [[0] * (n + 1) for _ in range(m + 1)]
        # Count matrix to store number of optimal alignments
        count = [[0] * (n + 1) for _ in range(m + 1)]
        
        # Initialize matrices
        F[0][0] = 0
        count[0][0] = 1
        for i in range(1, m + 1):
            F[i][0] = F[i-1][0] - self.GAP_PENALTY
            count[i][0] = 1
        for j in range(1, n + 1):
            F[0][j] = F[0][j-1] - self.GAP_PENALTY
            count[0][j] = 1
            
        # Fill matrices
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                score_match = F[i-1][j-1] + (self.MATCH_SCORE if X[i-1] == Y[j-1] else self.MISMATCH_SCORE)
                score_up = F[i-1][j] - self.GAP_PENALTY
                score_left = F[i][j-1] - self.GAP_PENALTY
                
                # Find maximum score
                F[i][j] = max(score_match, score_up, score_left)
                score_ls = np.array([score_match, score_up, score_left])
                count_ls = np.array([count[i-1][j-1], count[i-1][j], count[i][j-1]])
                mask = score_ls == score_ls.max()
                count[i][j] = np.sum(count_ls[mask])
                    
        return count[m][n]

    def get_all_optimal_alignments(self, X: str, Y: str):
        """Generate all optimal alignments"""
        def backtrack(i, j, aligned_X, aligned_Y, alignments):
            if i == 0 and j == 0:
                alignments.append((aligned_X[::-1], aligned_Y[::-1]))
                return
                
            current_score = F[i][j]
            
            # Check diagonal move
            if i > 0 and j > 0:
                score = self.MATCH_SCORE if X[i-1] == Y[j-1] else self.MISMATCH_SCORE
                if F[i-1][j-1] + score == current_score:
                    backtrack(i-1, j-1, aligned_X + X[i-1], aligned_Y + Y[j-1], alignments)
                    
            # Check up move
            if i > 0:
                if F[i-1][j] - self.GAP_PENALTY == current_score:
                    backtrack(i-1, j, aligned_X + X[i-1], aligned_Y + '-', alignments)
                    
            # Check left move
            if j > 0:
                if F[i][j-1] - self.GAP_PENALTY == current_score:
                    backtrack(i, j-1, aligned_X + '-', aligned_Y + Y[j-1], alignments)

        m, n = len(X), len(Y)
        # Score matrix
        global F
        F = [[0] * (n + 1) for _ in range(m + 1)]
        
        # Initialize matrix
        for i in range(m + 1):
            F[i][0] = -i * self.GAP_PENALTY
        for j in range(n + 1):
            F[0][j] = -j * self.GAP_PENALTY
            
        # Fill matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = F[i-1][j-1] + (self.MATCH_SCORE if X[i-1] == Y[j-1] else self.MISMATCH_SCORE)
                delete = F[i-1][j] - self.GAP_PENALTY
                insert = F[i][j-1] - self.GAP_PENALTY
                F[i][j] = max(match, delete, insert)
        
        alignments = []
        backtrack(m, n, "", "", alignments)
        return alignments

    def output(self, X:str, Y:str):
        print(f"\n")
        print(X)
        print(Y)
        count = self.count_optimal_alignments(X, Y)
        print(f"\n1. Number of optimal alignments: {count}")
    
        if count < 10:
            alignments = self.get_all_optimal_alignments(X, Y)
            print("\n2. All optimal alignments:")
            for i, (align_x, align_y) in enumerate(alignments, 1):
                print(f"\nAlignment {i}:")
                print(align_x)
                print(align_y)


if __name__ == "__main__":
    seq = SequenceAlignment()
    X = "ATTA"
    Y = "ATTTTA"    
    seq.output(X, Y)
    
    X = "AAAAA"
    Y = "AAAABAAAA"
    seq.output(X, Y)    