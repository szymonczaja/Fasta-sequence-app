class SequenceAlignment:
    def __init__(self, seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
        self.seq1 = seq1 
        self.seq2 = seq2 
        self.match_score = match_score 
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty 

    def _score(self, n1, n2):
        return self.match_score if n1 == n2 else self.mismatch_score

    def align(self):
        n = len(self.seq1) 
        m = len(self.seq2)
        score_matrix = [[0] * (m+1) for _ in range(n+1)]
        traceback_matrix = [[""] * (m+1) for _ in range(n + 1)]
        for i in range(1, n+1):
            score_matrix[i][0] = i * self.gap_penalty
            traceback_matrix[i][0] = "up"
        for j in range(1, m+1):
            score_matrix[0][j] = j * self.gap_penalty
            traceback_matrix[0][j] = "left"
        for i in range(1, n + 1):
            for j in range(1, m + 1): 
                diag = score_matrix[i - 1][j - 1] + self._score(self.seq1[i-1], self.seq2[j-1])
                up = score_matrix[i -1][j] + self.gap_penalty
                left = score_matrix[i][j - 1] + self.gap_penalty
                max_score = max(diag, up, left) 
                score_matrix[i][j] = max_score 
                if max_score == diag:
                    traceback_matrix[i][j] = "diag" 
                elif max_score == up: 
                    traceback_matrix[i][j] = "up"
                else: 
                    traceback_matrix[i][j] = "left"
        aligned_seq1 = ""
        aligned_seq2 = ""
        i, j = n, m 
        while i > 0 or j > 0: 
            if traceback_matrix[i][j] == "diag": 
                aligned_seq1 = self.seq1[i - 1] + aligned_seq1
                aligned_seq2 = self.seq2[j - 1] + aligned_seq2
                i -= 1 
                j -= 1 
            elif traceback_matrix[i][j] == "up": 
                aligned_seq1 = self.seq1[i - 1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                i -= 1 
            elif traceback_matrix[i][j] == "left":
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = self.seq2[j - 1] + aligned_seq2
                j -= 1 
            else:
                break
        return aligned_seq1, aligned_seq2, score_matrix[n][m]

class ProteinAlignment(SequenceAlignment): 
    def __init__(self, seq1, seq2, blosum_file_path, gap_penalty=-1):
        super().__init__(seq1, seq2, gap_penalty=gap_penalty)
        self.blosum_file_path = blosum_file_path
        self.blosum_matrix = self._load_blosum_matrix()

    def _load_blosum_matrix(self):
        try:
            with open(self.blosum_file_path, "r") as f: 
                lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                headers = lines[0].split() 
                self.blosum_matrix = {} 
                for line in lines[1:]:
                    parts = line.split() 
                    row_aa = parts[0] 
                    scores = list(map(int, parts[1:]))
                    self.blosum_matrix[row_aa] = dict(zip(headers, scores))
        except FileNotFoundError:
            print(f"Error: '{self.blosum_file_path}' not exists!")
        except ValueError:
            print(f"Error: Format is not correct!")
        return self.blosum_matrix
    
    def _score(self, p1, p2): 
        return self.blosum_matrix[p1][p2]


class NucleotideAlignment(SequenceAlignment):
    def __init__(self, seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
        super().__init__(seq1, seq2, match_score, mismatch_score, gap_penalty)
