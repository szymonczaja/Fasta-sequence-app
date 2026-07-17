import os 
from Bio import Entrez
from dotenv import load_dotenv
load_dotenv()

class Sequence:
    def __init__(self, sequence, header=None): 
        self.__sequence = sequence 
        self.__header = header

    @property
    def seq(self): 
        return self.__sequence 
    
    @property
    def header(self): 
        return self.__header

    def __repr__(self):
        return f"{type(self).__name__}(seq='{self.seq}', header={self.header!r})"  
      
    def __len__(self): 
        return len(self.__sequence) 
    
    def __getitem__(self, position):
        return self.__sequence[position]
    
    def __contains__(self, item):
        return item in self.__sequence
    
    def calculate_gc_pctg(self): 
        gc_count = len([x for x in self.__sequence if x == 'G' or x == 'C'])
        return f"{(gc_count / len(self.__sequence)) * 100}%" 
    
    def find_mutation(self, other): 
        if len(self.__sequence) != len(other): 
            raise ValueError('Podane sekwenje muszą być tej samej długości!')
        seq_muts = [] 
        for idx, (n1, n2) in enumerate(zip(self.__sequence, other)): 
            if n1 != n2:
                x = (idx, f"{n1}-->{n2}")
                seq_muts.append(x) 
        return seq_muts
    

class DNASequence(Sequence):
    VALID_NUCLEOTIDE = set('ATCG')
    COMPLEMENT = {
                'A':'T',
                'T':'A',
                'G':'C',
                'C':'G'
                }
    

    def __init__(self, seq, header=None):
        super().__init__(seq.upper(), header)
        self._validate()

    def _validate(self):
        if not set(self.seq).issubset(self.VALID_NUCLEOTIDE):
            raise ValueError('Sekwencja DNA zawiera nie dozwolone nukleotydy! Dozwolone: "ATCG"')
        
    def complement(self): 
        return "".join(self.COMPLEMENT[n] for n in self.seq)
    
    @classmethod
    def from_ncbi(cls, accesion_id):
        Entrez.email = os.getenv('NCBI_EMAIL')
        Entrez.api_key = os.getenv('NCBI_API_KEY')
        handle = Entrez.efetch(db='nucleotide', id=accesion_id, rettype='fasta', retmode='text')
        result = handle.read() 
        handle.close()
        if result is None: 
            raise ValueError(f'Nie znaleziono sekwencji w NCBI dla: {accesion_id}')
        header = None
        seq_lines = []
        for line in result.split('\n'):
            if line.startswith('>'):
                header = line
            else: 
                seq_lines.append(line)
        return cls(seq="".join(seq_lines), header=header)


class RNASequence(Sequence):
    VALID_NUCLEOTIDE = set('AUCG')
    CODON_TABLE = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y",
    "UAA": "*", "UAG": "*",  
    "UGU": "C", "UGC": "C",
    "UGA": "*",              
    "UGG": "W",

    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",              
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

    def __init__(self, seq, header=None):
        super().__init__(seq.upper(), header)
        self._validate()

    def _validate(self):
        if not set(self.seq).issubset(self.VALID_NUCLEOTIDE):
            raise ValueError('Sekwencja RNA zawira błędy! Dozwolone: "AUCG"')
        
    def translate(self):
        if len(self.seq) % 3 != 0:
            raise ValueError('Sekwencja RNA nie jest prawidłowa!')
        protein = [] 
        for i in range(0, len(self.seq), 3): 
            codon = self.seq[i:i+3]
            if self.CODON_TABLE[codon] == '*':
                break
            protein.append(self.CODON_TABLE[codon])
        return "".join(protein)
        
class FastaParser: 
    def __init__(self, file_path): 
        self.file_path = file_path
        self._check_file() 
        self.sequences = [] 

    def _check_file(self): 
        if not os.path.exists(self.file_path): 
            raise FileNotFoundError(f"Plik: {self.file_path} nie istnieje!")
        if not os.path.isfile(self.file_path): 
            raise ValueError(f"{self.file_path} nie jest plikiem!")
        
    def _detect_type(self, seq, header):
        seq_set = set(seq)
        if 'U' in seq_set:
            return RNASequence(seq, header=header) 
        elif 'T' in seq_set:
            return DNASequence(seq, header=header) 
        else: 
            raise ValueError('Nie da się okreslic typu!')

    def parse(self): 
        header = None 
        seq_lines = [] 
        with open(self.file_path, 'r') as f: 
            for line in f: 
                line = line.strip() 
                if not line: 
                    continue 
                if line.startswith('>'):
                    if header is not None: 
                        self.sequences.append(self._detect_type(seq="".join(seq_lines), header=header))
                    header = line
                    seq_lines = [] 
                else: 
                    seq_lines.append(line) 
        if header is not None: 
            self.sequences.append(self._detect_type(seq="".join(seq_lines), header=header))
        return self.sequences
    
    def to_fasta(self, file_path):
        with open(file_path, 'w') as f: 
            for seq_obj in self.sequences: 
                f.write(f"{seq_obj.header}\n{seq_obj.seq}\n")
        return f"Sekwencja zapisana w pliku: {file_path}"
    
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

class SequenceAnalyzer:
    def __init__(self, sequences):
        self.sequences = sequences 

    def gc_content(self): 
        result = {} 
        for obj in self.sequences: 
            gc_pctg = obj.calculate_gc_pctg()
            result[obj.header] = gc_pctg
        return result 
    
    def align(self, seq1, seq2, blosum_file_path=None):
        if set(seq1).issubset(set('AUCGT')) and set(seq2).issubset(set('AUCGT')):
            alignment = NucleotideAlignment(seq1, seq2).align()
        else:
            alignment = ProteinAlignment(seq1, seq2, blosum_file_path).align()
        return alignment

    def find_mutations(self, seq1, seq2): 
        mutations = Sequence(seq1).find_mutation(seq2)
        return mutations
    
    def translate_all(self):
        result = {} 
        for obj in self.sequences: 
            if isinstance(obj, RNASequence):
                proteins = obj.translate() 
                result[obj.header] = proteins 
        return result
