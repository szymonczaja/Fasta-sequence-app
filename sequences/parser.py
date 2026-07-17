from .rna import RNASequence
from .dna import DNASequence
import os 

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