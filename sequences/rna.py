from .base import Sequence

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