from Bio import Entrez
from .base import Sequence
import os 
from dotenv import load_dotenv
load_dotenv()

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

