from .alignment import NucleotideAlignment, ProteinAlignment
from .base import Sequence
from .rna import RNASequence


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
