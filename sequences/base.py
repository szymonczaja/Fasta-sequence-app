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