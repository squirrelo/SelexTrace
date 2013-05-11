from alignment import RNAAlignment, RNAStructureAlignment

class BayesWrapper:
    def __init__(self,name,labels,seqs):
        self.name = name
        self.temperature = '37'
        self.primer3 = ''
        self.primer5 = ''
        self.labels = labels
        self.sequences = seqs
        self.mappings = []
        self.id = '12345'


if __name__ == "__main__":
    test = BayesWrapper('test', ['a','b'], ['GUCAAAAGAC','GUCAAAAGAC'])
    temperature = 37
    sequences = RNAAlignment(
            sequences=['GUCAAAAGAC','GUCAAAAGAC'])
    structures = sequences.fold(temperature, 2,
                100) 
    alignment = RNAStructureAlignment(sequences,structures,temperature).Structures
    print str(alignment).split("\n")[0]