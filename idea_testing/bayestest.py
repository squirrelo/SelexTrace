from Bayes.alignment import RNAAlignment, RNAStructureAlignment
from sys import argv
from cogent.parse.fasta import MinimalFastaParser

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
    headers = []
    seqs = []
    for header, seq in MinimalFastaParser(open(argv[1])):
        if header == "SS_struct":
            continue
        headers.append(header)
        seqs.append(seq)
    test = BayesWrapper('test', headers, seqs)
    temperature = 37
    sequences = RNAAlignment(sequences=seqs)
    structures = sequences.fold(temperature, 2,
                100) 
    alignment = RNAStructureAlignment(sequences,structures,temperature).Structures
    print str(alignment)  # .split("\n")[0]