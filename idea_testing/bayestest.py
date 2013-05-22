from Bayes.bayes import BayesCalculation
from cogent.parse.fasta import MinimalFastaParser
from sys import argv
from cogent import LoadSeqs, RNA

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
    aln = LoadSeqs(argv[1], moltype=RNA)
    headers = aln.getSeqNames()
    seqs = map(str, aln.iterSeqs())
    print headers
    print seqs
    test = BayesWrapper('test', headers, seqs)
    temperature = 37
    bayescalc = BayesCalculation(test)
    bayescalc.run()
    print str(bayescalc.Alignment.Structures).split()[1]
    #sequences = RNAAlignment(sequences=seqs)
    #structures = sequences.fold(temperature, 2, 100) 
    #alignment = RNAStructureAlignment(sequences,structures,temperature).Structures
    #print str(alignment).split("\n")[0]