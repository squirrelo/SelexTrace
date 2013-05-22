from cogent.app.muscle_v38 import align_unaligned_seqs
from Bayes.bayes import BayesCalculation
from cogent import RNA


class BayesInputWrapper:
    def __init__(self, name, labels, seqs, temp='37'):
        self.name = name
        self.temperature = temp
        self.primer3 = ''
        self.primer5 = ''
        self.labels = labels
        self.sequences = seqs
        self.mappings = []
        self.id = '12345'


def bayesfold(seqs, temperature=37):
    '''Takes in a list of unaligned tuples [(header, seq)] and returns
    most likely structure from bayesfold'''
    aln = align_unaligned_seqs(seqs, RNA)
    test = BayesInputWrapper('test', aln.getSeqNames(),
        map(str, aln.iterSeqs()), str(temperature))
    bayescalc = BayesCalculation(test)
    bayescalc.run()
    return str(bayescalc.Alignment.Structures).split()[1]
