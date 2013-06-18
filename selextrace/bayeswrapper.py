from cogent.app.muscle_v38 import align_unaligned_seqs
from Bayes.bayes import BayesCalculation
from cogent import LoadSeqs, RNA


class BayesInputWrapper:
    def __init__(self, labels, seqs, temp='37'):
        self.name = 'name'
        self.temperature = temp
        self.primer3 = ''
        self.primer5 = ''
        self.labels = labels
        self.sequences = seqs
        self.mappings = []
        self.id = '12345'


def bayesfold(seqsin, temperature=37, params=None):
    '''Takes in LoadSeqs readable sequence data and returns
    most likely structure from bayesfold'''
    seqs = LoadSeqs(data=seqsin, moltype=RNA, aligned=False)
    if params == None:
        params = {}
    aln = align_unaligned_seqs(seqs, RNA, params=params)
    bayesinput = BayesInputWrapper(aln.getSeqNames(),
        map(str, aln.iterSeqs()), str(temperature))
    bayescalc = BayesCalculation(bayesinput)
    bayescalc.run()
    return aln, str(bayescalc.Alignment.Structures).split()[1]
