from cogent.struct.rna2d import ViennaStructure
from cogent.parse.fasta import MinimalFastaParser
from subprocess import Popen, PIPE
from StringIO import StringIO
import pylab as p
from sys import argv

if __name__ == "__main__":
    seqslen = []
    fin = open(argv[1], 'U')
    for head, seq in MinimalFastaParser(fin):
        seqslen.append(len(seq))
    seqsfin = []
    print min(seqslen)-1, max(seqslen)+1
    for length in range(min(seqslen)-1, max(seqslen)+1):
        count = sum(1 for x in seqslen if x == length)
        seqsfin.append(count)
fig = p.figure()
ax = fig.add_subplot(1,1,1)
ax.bar(range(min(seqslen)-1, max(seqslen)+1), seqsfin)
p.show()

