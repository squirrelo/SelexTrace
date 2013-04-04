from cogent.parse.fasta import MinimalFastaParser
from sys import argv

def write_fasta_list(lst, filename):
    '''writes MinimalFastaParser formatted list [(header,sequence)] to fasta file filename'''
    fileout = open(filename, 'w')
    for seq in lst:
        fileout.write('>%s\n%s\n' % (seq[0], seq[1]))
    fileout.close()

if __name__ == "__main__":
    seqlen = int(argv[3])
    seqs = [(header, seq) for header, seq in MinimalFastaParser(open(argv[1]))]
    newseqs = []
    for seq in seqs:
        newseqs.append((seq[0],seq[1][:seqlen]))
    write_fasta_list(newseqs, argv[2])
        