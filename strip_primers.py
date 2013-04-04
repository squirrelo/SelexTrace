from qiime.split_libraries import local_align_primer_seq
from cogent.parse.fasta import MinimalFastaParser

def strip_primer(seqfile, primer, maxmismatch = 0):
    '''strips 3 prime primer from sequences in fasta file and returns MinimalFastaParser
    formatted arrays for stripped and not stripped sequences'''
    seqs = MinimalFastaParser(open(seqfile, 'rU'))
    nostripped = []
    stripped = []
    for seq in seqs:
            #code adapted from truncate_reverse_primers.py in qiime
            rev_primer_mm, rev_primer_index =\
             local_align_primer_seq(primer, seq[1])
            if rev_primer_mm > maxmismatch:
                nostripped.append(seq)
            else:
                seqnew = seq[1][:rev_primer_index]
                stripped.append((seq[0], seqnew))
    #end for
    return stripped, nostripped

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Strips primers from sequences in a passed FASTA file.")
    parser.add_argument('-i', required=True, help="FASTA file input for cleaning")
    parser.add_argument('-o', required=True, help="output path for files")
    parser.add_argument('-m', default=2, help="max number of mismatches with primer allowed (default: 2)")
    parser.add_argument('-p', default="", required=True, help="3' primer sequence to strip")

    args = parser.parse_args()

    filein = args.i
    filenamebase = filein[filein.rfind('/') + 1: -6]
    if args.o[:-1] != '/':
        args.o += "/"
    maxmismatch = int(args.m)

    print "===================="
    print "Filein: " + args.i
    print "Save directory: " + args.o
    print "3' primer: " + args.p
    print "Max mismatches: " + str(maxmismatch)
    print ""
    kept, removed = strip_primer(args.i, args.p, args.m)
    fileout = open(args.o + filenamebase + "-Stripped.fasta", 'w')
    for seq in kept:
        fileout.write('>%s\n%s\n' % (seq[0], seq[1]))
    fileout.close()
    leftout = open(args.o + filenamebase + "-NotStripped.fasta", 'w')
    for seq in removed:
        leftout.write('>%s\n%s\n' % (seq[0], seq[1]))
    leftout.close()

    print str(kept) + " sequences stripped of primers"
    print str(left) + " sequences without primers"
