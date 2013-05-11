from cogent.parse.fasta import MinimalFastaParser
from sys import argv

if __name__ == "__main__":
    filein = open(argv[1])
    fileout = open("min2.fasta", 'w')
    for header, seq in MinimalFastaParser(filein):
        if int(header.split("_")[1]) > 1:
            fileout.write(">%s\n%s\n" % (header, seq))