
from strip_primers import strip_primer
from os.path import exists
from os import mkdir
from sys import exit
from cogent.parse.fasta import MinimalFastaParser
from cogent import LoadSeqs
from remove_duplicates import remove_duplicates
from tempfile import mktemp
from time import clock
import argparse


def write_fasta_list(lst, filename):
    '''writes MinimalFastaParser formatted list [(header,sequence)] to fasta file filename'''
    fileout = open(filename, 'w')
    for seq in lst:
        fileout.write('>%s\n%s\n' % (seq[0], seq[1]))
    fileout.close()

def write_fasta_dict(dct, filename):
    '''writes MinimalFastaParser formatted dict {header: sequence} to fasta file filename'''
    fileout = open(filename, 'w')
    for header, seq in dct:
        fileout.write('>%s\n%s\n' % (header, seq))
    fileout.close()


def rem_N_short(seqs, minlen=0):
    '''Takes in a MinimalFastaParser formatted list of sequences and returns list 
    with sequences containing Ns or shorter than minlen removed'''
    rem = []
    for i, seq in enumerate(seqs):
        if "N" in seq[1].upper() or len(seq[1]) < minlen:
            rem.append(i)
        rem.sort(reverse=True)
    for i in rem:
        seqs.pop(i)
    return seqs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cleans sequences by removing 3' \
    primers, duplicate sequences, and sequences with ambiguous bases.")
    parser.add_argument('-i', required=True, help="FASTA input file")
    parser.add_argument('-o', required=True, help="Output folder")
    parser.add_argument('-p', default="", required=True, help="3' primer \
    sequence to strip")
    parser.add_argument('-l', default = 0, type=int, help="minimum length of \
    sequence required to not be removed from pool (default 0)")
    parser.add_argument('-q', type=store_true, help="Input is fastq format \
    (default fasta)")

    args = parser.parse_args()
    if args.l < 0:
        print "ERROR: min sequence length must be greater than 1!"
        exit(1)

    #add trailing slash if necessary to folderin
    folderout = args.o
    if folderout[:-1] != '/':
        folderout += "/"

    filein = args.i
    basename = filein
    if "/" in filein:
        basename = args.i[args.i.rfind("/"):]
    basename = basename[:basename.rfind(".")]

    #make directory to store cleaned sequences in
    if not exists(folderout):
        mkdir(folderout)

    currfolder = folderout + basename

    #convert fastq to fasta if needed
    if args.q:
        seqs = LoadSeqs(args.i, format='fastq', aligned=False)
        f = open(currfolder + ".fasta", 'w')
        f.write(seqs.toFasta())
        f.close()
        filein = currfolder + ".fasta"


    print "===================="
    print "File in: " + args.i
    print "Output Folder: " + args.o
    print "3' primer: " + args.p
    print "Min length: " + str(args.l)
    print "====================\n"
    print "==Cleaning input sequences=="

    log = open(currfolder + "-cleanup.log", 'w')
    log.write("====================\nFile in: " + args.i + "\nOutput Folder: " + args.o + \
    "\n3' primer: " + args.p + "\nMin length: " + str(args.l) + "\n====================\n")
    #strip primers from sequences, print out not stripped to be safe
    #allowing up to 2 mismatches in the primer
    print "Primer stripping"
    secs = clock()
    kept, rem = strip_primer(filein, args.p, 2)
    log.write("Primer stripping\n" + str(len(kept)) + " sequences left, " + \
    str(len(rem)) + " sequences removed")
    print str(len(kept)) + " sequences left, " + \
    str(len(rem)) + " sequences removed. " + str((clock() - secs)/60) + " minutes\n"
    write_fasta_list(rem, currfolder + "-NotStripped.fasta")
    rem = 0
    #remove all sequences with Ns and short sequences
    print "Remove short and ambiguous sequences"
    secs = clock()
    kept = rem_N_short(kept, args.l)
    log.write("Remove short and ambiguous sequences\n" + str(len(kept)) + " sequences left\n")
    print str(len(kept)) + " sequences left. " + str((clock() - secs)/60) + " minutes"
    write_fasta_list(kept, currfolder + "-CleanStripped.fasta")

    #remove duplicate sequences from the fasta file and store for later
    print "Remove duplicates"
    secs = clock()
    kept, headers = remove_duplicates(kept)
    write_fasta_list(kept, currfolder + "-Unique.fasta")
    #write out file holding headers keyed to a sequence
    keyfile = open(currfolder + "-seqtoheaders.txt", 'w')
    for key in headers:
        keyfile.write(key + "\t")
        for item in headers[key]:
            keyfile.write(item + ",")
        keyfile.write("\n")
    keyfile.close()
    log.write("Remove duplicates\n" + str(len(kept)) + " sequences left")
    print str(len(kept)) + " sequences left. " + str((clock() - secs)/60) + " minutes\n"
    log.close()
