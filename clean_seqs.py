from strip_primers import strip_primer
from os.path import exists
from os import mkdir, walk
from sys import exit
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from cogent import LoadSeqs
from remove_duplicates import remove_duplicates
from tempfile import mktemp
from time import time
import argparse
from stutils import write_fasta_list, write_fasta_dict


def rem_N_short(seqs, minlen=1):
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
    parser.add_argument('-i', required=True, help="Input folder")
    parser.add_argument('-p', default="", required=True, help="3' primer \
    sequence to strip")
    parser.add_argument('-o', default = "", help="Output folder (default same as input)")
    parser.add_argument('-l', default = 1, type=int, help="minimum length of \
    sequence to keep (default 1)")
    parser.add_argument('-d', default = 2, type=int, help="minimum number of duplicates\
    needed to keep (default 2)")
    parser.add_argument('-q', action='store_true', default = False, help="Input is fastq format \
    (default fasta)")

    args = parser.parse_args()
    if args.l < 1:
        print "ERROR: min sequence length must be greater than 1!"
        exit(1)
    if args.d < 1:
        print "ERROR: min number of duplicates must be greater than 1!"
        exit(1)

    #add trailing slash if necessary to folderin
    folderin = args.i
    if folderin[:-1] != '/':
        folderin += "/"

    if args.o == "":
        folderout = folderin
    else:
        folderout = args.o
        if folderout[:-1] != '/':
            folderout += "/"
    print "===================="
    print "Folder in: " + folderin
    print "Output Folder: " + folderout
    print "3' primer: " + args.p
    print "Min length: " + str(args.l)
    print "====================\n"
    for filein in walk(args.i).next()[2]:
        #skip if not fastq or fasta file
        if filein.split(".")[-1] != "fastq" and filein.split(".")[-1] != "fasta" and filein.split(".")[-1] != "fas":
            continue
        basename = filein
        basename = basename[:basename.rfind(".")]

        #make directory to store cleaned sequences in
        if not exists(folderout + basename):
            mkdir(folderout + basename)

        currfolder = folderout + basename + "/" + basename
        print "==ROUND " + basename + "=="
        #convert fastq to fasta if needed
        if args.q:
            print "==Converting to FASTA=="
            f = open(folderout + basename + ".fasta", 'w')
            for header, seq, qual in MinimalFastqParser(folderin+filein, strict=False):
                f.write(''.join([">", header, '\n', seq, '\n']))
            f.close()
            filein = folderout + basename + ".fasta"
            seqs = 0

        print "==Cleaning input sequences=="

        log = open(currfolder + "-cleanup.log", 'w')
        log.write("====================\nFile in: " + folderin + filein + "\nOutput Folder: " + currfolder + \
        "\n3' primer: " + args.p + "\nMin length: " + str(args.l) + "\nMin duplicates: " + str(args.d) + "\n====================\n")
        #strip primers from sequences, print out not stripped to be safe
        #allowing up to 2 mismatches in the primer
        print "Primer stripping"
        secs = time()
        kept, rem = strip_primer(filein, args.p, 2)
        log.write("Primer stripping\n" + str(len(kept)) + " sequences left, " + \
        str(len(rem)) + " sequences removed")
        print str(len(kept)) + " sequences left, " + \
        str(len(rem)) + " sequences removed. " + str((time() - secs)/60) + " minutes\n"
        write_fasta_list(kept, currfolder + "-Stripped.fasta")
        write_fasta_list(rem, currfolder + "-NotStripped.fasta")
        rem = 0
        #remove all sequences with Ns and short sequences
        print "Remove short and ambiguous sequences"
        secs = time()
        kept = rem_N_short(kept, args.l)
        log.write("Remove short and ambiguous sequences\n" + str(len(kept)) + " sequences left\n")
        print str(len(kept)) + " sequences left. " + str((time() - secs)/60) + " minutes"
        write_fasta_list(kept, currfolder + "-CleanStripped.fasta")

        #remove duplicate sequences from the fasta file and store for later
        print "Remove duplicates"
        secs = time()
        kept, headers = remove_duplicates(kept)
        #parse out only sequences with enough duplicates
        if args.d > 1:
            newkept = [seq for seq in kept if int(seq[0].split("_")[1]) >= args.d]
            kept = newkept
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
        print str(len(kept)) + " sequences left. " + str((time() - secs)/60) + " minutes\n"
        log.close()
