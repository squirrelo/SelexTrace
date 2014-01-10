from os.path import exists, splitext
from os import mkdir, walk
from sys import exit
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from time import time
import argparse
from selextrace.stutils import write_fasta_list, strip_primer,\
    remove_duplicates

def rem_N_short(seqs, minlen=1):
    '''Takes in a [(header, seq)] formatted list of sequences and returns list 
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
    parser.add_argument('-ep', default="", help="3' primer \
    sequence to strip")
    #parser.add_argument('-sp', default="", help="5' primer \
    #sequence to strip")
    parser.add_argument('-o', default = "", help="Output folder (default same as input)")
    parser.add_argument('-l', default = 1, type=int, help="minimum length of \
    sequence to keep (default 1)")
    parser.add_argument('-d', default = 1, type=int, help="minimum number of duplicates\
    needed to keep (default 1)")
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
    print "Folder in:", folderin
    print "Output Folder:", folderout
    print "3' primer:", args.ep
    #print "5' primer:", args.sp
    print "Min length:", args.l
    print "Min duplicates:", args.d
    print "====================\n"
    for filein in walk(folderin).next()[2]:
        fext = splitext(filein)[1]
        #skip if not fastq or fasta file
        if fext != ".fastq" and fext != ".fasta" and fext != ".fas" and fext != ".fna":
            continue

        basename = splitext(filein)[0]

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

        print "==Cleaning input sequences=="

        log = open(currfolder + "-cleanup.log", 'w')
        log.write(''.join(["====================\nFile in: ", folderin, filein,
            "\nOutput Folder: ", currfolder, "\n3' primer: ", args.ep, 
            "\nMin length: " + str(args.l), "\nMin duplicates: ", str(args.d), 
            "\n====================\n"]))
        #remove all underscores from headers during load for compatability reasons
        seqs = []
        seqsin = open(folderin + filein, 'rU')
        for header, seq in MinimalFastaParser(seqsin):
            seqs.append((header, seq))
        print len(seqs), "starting sequences"
        log.write(len(seqs), "starting sequences\n")
        #strip primers from sequences, print out not stripped to be safe
        #allowing up to 2 mismatches in the primer
        if args.ep != '': #or args.sp != '':
            print "Primer stripping"
            secs = time()
            kept, rem = strip_primer(seqs, args.ep, maxmismatch=2, 
                keep_primer=True)
            del seqs
            log.write("Primer stripping\n" + str(len(kept)) + " sequences left\n")
            print str(len(kept)) + " sequences left, " + \
            str(len(rem)) + " sequences removed. " + str((time() - secs)/60) + " minutes\n"
            write_fasta_list(kept, currfolder + "-Stripped.fasta")
            write_fasta_list(rem, currfolder + "-NotStripped.fasta")
            del rem
        else:
            kept = seqs
            del seqs
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
        kept = remove_duplicates(kept)
        #parse out only sequences with enough duplicates
        if args.d > 1:
            stop = 0
            #already sorted most->least, so first hit with more means we have
            #found the cut point
            while kept[stop][1] >= args.d:
                stop += 1
            kept = kept[:stop]

        fout = open(currfolder + "-Unique.fasta", 'w')
        for num, seqinfo in enumerate(kept):
            seq = seqinfo[0]
            count = str(seqinfo[1])
            fout.write('>seq%s count:%s\n%s\n' % (str(num), count, seq))
        fout.close()
        log.write("Remove duplicates\n" + str(len(kept)) + " sequences left")
        print str(len(kept)) + " sequences left. " + str((time() - secs)/60) + " minutes\n"
        log.close()
