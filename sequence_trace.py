from cogent.parse.fasta import MinimalFastaParser
from cogent.align.algorithm import nw_align
import argparse
#requires http://code.google.com/p/pylevenshtein/
from Levenshtein import distance
from multiprocessing import Pool, Manager


def find_ends(seq):
    '''Finds where the gaps end in the passed sequence and returns positions'''
    left = -1
    right = -1
    for pos in range(len(seq)):
        if seq[pos] != "-":
            left = pos
            break
    for pos in reversed(range(len(seq))):
        if seq[pos] != "-":
            right = pos
            break
    return left, right


def compare_seqs(knownseq, testseq, header):
    #keep if exact match of substrings
    if knownseq in testseq or testseq in knownseq:
        return ">%s\n%s" % (header, testseq)
    #START: keep if one difference between each sequence
    alnknown, alnseq = nw_align(knownseq, testseq)
    #find where to start comparison by trimming gaps at ends
    knownstart, knownend = find_ends(alnknown)
    seqstart, seqend = find_ends(alnseq)
    start = min(seqstart, knownstart)
    end = max(seqend, knownend)
    if start < 0 or end < 0:
        raise ValueError("Start and end must be greater than 0!")
    #if only one difference, print out sequence
    if distance(knownseq[start:end], testseq[start:end]) <= 1:
        return ">%s\n%s" % (header, testseq)
    return ""


def multi_wrapper(knownseq, testseq, header, lock, fout):
    try:
        comp = compare_seqs(knownseq, testseq, header)
        if comp != "":
            lock.acquire()
            fout.write(comp + "\n")
            lock.release()
    finally:
        lock.release()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Finds instances of ")
    parser.add_argument('-s', required=True, help="Sequence being serched for \
        over all rounds.")
    parser.add_argument('-f', required=True, help="Base folder holding FASTA \
        files for all unique sequences in rounds of selection ")
    parser.add_argument('-r', required=True, type=int, help="Round of selection \
        the sequence comes from")
    parser.add_argument('-c', default=1, type=int, help="Number of CPUs \
        (Default 1)")

    args = parser.parse_args()
    basefolder = args.f
    if basefolder[-1] != "/":
        basefolder += "/"
    knownseq = args.s
    rnd = args.r
    if args.c < 1:
        raise ValueError("CPUs must be greater than 0!")

    #loop over all the rounds to find sequence matches
    fout = open("matches.fasta", 'w')
    fout.write(">%s\n%s\n" % ("SeedSequence", knownseq))
    for currrnd in range(1, rnd+1):
        rndname = "R" + str(currrnd)
        print rndname
        fout.write(">%s\n%s\n" % (rndname, rndname))
        #multiprocess each round
        manager = Manager()
        hold = manager.dict()
        pool = Pool(processes=args.c)
        lock = manager.Lock()
        #run the pool over all sequences for comparison
        rndfile = open(basefolder + rndname + "/" + rndname + "-Unique.fasta")
        for header, seq in MinimalFastaParser(rndfile):
            #add another for loop here to go through multiple knownseqs
            pool.apply_async(func=multi_wrapper, args=(knownseq, seq, header, lock, fout))
        rndfile.close()
        pool.close()
        pool.join()
    fout.close()
