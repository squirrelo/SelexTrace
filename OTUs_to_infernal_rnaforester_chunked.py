#dupestest.py /path/to/OTU.fasta /path/to/folder/out CPUs
import sys
from os.path import exists
from os import mkdir, remove, walk
from cogent import LoadSeqs, RNA
from cogent.app.locarna import create_locarnap_alignment
from cogent.app.infernal_v11 import cmsearch_from_file, cmbuild_from_file, calibrate_file
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from remove_duplicates import remove_duplicates
from subprocess import Popen, PIPE
from time import time, sleep
from datetime import datetime
from multiprocessing import Pool, Lock, Manager
import argparse
from math import floor, ceil
from random import shuffle


def count_lines(f):
    '''returns number of lines in a file'''
    return sum(1 for line in open(f, 'U'))


def fold_clusters(lock, known, otu, otufolder):
    '''Function for multithreading.
    Computes structure for a cluster and writes it to file'''
    #assuming that the fasta has 10 or more sequences in it. Safe assumption
    #if this is a significant cluster
    #only using 10 because this is initial structure calc so needs to be fast
    #and this gives a very good approximation for the initial rnaforester grouping
    try:
        struct = ""
        if otu in known:
            return 0
        else:
            filein = open(otus[otu], 'rU')
            seqs = [(header, seq) for header,seq in MinimalFastaParser(filein)]
            filein.close()
            aln, struct = run_locarnap(seqs, 10, foldless=True)
            #write structure out to file
            lock.acquire()
            cfo = open(otufolder + "cluster_structs.fasta", 'a')
            cfo.write(">" + otu + "\n" + struct+ "\n")
            cfo.close()
            lock.release()
    finally:
        lock.release()


def fold_groups(otus, clusters, struct, hold):
    '''Function for multithreading.
    Recomputes structure for a group'''
    #compute new structures for each group
    aln, currstruct = run_locarnap_groups(otus, clusters, 25, cpus=2, foldless=True)
    aln = 0
    if currstruct == "":
        currstruct = struct
    hold[currstruct] = structgroups[struct]
    count += 1


def run_locarnap(seqsin, numkept, cpus=1,foldless=False):
    '''Runs locarna-p on a set of sequences in format
    [(header, seq), (header, seq)] and returns alignment and structure'''
    seqs, headers = remove_duplicates(seqsin)
    #blank headers to save memory
    headers = 0
    #make sure group has enough sequences before continuing
    if len(seqs) < numkept and not foldless:
        return "", ""
    #headers come out in format Header_# so split to get # and sort by abundance
    seqs.sort(reverse=True, key=lambda count:int(count[0].split('_')[1]))
    #cut to numkept most abundant sequences
    if len(seqs) > numkept:
        seqs = seqs[:numkept]
    return create_locarnap_alignment(seqs, RNA, struct=True, params={'--cpus':cpus})


def run_locarnap_groups(otulist, clusters, numkept, cpus=1, foldless=False):
    '''runs locarna-p over groups of clusters to get alignment and 
    consensus secondary structure'''
    seqs = []
    for cluster in clusters:
        filein = open(otulist[cluster], 'rU')
        for header,seq in MinimalFastaParser(filein):
            seqs.append((header, seq))
        filein.close()

    return run_locarnap(seqs, 30, cpus, foldless)

def run_locarnap_for_infernal(currgroup, clusters, otus, basefolder):
    '''Function for multithreading
    creates the final locarna-p alignment and writes to files, then r2r struct'''
    #run locana-p on the superclusters to get the alignment and consensus structure
    #skip if already run and program just crashsed or whatever
    currotufolder = basefolder + "group_" + str(currgroup)
    if exists(currotufolder):
        return ""
    seqs = []
    out = "group " + str(currgroup)+ ": "
    for cluster in clusters:
        out += cluster+ " "
        for header,seq in MinimalFastaParser(open(otus[cluster], 'rU')):
            seqs.append((header, seq))
    out += "\n" + str(len(seqs)) + " sequences\n"
    #make sure group has enough sequences before continuing
    #run locarna-p on the at most 50 most abundant sequences in the group
    aln, struct = run_locarnap(seqs, 50, cpus=2, foldless=True)
    if(aln.getNumSeqs() < 50):
        out += str(aln.getNumSeqs()) + " unique sequences\n"
    else:
        s, h = remove_duplicates(seqs)
        out += str(len(s)) + " unique sequences\n"
        s = 0
        h = 0
    out += "Structure: " + struct + "\n"

    #write out alignment and structure in fasta and stockholm formats
    #create output folder for group
    mkdir(currotufolder)
    #write that shit
    logout = open(currotufolder + "/log.txt", 'w')
    logout.write(out)
    logout.close()
    alnout = open(currotufolder + "/locarnap-aln.fasta", 'w')
    alnout.write(aln.toFasta() + "\n>SS_struct\n" + struct + "\n")
    alnout.close()
    alnout = open(currotufolder + "/locarnap-aln.sto", 'w')
    struct_dict = {'SS_cons':struct}
    alnout.write(stockholm_from_alignment(aln,GC_annotation=struct_dict))
    alnout.close()
    #make R2R secondary structure for alignment
    make_r2r(currotufolder + "/locarnap-aln.sto", currotufolder, "group_" + str(currgroup))


def score_rnaforester(struct1, struct2):
    '''returns relative score of two structures'''
    p = Popen(["rnaforester", "--score", "-r"], stdin=PIPE, stdout=PIPE)
    p.stdin.write(''.join([struct1, "\n", struct2, "\n&"]))
    return float(p.communicate()[0].split("\n")[-2])


def group_by_forester(structgroups, foresterscore):
    topop = set([])
    for currstruct in structgroups: #for each structure
        if currstruct in topop: #skip if we already grouped this one
            continue
        for teststruct in structgroups: #else compare it to all other structures
            #only compare if not equal or not already grouped
            if teststruct != currstruct and teststruct not in topop:
                score = score_rnaforester(currstruct, teststruct)
                if score > foresterscore:
                    structgroups[currstruct].extend(structgroups[teststruct])
                    structgroups[teststruct] = 0
                    topop.add(teststruct)
    for key in topop: #pop all simmilar structures found
        structgroups.pop(key)
    return structgroups


def make_r2r(insto, outfolder, group):
    '''generates R2R secondary structure pdf with default colorings'''
    p = Popen(["r2r", "--GSC-weighted-consensus", insto, outfolder + "/" + group + ".sto", \
        "3", "0.97", "0.9", "0.75", "4", "0.97", "0.9", "0.75", "0.5", "0.1"])
    sleep(1)
    p = Popen(["r2r", outfolder + "/" + group + ".sto", outfolder + "/" + group + ".pdf"], stdout = PIPE)


def run_infernal(lock, cmfile, rnd, seqs, outfolder, cpus=1, score=0.0):
    try:
        result = cmsearch_from_file(cmfile, seqs, RNA, cutoff=score, params={'-g' : True, '--notrunc' : True, '--nohmm': True, '--toponly' : True, '--cpu' : cpus})
        fout = ""
        lock.acquire()
        if exists(outfolder + "/R" + str(rnd) + "hits.txt"):
            fout = open(outfolder + "/R" + str(rnd) + "hits.txt", 'a')
        else:
            fout = open(outfolder + "/R" + str(rnd) + "hits.txt", 'w')
        fout.write(str(len(result)) + " hits\nheader,bitscore,e-value\n")
        for hit in result:
            fout.write(hit[0] + "," + str(hit[14]) + "," + str(hit[15]) + "\n")
        fout.close()
        lock.release()
    except Exception, e:
        lock.acquire()
        print "ERROR:"
        print sys.exc_info()[0]
        sys.stdout.flush()
        lock.release()
    finally:
        lock.release()

def split_fasta(fasta, numsplit, mtype):
    '''Splits a fasta file into numsplit chunks and returns a list of SequenceCollection
    objects for each chunk'''
    count = 0
    totalsplit = 1
    seqdict = {}
    chunks = []
    fastafile = open(fasta, 'U')
    numseqs = count_lines(fasta) / 2
    splitat =  numseqs / numsplit
    fastafile.close()
    fastafile = open(fasta, 'U')
    for header, seq in MinimalFastaParser(fastafile):
        seqdict[header] = seq
        count += 1
        #need totalsplit != numsplit in case there are extra sequences 
        if count == splitat and totalsplit != numsplit:
            count = 0
            chunks.append(LoadSeqs(data=seqdict, moltype=mtype, aligned=False))
            seqdict = {}
            totalsplit += 1
    fastafile.close()
    chunks.append(LoadSeqs(data=seqdict, moltype=mtype, aligned=False))
    return chunks



if __name__ == "__main__":
    starttime = time()
    parser = argparse.ArgumentParser(description="Runs structural clustering \
    and infernal over all rounds of a SELEX selection")
    parser.add_argument('-i', required=True, help="TXT file of clusters and the\
    paths to their FASTA files")
    parser.add_argument('-f', required=True, help="Base folder holding FASTA files for all \
    unique sequences in the rounds of selection (for infernal)")
    parser.add_argument('-o', required=True, help="Base folder to output all data to")
    parser.add_argument('-r', required=True, type=int, help="Current round of selection\
    clusters come from")
    parser.add_argument('--rsc', type=float, default=0.5, help="Relative score cutoff for rnaforester \
    (Default 0.5)")
    parser.add_argument('--isc', type=float, default=0.0, help="Score cutoff for infernal.\
    (Default 0.0)")
    parser.add_argument('-c', type=int, default=1, help="Number of CPUs to use \
    (Default 1)")
    parser.add_argument('-inf', action="store_true", help="Skip straight to infernal (Default False)")

    args = parser.parse_args()
    if args.r < 1:
        print "ERROR: round must be at least 1!"
        sys.exit(1)
    if args.c < 1:
        print "ERROR: CPU count must be at least 1!"
        sys.exit(1)
    if args.isc < 0:
        print "ERROR: Infernal score cutoff must be greater than 0!"
        sys.exit(1)
    if args.rsc < 0 or args.rsc > 1:
        print "ERROR: RNAforester score cutoff must be between 0 and 1!"
        sys.exit(1)

    foresterscore = args.rsc
    infernalscore = args.isc
    otufolder = args.o

    if args.f[:-1] != "/":
        args.f += "/"

    if otufolder[:-1] != "/":
        otufolder += "/"
    if not exists(otufolder):
        mkdir(otufolder)

    #Get list of clusters from file, populate dictionary
    #path to cluster's fasta keyed to cluster name
    print "Program started ", datetime.now()
    otus = {}
    fn = open(args.i, 'rU')
    for line in fn:
        lineinfo = line.split()
        if not exists(lineinfo[1].strip()):
            print lineinfo[0] + ": OTU FASTA file missing!"
            exit(1)
        otus[lineinfo[0]] = lineinfo[1].strip()
    fn.close()
    known = {}
    if exists(otufolder + "cluster_structs.fasta"):
            #read in known structures
            csf = open(otufolder + "cluster_structs.fasta")
            for header, struct in MinimalFastaParser(csf):
                known[header] = struct
            csf.close()
    else:
        #create file to write to if not already there
        cfo = open(otufolder + "cluster_structs.fasta", 'w')
        cfo.close()
    print "==Clustering sequences by structure=="
    print "Running locarna-p over "+ str(len(otus)) +" clusters"
    secs = time()

    #make a pool of workers, one for each cpu available
    manager = Manager()
    pool = Pool(processes=args.c)
    lock = manager.Lock()
    #run the pool over all clusters to get file of structures
    for otu in otus:
        pool.apply_async(func=fold_clusters, args=(lock, known, otu, otufolder))
    pool.close()
    pool.join()
    #read in all structures now that they are folded
    structgroups = {}
    cfo = open(otufolder + "cluster_structs.fasta", 'rU')
    for cluster, struct in MinimalFastaParser(cfo):
        if struct in structgroups:
            structgroups[struct].append(cluster)
        else:
            structgroups[struct] = [cluster]
    cfo.close()

    print "Runtime: " + str((time() - secs)/60) + " min"
    print "==Grouping clusters by secondary structure=="
    #GROUP THE SECONDARY STRUCTURES BY RNAFORESTER
    #logic: first item in list is always unique, so use as base
    #for comparison. That way when list empty, all are grouped
    secs = time()
    skipiter = False
    #if groups already exist, read them in
    if exists(otufolder + "groups.txt"):
        structgroups = {}
        gin = open(otufolder + "groups.txt", 'rU')
        foresterscore = float(gin.readline().strip())
        for line in gin:
            if line.strip() == "FINAL": 
                skipiter = True
            else:
                lineinfo = line.split(':')
                structgroups[lineinfo[0]] = lineinfo[1].strip().split()
        gin.close()
    print "RNAforester score threshold: " + str(foresterscore)
    #Now need to iteratively refine the groups down
    #check to make sure we need to first
    if not args.inf and not skipiter:
        startcount = 1
        endcount = 0
        iteration = 0
        secs = time()
        print "iteration " + str(iteration) + ": " + str(len(structgroups)) + " initial groups"
        #initial clustering by structures generated in first folding
        structgroups = group_by_forester(structgroups, foresterscore)
        iteration += 1

        while startcount != endcount: #keep refining while we are still grouping structs
            startcount = len(structgroups)
            print "iteration " + str(iteration) + ": " + str(len(structgroups)) + " initial groups"
            #Refold all the groups to get new consensus secondary structure
            #make a pool of workers, one for each cpu available
            manager = Manager()
            hold = manager.dict()
            #make sure to check if we have only 1 cpu
            procs = int(floor(args.c/2))
            if procs == 0:
                procs = 1
            pool = Pool(processes=procs)
            #run the pool over all groups to get structures for new groups
            for struct in structgroups:
                pool.apply_async(func=fold_groups, args=(otus, structgroups[struct], struct, hold))
            pool.close()
            pool.join()
            
            #need to turn manager dict to regular dict
            structgroups = {}
            for key in hold.keys():
                structgroups[key] = hold[key]

            structgroups = group_by_forester(structgroups, foresterscore)
            endcount = len(structgroups)
            #write out file with groups
            gout = open(otufolder + "groups.txt", 'w')
            gout.write(str(foresterscore) + "\n")
            for group in structgroups:
                gout.write(group + ":")
                for g in structgroups[group]:
                    gout.write(g+" ")
                gout.write("\n")
            gout.close()
            iteration += 1
        #end while

        print str(len(structgroups)) + " final groups (" + str((time() - secs)/60)+ "m)"
        gout = open(otufolder + "groups.txt", 'a')
        gout.write("FINAL\n")
        gout.close()
    else:
        print str(len(structgroups)) + " final groups"

    print "==Creating CM and r2r structures=="
    secs = time()
    pool = Pool(processes=int(ceil(args.c/2)))
    #run the pool over all clusters to get file of structures
    for group,currstruct in enumerate(structgroups):
        pool.apply_async(func=run_locarnap_for_infernal, args=(group+1, structgroups[currstruct], otus, otufolder))
    pool.close()
    pool.join()
    print "Runtime: " + str((time() - secs)/60)+ "m"
    print "==Running Infernal for all groups=="
    print "Infernal score cutoff: " + str(infernalscore)
    #create teh csv file for holding all teh hit counts
    ihits = open(otufolder + "infernalhits.csv", 'w')
    ihits.write(",")
    for i in range(1, args.r+1):
        ihits.write("round" + str(i) + ",")
    ihits.write("\n")

    for group in walk(otufolder).next()[1]:
        secs = time()
        print "GROUP: " + group
        logfile = open(otufolder +  group + "/log.txt")
        log = logfile.readlines()
        logfile.close()
        print ''.join(log)
        seqs = int(log[1].split()[0])
        skip = False
        if exists(otufolder +  group + "/R1hits.txt"):
            skip = True
        #only run infernal if there were more than 100 total sequences in group
        if seqs > 99 and not skip:
            currotufolder = otufolder + group
            #create the cm file and calibrate it 
            if not exists(otufolder +  group + "/infernal.cm"):
                cmfile = open(otufolder +  group + "/infernal.cm", 'w')
                cmfile.write(cmbuild_from_file(otufolder +  group + "/locarnap-aln.sto"))
                cmfile.close()
                cmfile = otufolder +  group + "/infernal.cm"
                calibrate_file(cmfile)

            cmfile = otufolder +  group + "/infernal.cm"
            #Run all rounds of selection through infernal at once
            #make a pool of workers, one for each cpu available
            manager = Manager()
            lock = manager.Lock()
            #calculate maximal amount of CPU power and theads we can use
            procs = int(floor(args.c/args.r))
            if procs == 0:
                procs = 1
            poolsize = args.r
            if args.c < args.r:
                poolsize = args.c
            pool = Pool(processes=poolsize)
            extracpus = args.c - (procs*poolsize)
            for rnd in range(args.r, 0, -1):
                #check if previous run has removed some sequences, load correct file
                if exists(args.f + "R" + str(rnd) + "/R" + str(rnd) + "-Unique-Remaining.fasta"):
                    seqs = args.f + "R" + str(rnd) + "/R" + str(rnd) + "-Unique-Remaining.fasta"
                else:
                    seqs = args.f + "R" + str(rnd) + "/R" + str(rnd) + "-Unique.fasta"
                #chunk the sequences so infernal will run faster
                #if there is extra cpu power, apply it!
                if extracpus > 0:
                    chunks = split_fasta(seqs, procs+1, RNA)
                    extracpus -= 1
                else:
                    chunks = split_fasta(seqs, procs, RNA)
                for chunk in chunks:
                #run cmsearch over every round of SELEX
                    pool.apply_async(func=run_infernal, args=(lock, cmfile, rnd, chunk, currotufolder, 1, infernalscore))
            pool.close()
            pool.join()
            #count the number of hits per round and write it out to files
            roundhits = []
            for i in range(1, args.r+1):
                roundhits.append(count_lines(open(currotufolder + "/R" + str(i) + "hits.txt", 'U'))-1)
            logfile = open(otufolder +  group + "/log.txt", 'a')
            for rnd, hits in enumerate(roundhits):
                print "Round" + str(rnd+1) + ": " +  str(hits) + " hits"
                logfile.write("Round" + str(rnd+1) + ": " +  str(hits) + " hits")
            hitscsv = open(otufolder + "infernalhits.csv", 'a')
            hitscsv.write(group + ",")
            for hits in roundhits:
                hitscsv.write(str(hits) + ",")
            hitscsv.write("\n")
            print "Runtime: " + str((time() - secs)/60)+ "m"
        elif skip:
            print "Group already run"
            logfile = open(otufolder +  group + "/log.txt")
            roundhits = logfile.readlines()[4:]
            logfile.close()
            roundhits.sort()
            hitscsv = open(otufolder + "infernalhits.csv", 'a')
            hitscsv.write(group + ",")
            for r in roundhits:
                hitscsv.write(r.split()[2] + ",")
            hitscsv.write("\n")
        #skip group if less than 100 sequences
        else:
            print "Group has less than 100 sequences, skipping infernal run"
print "Program ended ", datetime.now(), "   Runtime: " + str((time() - starttime)/3600)+ "h" 
     
