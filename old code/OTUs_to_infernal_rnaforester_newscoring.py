#dupestest.py /path/to/OTU.fasta /path/to/folder/out CPUs
from sys import exit
from sys import stdout
from os.path import exists
from os import mkdir, remove
from cogent import LoadSeqs, RNA
from cogent.app.locarna import create_locarnap_alignment
from cogent.app.infernal_v11 import cmsearch_from_file, cmbuild_from_alignment
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from remove_duplicates import remove_duplicates
from subprocess import Popen, PIPE
from time import time, sleep
from multiprocessing import Pool, Lock, Manager
import argparse
from math import floor

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
            #print otu + ": " + struct
            #stdout.flush()
            lock.release()
    finally:
        lock.release()


def fold_groups(otus, clusters, struct, hold):
    '''Function for multithreading.
    Recomputes structure for a group'''
    #no need to refold if only one cluster in group
    if len(clusters) == 1:
        hold[struct] = structgroups[struct]
    else:
        #compute new structures for each group
        aln, currstruct = run_locarnap_groups(otus, clusters, 25, cpus=2)
        aln = 0
        if currstruct == "":
            currstruct = struct
        hold[currstruct] = structgroups[struct]
        count += 1


def run_locarnap(seqsin, numkept, cpus=1,foldless=False):
    '''Runs locarna-p on a set of sequences in MinimalFastaParser format
    [(header, seq), (header, seq)] and retgurns alignemtn and structure'''
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


def run_rnaforester(struct1, struct2):
    '''returns global optimal score of two structures'''
    p = Popen(["rnaforester", "--score"], stdin=PIPE, stdout=PIPE)
    p.stdin.write(''.join([struct1, "\n", struct2, "\n&"]))
    return int(p.communicate()[0].split("\n")[-2])

def make_r2r(insto, outfolder, group):
    '''generates R2R secondary structure pdf with dfault colorings'''
    p = Popen(["r2r", "--GSC-weighted-consensus", insto, outfolder + "/" + group + ".sto", \
        "3", "0.97", "0.9", "0.75", "4", "0.97", "0.9", "0.75", "0.5", "0.1"])
    sleep(1)
    p = Popen(["r2r", outfolder + "/" + group + ".sto", outfolder + "/" + group + ".pdf"], stdout = PIPE)


def cluster_rnaforester(structgroups, foresterscore):
    hold = structgroups.keys() #now hold contains list of structures to compare
    #list of positions of maximal match to the same structure in hold position
    maxscores = []
    #populate the score list
    for row, currstruct in enumerate(hold): #for each structure
        currmax = -1
        bestpos = -1
        for col, teststruct in enumerate(hold): #compare it to all other structures
            score = run_rnaforester(currstruct, teststruct)
            #keep the highest score as the best match
            if row == col:
                continue
            if score >= foresterscore and score > currmax:
                currmax = score
                bestpos = col
        if currmax == -1:
            print "APPEND ROW: " + str(row)
            maxscores.append(row)
        else:
            print "APPEND BEST: " + str(bestpos)
            maxscores.append(bestpos)
    #create new dictionary of structgroups
    newstructgroups = {}
    #for each structure, append it's max score match to the group
    for row, col in enumerate(maxscores):
        if not hold[col] in newstructgroups:
            newstructgroups[hold[col]] = structgroups[hold[row]] 
        else:
            newstructgroups[hold[col]].extend(structgroups[hold[row]])
    print newstructgroups
    return newstructgroups


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs structural clustering \
    and infernal over all rounds onf a SELEX selection")
    parser.add_argument('-i', required=True, help="TXT file of clusters and the\
    paths to their FASTA files")
    parser.add_argument('-f', required=True, help="Base folder holding FASTA files for all \
    unique sequences in the rounds of selection (for infernal)")
    parser.add_argument('-o', required=True, help="Base folder to output all data to")
    parser.add_argument('-r', required=True, type=int, help="Current round of selection\
    clusters come from")
    parser.add_argument('--rsc', type=int, default=-10000, help="Score cutoff for rnaforester \
    (Default 2 * max sequence length)")
    parser.add_argument('--isc', type=float, default=0.0, help="Score cutoff for infernal.\
    (Default 0.0)")
    parser.add_argument('-c', type=int, default=1, help="Number of CPUs to use \
    (Default 1)")
    parser.add_argument('-inf', action="store_true", help="Skip straight to infernal (Default False)")

    args = parser.parse_args()
    if args.r < 1:
        print "ERROR: round must be at least 1!"
        exit(1)
    if args.c < 1:
        print "ERROR: CPU count must be at least 1!"
        exit(1)
    if args.isc < 0:
        print "ERROR: Infernal score cutoff must be greater than 0!"
        exit(1)
    compscore = False
    if args.rsc == -10000:
        compscore = True

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
    count = 1
    secs = time()

    #print otu + " >> START"

    #make a pool of workers, one for each cpu available
    manager = Manager()
    pool = Pool(processes=args.c)
    lock = manager.Lock()
    #run the pool over all clusters to get file of structures
    for otu in otus:
        pool.apply_async(func=fold_clusters, args=(lock, known, otu, otufolder))
    pool.close()
    pool.join()

    structgroups = {}
    #read in all structures and populate structgroups dict
    cfo = open(otufolder + "cluster_structs.fasta", 'rU')
    for cluster, struct in MinimalFastaParser(cfo):
        if struct in structgroups:
            structgroups[struct].append(cluster)
        else:
            structgroups[struct] = [cluster]
        if compscore and len(struct) > foresterscore:
            foresterscore = len(struct)
        count += 1
    cfo.close()

    print "Runtime: " + str((time() - secs)/60) + " min"
    print "==Grouping clusters by secondary structure=="
    #sequences must match 1.5 * longest seq length
    foresterscore = 1.5 * foresterscore
    #GROUP THE SECONDARY STRUCTURES BY RNAFORESTER
    #logic: first item in list is always unique, so use as base
    #for comparison. That way when list empty, all are grouped
    secs = time()
    skipiter = False
    #if groups already exist, read them in
    if exists(otufolder + "groups.txt"):
        gin = open(otufolder + "groups.txt", 'rU')
        #get score from file
        foresterscore = float(gin.readline().strip())
        for line in gin:
            if line.strip() == "FINAL": 
                skipiter = True
            else:
                lineinfo = line.split(':')
                structgroups[lineinfo[0]] = lineinfo[1].strip().split()
        gin.close()
        print "RNAforester score threshold: " + str(foresterscore)
    else:
        print "RNAforester score threshold: " + str(foresterscore)
        structgroups = cluster_rnaforester(structgroups, foresterscore)
        #write out file with groups
        gout = open(otufolder + "groups.txt", 'w')
        gout.write(str(foresterscore) + "\n")
        for group in structgroups:
            gout.write(group + ":")
            for g in structgroups[group]:
                gout.write(g+" ")
            gout.write("\n")
        gout.close()
        print "Initial grouping (" + str((time() - secs)/60)+ "m)"
    #Now need to iteratively refine the groups down
    #check to make sure we need to first
    if not args.inf and not skipiter:
        startcount = 1
        endcount = 0
        iteration = 1
        secs = time()

        while startcount != endcount: #keep refining while we are still grouping structs
            startcount = len(structgroups)
            print "iteration " + str(iteration) + ": " + str(startcount) + " groups"

            #Refold all the groups to get new consensus secondary structure
            #make a pool of workers, one for each cpu available
            manager = Manager()
            hold = manager.dict()
            #make sure to check if we have only 1 cpu
            procs = int(floor(args.c/2))
            if procs == 0:
                procs = 1
            pool = Pool(processes=procs)
            #run the pool over all groups to get structures
            for struct in structgroups:
                pool.apply_async(func=fold_groups, args=(otus, structgroups[struct], struct, hold))
            pool.close()
            pool.join()
            
            #need to turn manager dict to regular dict
            structgroups = {}
            for key in hold:
                structgroups[key] = hold[key]

            structgroups = cluster_rnaforester(structgroups, foresterscore)

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

        print str(len(structgroups)) + " final groups (" + str((time() - secs)/60)+ "m)"
        gout = open(otufolder + "groups.txt", 'w')
        gout.write("FINAL\n")
        gout.close()
        #end while
    else:
        print str(len(structgroups)) + " final groups"
    print "==Creating CM and running Infernal over all rounds=="
    print "Infernal score cutoff: " + str(infernalscore)
    currgroup = 1
    for currstruct in structgroups:
        secs = time()
        #run locana-p on the superclusters to get the alignment and consensus structure
        seqs = []
        out = "group " + str(currgroup)+ ": "
        for cluster in structgroups[currstruct]:
            out += cluster+ " "
            for header,seq in MinimalFastaParser(open(otus[cluster], 'rU')):
                seqs.append((header, seq))
        print out
        print str(len(seqs)) + " sequences"
        #make sure group has enough sequences before continuing
        #run locarna-p on the at most 50 most abundant sequences in the group
        aln, struct = run_locarnap(seqs,50, cpus=args.c, foldless=True)
        if(aln.getNumSeqs() < 50):
            print str(aln.getNumSeqs()) + " unique sequences"
        else:
            s, h = remove_duplicates(seqs)
            print str(len(s)) + " unique sequences"
            s = 0
            h = 0
        print "Structure: " + struct

        #print out alignment and structure in fasta and stockholm formats
        #create output folder for group
        currotufolder = otufolder + "group_" + str(currgroup)
        if not exists(currotufolder):
            mkdir(currotufolder)
        #print that shit
        alnout = open(currotufolder + "/locarnap-aln.fasta", 'w')
        alnout.write(aln.toFasta() + "\n>SS_struct\n" + struct + "\n")
        alnout.close()
        alnout = open(currotufolder + "/locarnap-aln.sto", 'w')
        struct_dict = {'SS_cons':struct}
        alnout.write(stockholm_from_alignment(aln,GC_annotation=struct_dict))
        alnout.close()
        #make R2R secondary structure for alignment
        make_r2r(currotufolder + "/locarnap-aln.sto", currotufolder, "group_" + str(currgroup))
        #only run infernal if there were more than 100 total sequences in group
        #essentially need at least 0.1% of pool in cluster. save time this way
        if len(seqs) > 99:
            #create the cm file. Could call cmsearch_from_alignment but dont want to build
            #cm file multiple times since is time consuming and processor intensive
            cmfile = cmbuild_from_alignment(aln, struct,calibrate=True)
            for i in range(7, 0, -1):
                #run cmsearch over every round of SELEX
                #Only search unique sequences to save time
                seqs = 0
                #check if previous run has removed some sequences, load correct file
                if exists(args.f + "R" + str(i) + "/R" + str(i) + "-Unique-Remaining.fasta"):
                    print "Previous round run detected, runnning over remaining seqs"
                    seqs = LoadSeqs(args.f + "R" + str(i) + "/R" + str(i) + "-Unique-Remaining.fasta", moltype=RNA, aligned=False)
                else:
                    seqs = LoadSeqs(args.f + "R" + str(i) + "/R" + str(i) + "-Unique.fasta", moltype=RNA, aligned=False)
                result = cmsearch_from_file(cmfile.name, seqs, RNA, cutoff=infernalscore, params={'--toponly' : True, '--cpu' : args.c})
                print "R" + str(i) + ": " + str(len(result)) + " hits"
                fout = open(currotufolder + "/R" + str(i) + "hits.txt", 'w')
                fout.write("header,bitscore,e-value\n")
                for hit in result:
                    fout.write(hit[0] + "," + str(hit[14]) + "," + str(hit[15]) + "\n")
                fout.close()
            #clean up by removing cm file
            remove(cmfile.name)
                #remove found sequences from the round files
            print "Runtime: " + str(time() - secs)+ "m"
            currgroup += 1
        #skip group if less than 100 sequences
        else:
            print "Group has less than 100 sequences, skipping infernal run"
            currgroup += 1
     
