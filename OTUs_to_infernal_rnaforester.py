#dupestest.py /path/to/OTU.fasta /path/to/folder/out CPUs
import sys
from os.path import exists
from os import mkdir, walk
from clean_seqs import write_fasta_list
from cogent import LoadSeqs, RNA
from cogent.app.locarna import create_locarnap_alignment
from cogent.app.infernal_v11 import cmsearch_from_file, cmbuild_from_file, calibrate_file
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from remove_duplicates import remove_duplicates
from subprocess import Popen, PIPE
from time import time
from datetime import datetime
from multiprocessing import Pool, Manager
import argparse
from math import floor, ceil


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
            seqs = [(header, seq) for header, seq in MinimalFastaParser(filein)]
            filein.close()
            aln, struct = run_locarnap(seqs, 10, foldless=True)
            #write structure out to file
            lock.acquire()
            cfo = open(otufolder + "cluster_structs.fasta", 'a')
            cfo.write(">" + otu + "\n" + struct + "\n")
            cfo.close()
            lock.release()
    finally:
        lock.release()


def fold_groups(otus, clusters, struct, hold):
    '''Function for multithreading.
    Recomputes structure for a group'''
    #compute new structures for each group
    aln, currstruct = run_locarnap_groups(otus, clusters, 25, foldless=True)
    if currstruct == "":
        currstruct = struct
    hold[currstruct] = structgroups[struct]


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
    seqs.sort(reverse=True, key=lambda count: int(count[0].split('_')[1]))
    #cut to numkept most abundant sequences
    if len(seqs) > numkept:
        seqs = seqs[:numkept]
    aln, struct = create_locarnap_alignment(seqs, RNA, struct=True, params={'--cpus': cpus})
    struct = struct.replace('-', ".")
    return aln, struct


def run_locarnap_groups(otulist, clusters, numkept, cpus=1, foldless=False):
    '''runs locarna-p over groups of clusters to get alignment and
     consensus secondary structure'''
    seqs = []
    for cluster in clusters:
        filein = open(otulist[cluster], 'rU')
        for header, seq in MinimalFastaParser(filein):
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
    out = "group " + str(currgroup) + ": "
    for cluster in clusters:
        out += cluster + " "
        for header, seq in MinimalFastaParser(open(otus[cluster], 'rU')):
            seqs.append((header.split()[0], seq))
    out += "\n" + str(len(seqs)) + " sequences\n"
    #make sure group has enough sequences before continuing
    #run locarna-p on the at most 50 most abundant sequences in the group
    aln, struct = run_locarnap(seqs, 50, cpus=2, foldless=True)

    #create output folder for group
    mkdir(currotufolder)
    if(aln.getNumSeqs() < 50):
        out += str(aln.getNumSeqs()) + " unique sequences\n"
        fout = open(currotufolder + "/unique.fasta", 'w')
        fout.write(aln.toFasta())
        fout.close()
    else:
        s, h = remove_duplicates(seqs)
        out += str(len(s)) + " unique sequences\n"
        write_fasta_list(s, currotufolder + "/unique.fasta")
    out += "Structure: " + struct + "\n"

    #write out alignment and structure in fasta and stockholm formats
    #write that shit
    logout = open(currotufolder + "/log.txt", 'w')
    logout.write(out)
    logout.close()
    alnout = open(currotufolder + "/locarnap-aln.fasta", 'w')
    alnout.write(aln.toFasta() + "\n>SS_struct\n" + struct + "\n")
    alnout.close()
    alnout = open(currotufolder + "/locarnap-aln.sto", 'w')
    struct_dict = {'SS_cons': struct}
    alnout.write(stockholm_from_alignment(aln, GC_annotation=struct_dict))
    alnout.close()
    #make R2R secondary structure for alignment
    make_r2r(currotufolder + "/locarnap-aln.sto", currotufolder, "group_" + str(currgroup))


def score_rnaforester(struct1, struct2):
    '''returns relative score of two structures'''
    #return gigantically negative number if no structure for one struct
    if "(" not in struct1 or "(" not in struct2:
        raise ValueError(struct1 + "\n" + struct2 + "\nNo pairing in given structures!")
    p = Popen(["RNAforester", "--score", "-r"], stdin=PIPE, stdout=PIPE)
    p.stdin.write(''.join([struct1, "\n", struct2, "\n&"]))
    return float(p.communicate()[0].split("\n")[-2])


def group_by_forester(structgroups, foresterscore):
    topop = set([])
    for currstruct in structgroups:  # for each structure
        if currstruct in topop:  # skip if we already grouped this one
            continue
        for teststruct in structgroups:  # else compare it to all other structures
            #only compare if not equal and not already grouped
            if teststruct != currstruct and teststruct not in topop:
                score = score_rnaforester(currstruct, teststruct)
                if score > foresterscore:
                    #add matched groups clusters to current group
                    #then wipe test group to save memory and add to topop
                    structgroups[currstruct].extend(structgroups[teststruct])
                    structgroups[teststruct] = 0
                    topop.add(teststruct)
    for key in topop:  # pop all simmilar structures found
        structgroups.pop(key)
    return structgroups


def make_r2r(insto, outfolder, group):
    '''generates R2R secondary structure pdf with default colorings'''
    p = Popen(["r2r", "--GSC-weighted-consensus", insto, outfolder + "/" + group + ".sto",
        "3", "0.97", "0.9", "0.75", "4", "0.97", "0.9", "0.75", "0.5", "0.1"])
    p.wait()
    p = Popen(["r2r", outfolder + "/" + group + ".sto", outfolder + "/" + group + ".pdf"], stdout=PIPE)
    retcode = p.wait()
    #fix known r2r base-pair issue if PDF not created
    if retcode != 0:
        fin = open(outfolder + "/" + group + ".sto", 'U')
        sto = fin.readlines()
        fin.close()
        sto[-2] = "#=GF R2R SetDrawingParam autoBreakPairs true\n"
        sto[-1] = "//\n"
        fout = open(outfolder + "/" + group + ".sto", 'w')
        fout.write(''.join(sto))
        fout.close()
        p = Popen(["r2r", outfolder + "/" + group + ".sto", outfolder + "/" + group + ".pdf"], stdout=PIPE)
        p.wait()


def run_infernal(lock, cmfile, rnd, basefolder, outfolder, cpus=1, score=0.0, mpi=False):
    try:
        seqs = 0
        #Only search unique sequences to save time
        #check if previous run has removed some sequences, load correct file
        if exists(basefolder + "R" + str(rnd) + "/R" + str(rnd) + "-Unique-Remaining.fasta"):
            seqs = LoadSeqs(basefolder + "R" + str(rnd) + "/R" + str(rnd) + "-Unique-Remaining.fasta", moltype=RNA, aligned=False)
        else:
            seqs = LoadSeqs(basefolder + "R" + str(rnd) + "/R" + str(rnd) + "-Unique.fasta", moltype=RNA, aligned=False)
        params = {'--mid': True, '--Fmid': 0.0002, '--notrunc': True, '--toponly': True, '--cpu': cpus}  # '-g': True,
        if mpi:
            params['mpi'] = True
        result = cmsearch_from_file(cmfile, seqs, RNA, cutoff=score, params=params)
        fout = open(outfolder + "/R" + str(rnd) + "hits.txt", 'w')
        fout.write(str(len(result)) + " hits\nheader,bitscore,e-value\n")
        for hit in result:
            fout.write(hit[0] + "," + str(hit[14]) + "," + str(hit[15]) + "\n")
        fout.close()
        print "Round " + str(rnd) + ": " + str(len(result)) + " hits"
        sys.stdout.flush()
        lock.acquire()
        fout = open(outfolder + "/log.txt", 'a')
        fout.write("Round " + str(rnd) + ": " + str(len(result)) + " hits\n")
        fout.close()
        lock.release()
    except Exception, e:
        print str(e)
        lock.release()
    finally:
        lock.release()


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
    parser.add_argument('--inf', action="store_true", default=False, help="Skip straight to infernal \
    (Default False)")
    parser.add_argument('--mpi', action="store_true", default=False, help="Run infernal in MPI mode \
    (Default False)")

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
    print "Running locarna-p over " + str(len(otus)) + " clusters"
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
    nostructout = open(otufolder + "unstructured_clusters.txt", 'w')
    for cluster, struct in MinimalFastaParser(cfo):
        #remove unstructured groups and print to file for examination
        if "(" not in struct:
            print cluster + ": NO STRUCTURE PREDICTED"
            nostructout.write(cluster + "\t" + struct + "\n")
            continue
        if struct in structgroups:
            structgroups[struct].append(cluster)
        else:
            structgroups[struct] = [cluster]
    cfo.close()
    nostructout.close()
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
        gout = open(otufolder + "groups.txt", 'w')
        gout.write(str(foresterscore) + "\n")
        for group in structgroups:
            gout.write(group + ":")
            for g in structgroups[group]:
                gout.write(g+" ")
            gout.write("\n")
        gout.close()
        iteration += 1

        while startcount != endcount:  # keep refining while we are still grouping structs
            startcount = len(structgroups)
            print "iteration " + str(iteration) + ": " + str(len(structgroups)) + " initial groups"
            #Refold all the groups to get new consensus secondary structure
            #make a pool of workers, one for each cpu available
            manager = Manager()
            hold = manager.dict()
            pool = Pool(processes=args.c)
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

        print str(len(structgroups)) + " final groups (" + str((time() - secs) / 60) + "m)"
        gout = open(otufolder + "groups.txt", 'a')
        gout.write("FINAL\n")
        gout.close()
    else:
        print str(len(structgroups)) + " final groups"

    print "==Creating CM and r2r structures=="
    secs = time()
    pool = Pool(processes=int(ceil(args.c/2)))
    #run the pool over all clusters to get file of structures
    for group, currstruct in enumerate(structgroups):
        pool.apply_async(func=run_locarnap_for_infernal, args=(group+1, structgroups[currstruct], otus, otufolder))
    pool.close()
    pool.join()
    print "Runtime: " + str((time() - secs) / 60) + "m"
    print "==Running Infernal for all groups=="
    print "Infernal score cutoff: " + str(infernalscore)
    #create the csv file for holding all the hit counts
    ihits = open(otufolder + "infernalhits.csv", 'w')
    ihits.write(",")
    for i in range(1, args.r+1):
        ihits.write("round" + str(i) + ",")
    ihits.write("\n")
    #loop over each group and run infernal on it for all rounds
    for group in walk(otufolder).next()[1]:
        secs = time()
        print "GROUP: " + group
        logfile = open(otufolder + group + "/log.txt")
        log = logfile.readlines()
        logfile.close()
        print ''.join(log)
        seqs = int(log[1].split()[0])
        skip = False
        if exists(otufolder + group + "/R1hits.txt"):
            skip = True
        #only run infernal if there were more than 100 total sequences in group
        if seqs > 99 and not skip:
            currotufolder = otufolder + group
            #create the cm file and calibrate it
            cmfile = open(currotufolder + "/infernal_" + group + ".cm", 'w')
            cmfile.write(cmbuild_from_file(currotufolder + "/locarnap-aln.sto"))
            cmfile.close()
            cmfile = currotufolder + "/infernal_" + group + ".cm"
            calibrate_file(cmfile)

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
            for i in range(args.r, 0, -1):
                #run cmsearch over every round of SELEX
                #if there is extra cpu power, apply it!
                if extracpus > 0:
                    pool.apply_async(func=run_infernal, args=(lock, cmfile, i, args.f, currotufolder, procs+1, infernalscore, args.mpi))
                    extracpus -= 1
                else:
                    pool.apply_async(func=run_infernal, args=(lock, cmfile, i, args.f, currotufolder, procs, infernalscore, args.mpi))
            pool.close()
            pool.join()
            logfile = open(otufolder + group + "/log.txt")
            roundhits = logfile.readlines()[4:]
            logfile.close()
            roundhits.sort()
            hitscsv = open(otufolder + "infernalhits.csv", 'a')
            hitscsv.write(group + ",")
            for r in roundhits:
                hitscsv.write(r.split()[2] + ",")
            hitscsv.write("\n")
            print "Runtime: " + str((time() - secs) / 60) + "m"
        elif skip:
            print "Group already run"
            logfile = open(otufolder + group + "/log.txt")
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
print "Program ended ", datetime.now(), "   Runtime: " + str((time() - starttime)/3600) + "h"
