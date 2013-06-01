#dupestest.py /path/to/OTU.fasta /path/to/folder/out CPUs
from sys import stdout
from os.path import exists
from os import mkdir, walk
from cogent import LoadSeqs, RNA
from cogent.app.infernal_v11 import cmsearch_from_file, cmbuild_from_file, calibrate_file
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from subprocess import Popen, PIPE
from time import time
from datetime import datetime
from multiprocessing import Pool, Manager
import argparse
from math import floor, ceil, log
from random import shuffle
from stutils import cluster_seqs
from bayeswrapper import bayesfold

def fold_clusters(lock, cluster, seqs, otufolder):
    '''Function for multithreading.
    Computes structure for a cluster and writes it to file'''
    try:
        aln, struct = bayesfold(seqs)
        #write structure out to file
        lock.acquire()
        cfo = open(otufolder + "cluster_structs.fasta", 'a')
        cfo.write(">" + cluster + "\n" + struct + "\n")
        cfo.close()
        lock.release()
        #print cluster + ": " + struct
        #stdout.flush()
    except Exception, e:
        cluster, struct, "\nERROR!"
        stdout.flush()
        lock.release()
    finally:
        lock.release()


def group_by_shape(structs):
    fams = {}
    for struct in structs:
        famed = False
        #convert to shape
        gshape = get_shape(struct)
        for secshape in fams:
            #loop over all previously found shapes, see if it fits
            if gshape == secshape:
                fams[gshape].append(struct)
                famed = True
                break
        #if not fitted, create new group for this shape
        if not famed:
            fams[gshape] = [struct]
    return fams


def fold_groups(seqs, struct, hold, numkept=None):
    '''Function for multithreading.
    Recomputes structure for a group'''
    try:
        aln, currstruct = bayesfold(seqs)
        if currstruct == "":
            currstruct = struct
        hold[currstruct] = structgroups[struct]
    except Exception, e:
        print str(e)
        stdout.flush()


def run_fold_for_infernal(currgroup, groupfasta, basefolder, minseqs=1):
    '''Function for multithreading. Creates the final BayesFold alignment and 
    writes to files, then r2r struct'''
    try:
        #run locana-p on the superclusters to get alignment and structure
        #skip if already run and program just crashed or whatever
        currotufolder = basefolder + "group_" + str(currgroup)
        if exists(currotufolder):
            return ""
        seqs = []
        count = 0
        out = "group " + str(currgroup) + ": "
        for header, seq in MinimalFastaParser(open(groupfasta, 'rU')):
            seqs.append((header.split()[0] + "_" + header.split("_")[1], seq))
            count += int(header.split("_")[1])
        out += "\n" + str(count) + " sequences\n"
        if count < minseqs:
            return ""
        #run BayesFold on sequences in the group
        aln, struct = bayesfold(seqs)
        #create output folder for group
        mkdir(currotufolder)
        out += str(aln.getNumSeqs()) + " unique sequences\n"
        out += "Structure: " + struct + "\n"
        #write out alignment and structure in fasta and stockholm formats
        #write that shit
        logout = open(currotufolder + "/log.txt", 'w')
        logout.write(out)
        logout.close()
        alnout = open(currotufolder + "/bayesfold-aln.fasta", 'w')
        alnout.write(aln.toFasta() + "\n>SS_struct\n" + struct + "\n")
        alnout.close()
        alnout = open(currotufolder + "/bayesfold-aln.sto", 'w')
        struct_dict = {'SS_cons': struct}
        alnout.write(stockholm_from_alignment(aln, GC_annotation=struct_dict))
        alnout.close()
        #make R2R secondary structure for alignment
        make_r2r(currotufolder + "/bayesfold-aln.sto", currotufolder, "group_" + str(currgroup))
    except Exception, e:
        print str(e)
        stdout.flush()


def score_rnaforester(struct1, struct2):
    '''returns relative score of two structures'''
    #return gigantically negative number if no structure for one struct
    if "(" not in struct1 or "(" not in struct2:
        raise ValueError(struct1 + "\n" + struct2 + "\nNo pairing in structs!")
    p = Popen(["RNAforester", "--score", "-r"], stdin=PIPE, stdout=PIPE)
    p.stdin.write(''.join([struct1, "\n", struct2, "\n&"]))
    return float(p.communicate()[0].split("\n")[-2])


def build_reference(dictkeys, refsize):
    '''Creates a random list of references comprising percent of total list'''
    shuffle(dictkeys)
    #needed in case passed a float
    refsize = int(refsize)
    #return reference, nonreference by slicing list
    return dictkeys[:refsize], dictkeys[refsize:]


def group_to_reference(fulldict, reference, nonref, foresterscore):
    nogroup = []
    for currstruct in nonref:
        score = foresterscore
        bestref = ""
        for teststruct in reference:
            holdscore = score_rnaforester(currstruct, teststruct)
            if holdscore > score:
                score = holdscore
                bestref = teststruct
        if bestref != "":
            fulldict[bestref].extend(fulldict[currstruct])
            fulldict.pop(currstruct)
        else:
            nogroup.append(currstruct)
    return fulldict, nogroup


def group_denovo(fulldict, keys, foresterscore):
    topop = []
    for pos, currstruct in enumerate(keys):
        score = foresterscore
        bestref = ""
        for secpos in range(pos+1, len(keys)):
            holdscore = score_rnaforester(currstruct, keys[secpos])
            if holdscore > score:
                score = holdscore
                bestref = keys[secpos]
        if bestref != "":
            fulldict[bestref].extend(fulldict[currstruct])
            fulldict.pop(currstruct)
            topop.append(pos)
    topop.sort(reverse=True)
    for pos in topop:
        keys.pop(pos)
    return fulldict, keys


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
    #TURN INTO COGENT OPTION PARSING   cogent.util.option_parsing
    parser = argparse.ArgumentParser(description="Runs sequence clustering \
    and infernal over all rounds of a SELEX selection")
    parser.add_argument('-i', required=True, help="FASTA file of unique sequences sorted by abundance.")
    parser.add_argument('-f', required=True, help="Base folder holding FASTA files for all \
    unique sequences in the rounds of selection (for infernal)")
    parser.add_argument('-o', required=True, help="Base folder to output all data to")
    parser.add_argument('-r', required=True, type=int, help="Current round of selection\
    clusters come from")
    parser.add_argument('--sim', required=True, type=float, default=0.99, help="Simmilarity for uclust.\
    (Default 0.99)")
    parser.add_argument('--minseqs', type=int, default=100, help="Minimum number of sequences for \
    groups to be significant (Default 100)")
    parser.add_argument('--iter', type=int, default=10, help="Max number of interations during \
    RNAForester grouping (Default 10)")
    parser.add_argument('--rsc', type=float, default=0.5, help="Score cutoff for rnaforester \
    (Default 0.5)")
    parser.add_argument('--isc', type=float, default=0.0, help="Score cutoff for infernal.\
    (Default 0.0)")
    parser.add_argument('-c', type=int, default=1, help="Number of CPUs to use \
    (Default 1)")

    args = parser.parse_args()
    if args.r < 1:
        print "ERROR: round must be at least 1!"
        exit(1)
    if args.c < 1:
        print "ERROR: CPU count must be at least 1!"
        exit(1)
    if args.isc < 0.0:
        print "ERROR: Infernal score cutoff must be greater than 0!"
        exit(1)
    if args.sim < 0.0 or args.sim > 1.0:
        print "ERROR: Infernal score cutoff must be greater than 0!"
        exit(1)
    if args.iter < 1:
        print "ERROR: Must have at least one iteration!"
        exit(1)

    foresterscore = args.rsc
    infernalscore = args.isc
    otufolder = args.o

    if args.f[:-1] != "/":
        args.f += "/"

    if otufolder[:-1] != "/":
        otufolder += "/"
    if not exists(otufolder):
        mkdir(otufolder)

    print "Program started ", datetime.now()

    #print out run info to a file
    infofile = open(otufolder + "runparams.txt", 'w')
    infofile.write(''.join(["FASTA file:\t\t\t\t", args.i, "\n",
                    "Base folder:\t\t\t", args.f, "\n",
                    "Output folder:\t\t\t", args.o, "\n"
                    "Selection round:\t\t", str(args.r), "\n",
                    "Uclust simmilarity:\t\t", str(args.sim), "\n",
                    "Min seqs for group:\t\t", str(args.minseqs), "\n",
                    "Iterations for grouping:", str(args.iter), "\n",
                    "RNAForester min score:\t", str(args.rsc), "\n",
                    "Infernal min score:\t\t", str(args.isc), "\n",
                    "CPUs:\t\t\t\t\t", str(args.c), "\n"]))
    infofile.close()
    print "==Clustering sequences by primary sequence=="
    clusters = {}
    secs = time()
    if exists(otufolder + "clusters.txt"):
        print "sequences previously clustered"
        clustersin = open(otufolder + "clusters.txt")
        currclust = ""
        for header, seq in MinimalFastaParser(clustersin):
            if "cluster_" in header:
                currclust = header
                clusters[currclust] = []
            else:
                clusters[currclust].append((header, seq))
    else:
        print "Running uclust over sequences"
        #cluster the initial sequences by sequence simmilarity
        clusters = cluster_seqs(args.i, args.sim, folderout=args.o, gapopen='10.0', gapext='10.0')

        #print that shit to file
        cout = open(otufolder + "clusters.txt", 'w')
        hold = clusters.keys()
        hold.sort()
        for cluster in hold:
            cout.write(">%s\n%s\n" % (cluster, cluster))
            for seq in clusters[cluster]:
                cout.write(">%s\n%s\n" % seq)
        cout.close()
        print str(len(clusters)) + " clusters"
        print "Runtime: " + str((time() - secs)/60) + " min"

    if not exists(otufolder + "cluster_structs.fasta"):
        #create file to write to if not already there
        cfo = open(otufolder + "cluster_structs.fasta", 'w')
        cfo.close()

        print "Running BayesFold over " + str(len(clusters)) + " clusters"
        secs = time()
        #make a pool of workers, one for each cpu available
        manager = Manager()
        pool = Pool(processes=args.c)
        lock = manager.Lock()
        #run the pool over all clusters to get file of structures
        for cluster in clusters:
            pool.apply_async(func=fold_clusters, args=(lock, cluster, clusters[cluster], otufolder))
        pool.close()
        pool.join()

    #read in all structures now that they are folded
    structgroups = {}
    count = 0
    cfo = open(otufolder + "cluster_structs.fasta", 'rU')
    nostruct = []
    for cluster, struct in MinimalFastaParser(cfo):
        count += 1
        #remove unstructured groups
        if struct.count("(") < 5:
            nostruct.append(cluster)
            continue
        #create dictionary of sequences now keyed to structures
        if struct in structgroups:
            structgroups[struct].extend(clusters[cluster])
        else:
            structgroups[struct] = clusters[cluster]
    cfo.close()

    if count != len(clusters):
        raise AssertionError(str(count) + " structures, " + str(len(clusters)) + " clusters. Not all clusters folded!")

    print "Runtime: " + str((time() - secs)/60) + " min"
    print "==Grouping clusters by secondary structure=="
    #GROUP THE SECONDARY STRUCTURES BY RNAFORESTER
    secs = time()
    skipiter = False
    if exists(otufolder + "fasta_groups/"):
        skipiter = True

    if not skipiter:
        print "RNAforester score threshold: " + str(foresterscore)
        #Now need to iteratively refine the groups down
        iteration = 0
        secs = time()
        #creating reference subset, then de-novo grouping the references.
        #then group to those references. Should speed things up significantly. 
        print "start: " + str(len(structgroups)) + " initial groups"
        finishlen = len(structgroups) * 0.01
        reference, nonreference = build_reference(structgroups.keys(), finishlen)
        print len(reference), "KEYS!", str((time() - secs) / 60), " min"
        structgroups, reference = group_denovo(structgroups, reference, foresterscore)
        print len(reference), "GROUPED REF!", str((time() - secs) / 3600), " hrs"
        structgroups, ungrouped = group_to_reference(structgroups, reference, nonreference, foresterscore)
        print len(ungrouped), "UNGROUPED!", str((time() - secs) / 3600), " hrs"
        # keep refining while not at limit and are still grouping structs
        startungrouped = 0
        endungrouped = 1
        while len(ungrouped) > finishlen and startungrouped != endungrouped:
            startungrouped = len(ungrouped)
            iteration += 1
            print "iteration " + str(iteration) + ": " + str(len(structgroups)) + " initial groups"
            reference, nonreference = build_reference(ungrouped, finishlen)
            structgroups, reference = group_denovo(structgroups, reference, foresterscore)
            structgroups, ungrouped = group_to_reference(structgroups, reference, nonreference, foresterscore)
            print len(ungrouped), "UNGROUPED"
            endungrouped = len(ungrouped)
        #end while
        structgroups, reference = group_denovo(structgroups, ungrouped, foresterscore)
        #sort all structure sequences by count
        for struct in structgroups:
            structgroups[struct].sort(reverse=True, key=lambda count: int(count[0].split('_')[1]))
        print str(len(structgroups)) + " final groups (" + str((time() - secs) / 3600) + " hrs)"

        #write out fasta files for groups: header of each sequence in the group
        mkdir(otufolder+"fasta_groups")
        for num, group in enumerate(structgroups):
            gout = open(otufolder+"fasta_groups/group_" + str(num) + ".fasta", 'w')
            for g in structgroups[group]:
                gout.write(">%s\n%s\n" % g)
            gout.close()
    else:
        print "Previously grouped"

    #wipe out structgroups dict to save memory
    structgroups.clear()
    structgroups = 0

    print "==Creating CM and r2r structures=="
    secs = time()
    pool = Pool(processes=int(ceil(args.c/2)))
    #run the pool over all groups to get structures
    for group in walk(otufolder + "fasta_groups").next()[2]:
        groupnum = group.split("_")[-1].split(".")[0]
        pool.apply_async(func=run_fold_for_infernal, args=(groupnum, otufolder+"fasta_groups/" + group, otufolder, args.minseqs))
    pool.close()
    pool.join()
    print "Runtime: " + str((time() - secs) / 3600) + " hrs"

    #get sequence counts for each group
    infernalorder = []
    count = 0
    for group in walk(otufolder).next()[1]:
        if group == "fasta_groups":
            continue
        count += 1
        log = open(otufolder + group + "/log.txt")
        loginfo = log.readlines()
        log.close()
        infernalorder.append((group, int(loginfo[1].split()[0]), int(loginfo[2].split()[0])))
    #write out file of sequence counts
    infernalorder.sort(reverse=True, key=lambda x: x[1])
    groupsizefile = open(otufolder + "/group_sizes.txt", 'w')
    groupsizefile.write("Group\tTotal Seqs\tUnique Seqs\n")
    for info in infernalorder:
        groupsizefile.write("%s\t%s\t%s\n" % info)
    groupsizefile.close()

    print count, "final groups"
    print "==Running Infernal for all groups=="
    print "Infernal score cutoff: " + str(infernalscore)
    #create the csv file for holding all the hit counts
    ihits = open(otufolder + "infernalhits.csv", 'w')
    ihits.write(",")
    for i in range(1, args.r+1):
        ihits.write("round" + str(i) + ",")
    ihits.write("\n")
    ihits.close()
    #loop over each group and run infernal on it for all rounds
    for groupinfo in infernalorder:
        group = groupinfo[0]
        if group == "fasta_groups":
            continue
        secs = time()
        skip = False
        if exists(otufolder + group + "/R1hits.txt"):
            skip = True
        #only run infernal if there were more than 100 total sequences in group
        if not skip:
            currotufolder = otufolder + group
            #create the cm file and calibrate it
            cmfile = open(currotufolder + "/infernal_" + group + ".cm", 'w')
            cmfile.write(cmbuild_from_file(currotufolder + "/bayesfold-aln.sto"))
            cmfile.close()
            cmfile = currotufolder + "/infernal_" + group + ".cm"
            calibrate_file(cmfile, cpus=args.c)

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
                    pool.apply_async(func=run_infernal, args=(lock, cmfile, i, args.f, currotufolder, procs+1, infernalscore))
                    extracpus -= 1
                else:
                    pool.apply_async(func=run_infernal, args=(lock, cmfile, i, args.f, currotufolder, procs, infernalscore))
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
            print group + "\tRuntime: " + str((time() - secs) / 60) + " min"
        else:
            print group + "\talready run"
            logfile = open(otufolder + group + "/log.txt")
            roundhits = logfile.readlines()[4:]
            logfile.close()
            roundhits.sort()
            hitscsv = open(otufolder + "infernalhits.csv", 'a')
            hitscsv.write(group + ",")
            for r in roundhits:
                hitscsv.write(r.split()[2] + ",")
            hitscsv.write("\n")
print "Program ended ", datetime.now(), "   Runtime: " + str((time() - starttime)/3600) + "h"
