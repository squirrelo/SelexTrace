from sys import stdout
from os.path import exists
from os import mkdir, walk
from cogent.app.muscle_v38 import align_unaligned_seqs
from cogent.app.infernal_v11 import cmbuild_from_file, calibrate_file
from cogent.parse.fasta import MinimalFastaParser
from cogent import LoadSeqs, RNA
from time import time
from datetime import datetime
from multiprocessing import Pool, Manager
import argparse
from math import floor
from selextrace.stutils import cluster_seqs, write_fasta_list
from selextrace.bayeswrapper import bayesfold
from selextrace.ctilib import fold_clusters, group_by_shape, \
    run_fold_for_infernal, build_reference, group_to_reference, group_denovo, \
    group_by_forester, group_by_distance, make_r2r, run_infernal


if __name__ == "__main__":
    starttime = time()
    #TURN INTO COGENT OPTION PARSING   cogent.util.option_parsing
    parser = argparse.ArgumentParser(description="Runs sequence grouping \
        and infernal over all rounds of a SELEX selection")
    parser.add_argument('-i', required=True, 
        help="FASTA file of unique sequences sorted by abundance.")
    parser.add_argument('-f', required=True, 
        help="Base folder holding FASTA files for all rounds of selection")
    parser.add_argument('-o', required=True, 
        help="Base folder to output all data")
    parser.add_argument('-r', required=True, type=int, 
        help="Current round of selection clusters come from")
    parser.add_argument('--sim', required=True, type=float, default=0.99,
        help="Simmilarity for uclust. (Default 0.99)")
    parser.add_argument('--minseqs', type=int, default=100,
        help="Min number of seqs for group to be significant (Default 100)")
    parser.add_argument('--rsc', type=float, default=40, 
        help="Score cutoff for RNAdistance (Default 40)")
    parser.add_argument('--fsc', type=int, default=200,
        help="Score cutoff for RNAforester. (Default 200)")
    parser.add_argument('--isc', type=float, default=0.0,
        help="Score cutoff for Infernal. (Default 0.0)")
    parser.add_argument('-c', type=int, default=1,
        help="Number of CPUs to use (Default 1)")

    args = parser.parse_args()
    if args.r < 1:
        print "ERROR: round must be at least 1!"
        exit(1)
    if args.c < 1:
        print "ERROR: CPU count must be at least 1!"
        exit(1)
    if args.fsc < 0:
        print "ERROR: RNAforester score cutoff must be greater than 0!"
        exit(1)
    if args.sim < 0.0 or args.sim > 1.0:
        print "ERROR: Infernal score cutoff must be greater than 0!"
        exit(1)

    structscore = args.rsc
    infernalscore = args.isc
    otufolder = args.o

    if args.f[:-1] != "/":
        args.f += "/"

    if otufolder[:-1] != "/":
        otufolder += "/"
    if not exists(otufolder):
        mkdir(otufolder)

    date = str(datetime.now())
    print "Program started ", date

    #print out run info to a file
    infofile = open(otufolder + "runparams.txt", 'w')
    infofile.write(''.join(["Program started ", date, "\n",
                    "FASTA file:\t", args.i, "\n",
                    "Base folder:\t", args.f, "\n",
                    "Output folder:\t", args.o, "\n"
                    "Selection round:\t", str(args.r), "\n",
                    "Uclust simmilarity:\t", str(args.sim), "\n",
                    "Min seqs for group:\t", str(args.minseqs), "\n",
                    "RNAdistance score cutoff:\t", str(args.rsc), "\n",
                    "RNAforester score cutoff:\t", str(args.fsc), "\n",
                    "Infernal min score:\t", str(args.isc), "\n",
                    "CPUs:\t", str(args.c), "\n"]))
    infofile.close()
    print "==Clustering sequences by primary sequence=="
    clusters = {}
    secs = time()
    if exists(otufolder + "clusters.txt"):
        clustersin = open(otufolder + "clusters.txt")
        currclust = ""
        for header, seq in MinimalFastaParser(clustersin):
            if "cluster_" in header:
                currclust = header
                clusters[currclust] = []
            else:
                clusters[currclust].append((header, seq))
        clustersin.close()
        print "Sequences previously clustered,", len(clusters), "clusters"
    else:
        print "Running uclust over sequences"
        #cluster the initial sequences by sequence simmilarity
        clusters = cluster_seqs(args.i, args.sim, folderout=args.o,
            gapopen='10.0', gapext='10.0')

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
            pool.apply_async(func=fold_clusters, args=(lock, cluster,
                clusters[cluster], otufolder))
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
        raise AssertionError(str(count) + " structures, " + str(len(clusters))
            + " clusters. Not all clusters folded!")
    else:
        print "Clusters previously folded"

    clusters.clear()
    del clusters
    if not exists(otufolder + "clusters_aln.fasta"):
        startaln = time()
        #add refseq and align all sequences for later seq/struct comparisons
        refin = open(args.f + "refseq.fasta")
        crap, refseq = MinimalFastaParser(refin).next()
        refin.close()
        for struct in structgroups:
            structgroups[struct].append(("refseq", refseq))
            structgroups[struct] = align_unaligned_seqs(structgroups[struct], RNA)
            #hacky way to make sure refseq first, will need to fix later
            structgroups[struct].Names.remove("refseq")
            structgroups[struct].Names.insert(0, "refseq")
        #write that shit to a file to save time if rerun needed
        cout = open(otufolder + "clusters_aln.fasta", 'w')
        for struct in structgroups:
            cout.write(">%s\n%s\n" % ("newcluster", struct))
            cout.write(structgroups[struct].toFasta() + "\n")
        cout.close()
        print "Aligned all clusters: " + str((time() - startaln)/60) + " min"
    else:
        cin = open(otufolder + "clusters_aln.fasta")
        currclust = []
        currstruct = ""
        structgroups = {}
        first = True
        for header, seq in MinimalFastaParser(cin):
            if "cluster" in header:
                if not first:
                    structgroups[currstruct] = LoadSeqs(data=currclust, moltype=RNA, aligned=True)
                    structgroups[currstruct].Names.remove("refseq")
                    structgroups[currstruct].Names.insert(0, "refseq")
                first = False
                currstruct = seq
                currclust = []
            else:
                currclust.append((header, seq))
        print "Clusters previously aligned"

    print "Runtime: " + str((time() - secs)/60) + " min"
    print "==Grouping clusters by secondary structure=="
    #GROUP THE SECONDARY STRUCTURES BY RNADISTANCE
    secs = time()
    skipiter = False
    if exists(otufolder + "fasta_groups/"):
        skipiter = True

    if not skipiter:
        print "Abstract shape assignment"
        secs = time()
        manager = Manager()
        groups_shape = manager.dict()
        pool = Pool(processes=args.c)
        #run the pool over all groups to get structures
        for struct in structgroups:
            pool.apply_async(func=group_by_shape, args=(groups_shape, struct))
        pool.close()
        pool.join()
        print len(groups_shape), "shape groups ("+str((time()-secs)/60)+" min)"

        print "start: " + str(len(structgroups)) + " initial groups"
        #initial clustering by structures generated in first folding,
        #broken out by shapes. No comparison needed at first if not same shape,
        #as most likely not simmilar enough
        hold = {}
        pool = Pool(processes=args.c)
        #run the pool over all shape groups to get final grouped structgroups
        fout = open(otufolder + "shapesizes.txt", 'w')
        for shapegroup in groups_shape.keys():
            #create dict of structs in the group, then call grouping func on it
            groupinfo = {struct: structgroups[struct] for struct in groups_shape[shapegroup]}
            fout.write(shapegroup + "\t" + str(len(groupinfo)) + "\n")
            pool.apply_async(func=group_by_distance,
                args=(groupinfo, structscore), callback=hold.update)
            #hold.update(group_by_distance(groupinfo, structscore))
        fout.close()
        pool.close()
        pool.join()
        #memory saving wipe of structgroups, groups_shape, and groupinfo
        groups_shape.clear()
        del groups_shape
        groupinfo.clear()
        del groupinfo
        structgroups.clear()
        del structgroups
        #hold should now be the combined dictionaries from all calls of
        #group_by_forester, aka new structgroups
        #do one more grouping with all remaining structs regardless of shape
        structgroups = dict(hold)
        del hold
        structgroups = group_by_distance(structgroups, structscore)

        #sort all structure sequences by count, highest to lowest
        for struct in structgroups:
            structgroups[struct].Names.remove('refseq')
            structgroups[struct].Names.sort(reverse=True,
                key=lambda count: int(count.split('_')[1]))
        print str(len(structgroups))+" end groups ("+str((time()-secs)/60)+" min)"

        #write out fasta files for groups: header of each sequence in the group
        mkdir(otufolder+"fasta_groups")
        for num, group in enumerate(structgroups):
            structgroups[group].degap().writeToFile(otufolder+"fasta_groups/group_"+str(num)+".fasta", "fasta")

    else:
        print "Previously grouped"

    #wipe out structgroups dict to save memory
    structgroups.clear()
    del structgroups

    print "==Creating CM and r2r structures=="
    secs = time()
    #smaller pool for memeory savings
    #NEED TO FIX BAYESFOLD ARRAY2D BEING A MEMORY HOG
    pool = Pool(processes=int(args.c/2))
    #run the pool over all groups to get final structures
    for group in walk(otufolder + "fasta_groups").next()[2]:
        groupnum = group.split("_")[-1].split(".")[0]
        pool.apply_async(func=run_fold_for_infernal,
        args=(groupnum, otufolder+"fasta_groups/"+group, otufolder, args.minseqs))
    pool.close()
    pool.join()

    #get sequence counts for each group
    infernalorder = []
    groups = {}
    count = 0
    for group in walk(otufolder).next()[1]:
        if group == "fasta_groups":
            continue
        count += 1
        #read in group sequence counts
        log = open(otufolder + group + "/log.txt")
        loginfo = log.readlines()
        log.close()
        infernalorder.append((group, int(loginfo[1].split()[0]),
            int(loginfo[2].split()[0])))
        #read in group structure and build dict for families creation
        structin = open(otufolder + group + "/bayesfold-aln.fasta")
        struct = structin.read().split("\n")[-2].strip()
        structin.close()
        groups[struct] = [group]

    #write out file of sequence counts
    infernalorder.sort(reverse=True, key=lambda x: x[1])
    groupsizefile = open(otufolder + "/group_sizes.txt", 'w')
    groupsizefile.write("Group\tTotal Seqs\tUnique Seqs\n")
    for info in infernalorder:
        groupsizefile.write("%s\t%s\t%s\n" % info)
    groupsizefile.close()

    print count, "final groups"
    print "Runtime: " + str((time() - secs) / 3600) + " hrs"

    print "==Creating families from groups=="
    print "RNAforester score cutoff:", args.fsc
    print len(groups), "initial groups"
    #Now group by forester local aignment to get larger families grouped
    groups = group_by_forester(groups, args.fsc, args.c)
    fout = open(otufolder + "families.txt", 'w')
    for fam in groups:
        fout.write("\t".join(groups[fam]) + "\n")
    fout.close()
    print len(groups), "final groups"
    print "Runtime: " + str((time() - secs) / 3600) + " hrs"
    groups.clear()
    del groups

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
    inftime = 0.0
    secs = time()
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
            cmfile.write(cmbuild_from_file(currotufolder+"/bayesfold-aln.sto"))
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
                    pool.apply_async(func=run_infernal, args=(lock, cmfile, i,
                        args.f, currotufolder, procs+1, infernalscore))
                    extracpus -= 1
                else:
                    pool.apply_async(func=run_infernal, args=(lock, cmfile, i,
                        args.f, currotufolder, procs, infernalscore))
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
            inftime += (time() - secs)
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
    print "Runtime:", (time()-secs)/3600, " hrs   Average runtime: ", \
        (inftime/count)/60, " mins"
    print "Program ended ", datetime.now(), "   Runtime: " + \
        str((time() - starttime)/3600) + "h"
