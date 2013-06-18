from sys import stdout
from os.path import exists
from os import mkdir
from cogent import LoadSeqs, RNA
from cogent.app.infernal_v11 import cmsearch_from_file
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from subprocess import Popen, PIPE
from math import ceil
from random import shuffle
from bayeswrapper import bayesfold
from selextrace.stutils import get_shape


def fold_clusters(lock, cluster, seqs, otufolder):
    '''Function for multithreading.
    Computes structure for a cluster and writes it to file'''
    aln, struct = bayesfold(seqs, params={"-diags": True})
    #write structure out to file
    lock.acquire()
    cfo = open(otufolder + "cluster_structs.fasta", 'a')
    cfo.write(">" + cluster + "\n" + struct + "\n")
    cfo.close()
    #print cluster + ": " + struct
    #stdout.flush()
    lock.release()


def group_by_shape(shapedict, struct):
    '''Multithreading function
    Takes in shapedict manager dictionary and structure
    Fits structure's shape into dictionary or adds new shape'''
    try:
        famed = False
        #convert to shape
        gshape = get_shape(struct)
        for famshape in shapedict.keys():
            #loop over all previously found shapes, see if it fits
            if gshape == famshape:
                shapedict[gshape] += [struct]
                famed = True
                break
        #if not fitted, create new group for this shape
        if not famed:
            shapedict[gshape] = [struct]
    except Exception, e:
        print "ERROR:", str(e)
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
        print "group " + str(currgroup) + ":\t" + str(len(seqs)) + "\t" + str(count)
        stdout.flush()
        #hard limit of 3000 sequences to align and fold for memory reasons
        if len(seqs) > 3000:
            seqs = seqs[:3000]
        #run BayesFold on sequences in the group
        #maxiters set to 3 because should have huge amount of sequences for some groups
        aln, struct = bayesfold(seqs, params={"-diags": True, "-maxiters": 2})
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

#object wrapper so create Popen object once: saves tons of overhead
class ScoreStructures(object):
    def __init__(self):
        self.p = Popen(["RNAdistance"], stdin=PIPE, stdout=PIPE)
    def __call__(self, struct1, struct2):
        self.p.stdin.write(''.join([struct1, "\n", struct2,"\n"]))
        self.p.stdin.flush()
        self.p.stdout.flush()
        return float(self.p.stdout.readline().strip().split(":")[1])  # /min(len(struct2), len(struct1))
    #def __del__(self):
    #    self.p.terminate()


def build_reference(dictkeys, refsize):
    '''Creates a random list of references comprising percent of total list'''
    shuffle(dictkeys)
    #needed in case passed a float
    refsize = int(refsize)
    #return reference, nonreference by slicing list
    return dictkeys[:refsize], dictkeys[refsize:]


def group_to_reference(fulldict, reference, nonref, structscore):
    nogroup = []
    score_structures = ScoreStructures()
    for currstruct in nonref:
        score = structscore
        bestref = ""
        for teststruct in reference:
            holdscore = score_structures(currstruct, teststruct)
            if holdscore < score:
                score = holdscore
                bestref = teststruct
        if bestref != "":
            fulldict[bestref].extend(fulldict[currstruct])
            fulldict.pop(currstruct)
        else:
            nogroup.append(currstruct)
    return fulldict, nogroup


def group_denovo(fulldict, keys, structscore):
    topop = []
    score_structures = ScoreStructures()
    for pos, currstruct in enumerate(keys):
        score = structscore
        bestref = ""
        for secpos in range(pos+1, len(keys)):
            holdscore = score_structures(currstruct, keys[secpos])
            if holdscore < score:
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


def group_by_forester(structgroups, structscore, specstructs=None):
        '''Does grouping by way of de-novo reference creation and clustering
            structgroups - dictionary with ALL structures and the sequences that fall in them
            structscore - maximum score to consider grouping structures
            specstructs - a list of a subset of structures in structgroups to cluster
        '''
        #just de-novo group if 10 or less to save time and effort
        if len(structgroups) <= 10:
            if specstructs == None:
                structgroups, reference = group_denovo(structgroups, structgroups.keys(), structscore)
            else:
                structgroups, reference = group_denovo(structgroups, specstructs, structscore)
            return structgroups
        #for speed, get 1% as initial clustering. need at least 10 structs though
        finishlen = int(ceil(len(structgroups) * 0.01))
        if finishlen < 10:
            finishlen = 10
        #do initial ref grab by either all structures or specific ones passed
        if specstructs == None:
            reference, nonreference = build_reference(structgroups.keys(), finishlen)
        else:
            reference, nonreference = build_reference(specstructs, finishlen)
        startungrouped = len(structgroups)
        structgroups, reference = group_denovo(structgroups, reference, structscore)
        structgroups, ungrouped = group_to_reference(structgroups, reference, nonreference, structscore)
        endungrouped = len(structgroups)
        # keep refining while not at limit and are still grouping structs
        while len(ungrouped) > finishlen and startungrouped != endungrouped:
            startungrouped = len(structgroups)
            reference, nonreference = build_reference(ungrouped, finishlen)
            structgroups, reference = group_denovo(structgroups, reference, structscore)
            structgroups, ungrouped = group_to_reference(structgroups, reference, nonreference, structscore)
            endungrouped = len(structgroups)
        #end while
        structgroups, reference = group_denovo(structgroups, ungrouped, structscore)
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

