from sys import stdout
from os.path import exists
from os import mkdir
from cogent import LoadSeqs, RNA
from cogent.core.sequence import RnaSequence
from cogent.app.infernal_v11 import cmsearch_from_file
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from cogent.app.muscle_v38 import align_unaligned_seqs
from cogent.align.align import classic_align_pairwise
from subprocess import Popen, PIPE
from math import ceil
from random import shuffle
from bayeswrapper import bayesfold
from selextrace.stutils import get_shape
from multiprocessing import Pool

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


#object wrapper so create Popen object once: saves DAYS of overhead
class ScoreStructures(object):
    def __init__(self):
        self.p = Popen(["RNAdistance"], stdin=PIPE, stdout=PIPE)

    def __call__(self, struct1, struct2):
        self.p.stdin.write(''.join([struct1, "\n", struct2, "\n"]))
        self.p.stdin.flush()
        #self.p.stdout.flush()
        return float(self.p.stdout.readline().strip().split(":")[1])
    def end(self):
        self.p.kill()


def score_local_rnaforester(struct1, struct2):
    '''returns local aignment score of two structures'''
    #return gigantically negative number if no structure for one struct
    if "(" not in struct1 or "(" not in struct2:
        raise ValueError(struct1 + "\n" + struct2 + "\nNo pairing in given structures!")
    p = Popen(["RNAforester", "--score", "-l"], stdin=PIPE, stdout=PIPE)
    p.stdin.write(''.join([struct1, "\n", struct2, "\n&"]))
    return int(p.communicate()[0].split("\n")[-2])


def score_multi_forester(basestruct, checkstruct, foresterscore):
    if score_local_rnaforester(basestruct, checkstruct) > foresterscore:
        return checkstruct
    else:
        return ""

def group_by_forester(fulldict, foresterscore, cpus=1):
    structs = fulldict.keys()
    for pos, currstruct in enumerate(structs): # for each structure
        #skip if already grouped
        if currstruct not in fulldict:
            continue
        scores = set([])
        pool = Pool(processes=cpus)
        #compare everything as fast as possible using multiprocessing
        #comparisons end up as set of structs above threshold plus ""
        for teststruct in structs[pos+1:]:
            pool.apply_async(func=score_multi_forester,
                args=(currstruct, teststruct, foresterscore), callback=scores.add)
        pool.close()
        pool.join()
        #remove empty and add remaining structs to currgroup
        if "" in scores:
            scores.discard("")
        for struct in scores:
            if struct in fulldict:
                fulldict[currstruct].extend(fulldict[struct])
                fulldict.pop(struct)
    return fulldict


def build_reference(dictkeys, refsize):
    '''Creates a random list of references comprising percent of total list'''
    shuffle(dictkeys)
    #needed in case passed a float
    refsize = int(refsize)
    #return reference, nonreference by slicing list
    return dictkeys[:refsize], dictkeys[refsize:]


alnscores = {
    ('A', 'A'): 10, ('A', 'U'): -10, ('U', 'U'): 10, ('U', 'A'): -10, 
    ('C', 'A'): -10, ('C', 'U'): -10, ('G', 'G'): 10, ('U', 'C'): -10, 
    ('G', 'A'): -10, ('G', 'U'): -10, ('A', 'G'): -10, ('C', 'G'): -10, 
    ('C', 'C'): 10, ('U', 'G'): -10, ('G', 'C'): -10, ('A', 'C'): -10
}

def group_to_reference(fulldict, reference, nonref, structscore, norefseq=False):
    nogroup = []
    score_structures = ScoreStructures()
    for currstruct in nonref:
        strscore = structscore
        seqscore = 0
        bestref = ""
        #remove gaps from majority and set as seq1
        seq = fulldict[currstruct].majorityConsensus()
        seq1 = RnaSequence(''.join(seq).replace('-', ''))
        for teststruct in reference:
            holdscore = score_structures(currstruct, teststruct)
            if holdscore <= strscore:
                #remove gaps from majority and set as seq2
                seq = fulldict[teststruct].majorityConsensus()
                seq2 = RnaSequence(''.join(seq).replace('-', ''))
                #compare alignment score. subtract so lower is still better
                aln, alnscore = classic_align_pairwise(seq1, seq2, alnscores, -10, -10, False, return_score=True)
                if alnscore > seqscore:
                    strscore = holdscore
                    seqscore = alnscore
                    bestref = teststruct
        if bestref != "":
            #combine the two alignments into one alignment using reference sequence as guide
            #refseq must be ungapped to do this without realigning, hence the checks
            if norefseq:
                #realign all sequences since no refseq available
                combinedseqs = fulldict[bestref].degap().addSeqs(fulldict[currstruct].degap())
                fulldict[bestref] = align_unaligned_seqs(combinedseqs, RNA)
                continue
            #if one refseq is gapless, can easily combine them without realigning
            if not fulldict[currstruct].getGappedSeq("refseq").isGapped():
                fulldict[bestref].addFromReferenceAln(fulldict[currstruct])
            elif not fulldict[bestref].getGappedSeq("refseq").isGapped():
                fulldict[currstruct].addFromReferenceAln(fulldict[bestref])
                fulldict[bestref] = fulldict[currstruct]
            else:
                #realign all sequences since both refseqs have gaps
                #hacky but it works, need to fix later
                fulldict[bestref].Names.remove("refseq")
                combinedseqs = fulldict[bestref].degap().addSeqs(fulldict[currstruct].degap())
                fulldict[bestref] = align_unaligned_seqs(combinedseqs, RNA)
                fulldict[bestref].Names.remove("refseq")
                fulldict[bestref].Names.insert(0, "refseq")
            fulldict.pop(currstruct)
        else:
            nogroup.append(currstruct)
    score_structures.end()
    return fulldict, nogroup


def group_denovo(fulldict, keys, structscore, norefseq=False):
    topop = []
    score_structures = ScoreStructures()
    for pos, currstruct in enumerate(keys):
        strscore = structscore
        seqscore = 0
        bestref = ""
        #remove gaps from majority and set as seq1
        seq = fulldict[currstruct].majorityConsensus()
        seq1 = RnaSequence(''.join(seq).replace('-', ''))
        for secpos in range(pos+1, len(keys)):
            holdscore = score_structures(currstruct, keys[secpos])
            if holdscore <= strscore:
                #remove gaps from majority and set as seq2
                seq = fulldict[keys[secpos]].majorityConsensus()
                seq2 = RnaSequence(''.join(seq).replace('-', ''))
                #compare alignment score. Higher is better.
                aln, alnscore = classic_align_pairwise(seq1, seq2, alnscores, -10, -10, False, return_score=True)
                if alnscore > seqscore:
                    strscore = holdscore
                    seqscore = alnscore
                    bestref = keys[secpos]
        if bestref != "":
            if norefseq:
                #realign all sequences since no refseq available
                combinedseqs = fulldict[bestref].degap().addSeqs(fulldict[currstruct].degap())
                fulldict[bestref] = align_unaligned_seqs(combinedseqs, RNA)
                continue
            if not fulldict[currstruct].getGappedSeq("refseq").isGapped():
                fulldict[bestref].addFromReferenceAln(fulldict[currstruct])
            elif not fulldict[bestref].getGappedSeq("refseq").isGapped():
                fulldict[currstruct].addFromReferenceAln(fulldict[bestref])
                fulldict[bestref] = fulldict[currstruct]
            else:
                #realign all sequences since both refseqs have gaps
                #hacky but it works, need to fix later
                fulldict[bestref].Names.remove("refseq")
                combinedseqs = fulldict[bestref].degap().addSeqs(fulldict[currstruct].degap())
                fulldict[bestref] = align_unaligned_seqs(combinedseqs, RNA)
                fulldict[bestref].Names.remove("refseq")
                fulldict[bestref].Names.insert(0, "refseq")
            fulldict.pop(currstruct)
            topop.append(pos)
    topop.sort(reverse=True)
    for pos in topop:
        keys.pop(pos)
    score_structures.end()
    return fulldict, keys


def group_by_distance(structgroupsfile, structscore, specstructs=None, setfinishlen=None, norefseq=False):
        '''Does grouping by way of de-novo reference creation and clustering
            structgroups - dictionary with ALL structures and the Alignment
                           object keyed to them
            structscore - maximum score to consider grouping structures
            specstructs - a list of a subset of structures in structgroups
                          to cluster (optional)
            setfinishlen - Allows manual number of reference structures (default 1%  of dict)
            norefseq - boolean indicating reference sequence is not in each alignment (default False)
        '''
        #read in file info
        fin = open(structgroupsfile, 'rU')
        first = True
        currstruct = ""
        curraln = []
        structgroups = {}
        for head,seq in MinimalFastaParser(fin):
            if first and head == "newaln":
                first = False
                currstruct = seq
            elif head == "newaln":
                structgroups[currstruct] = LoadSeqs(data=curraln, moltype=RNA)
                structgroups[currstruct].Names.remove("refseq")
                structgroups[currstruct].Names.insert(0, "refseq")
                currstruct = seq
                curraln = []
            else:
                curraln.append((head,seq))
        structgroups[currstruct] = LoadSeqs(data=curraln, moltype=RNA)
        structgroups[currstruct].Names.remove("refseq")
        structgroups[currstruct].Names.insert(0, "refseq")
        fin.close()

        #fail if nothing to compare
        if len(structgroups) < 1:
            raise ValueError("Must have at least one structure to group!")
        #return the list directly if only one item (useful for breakout work)
        if len(structgroups) == 1:
            return structgroups
        #just de-novo group if 10 or less to save time and effort
        if len(structgroups) <= 10:
            if specstructs is None:
                structgroups, reference = group_denovo(structgroups, structgroups.keys(), structscore, norefseq)
            else:
                structgroups, reference = group_denovo(structgroups, specstructs, structscore, norefseq)
            return structgroups
        #for speed, get 1% as initial clustering or user defined. need at least 10 structs though
        if setfinishlen == "None":
            finishlen = int(ceil(len(structgroups) * 0.01))
        else:
            finishlen = setfinishlen
        if finishlen < 10 and len(structgroups) >= 10:
            finishlen = 10
        #do initial ref grab by either all structures or specific ones passed
        if specstructs == None:
            reference, nonreference = build_reference(structgroups.keys(), finishlen)
        else:
            reference, nonreference = build_reference(specstructs, finishlen)
        startungrouped = len(structgroups)
        structgroups, reference = group_denovo(structgroups, reference, structscore, norefseq)
        structgroups, ungrouped = group_to_reference(structgroups, reference, nonreference, structscore, norefseq)
        endungrouped = len(structgroups)
        # keep refining while not at limit and are still grouping structs
        while len(ungrouped) > finishlen and startungrouped != endungrouped:
            startungrouped = len(structgroups)
            reference, nonreference = build_reference(ungrouped, finishlen)
            structgroups, reference = group_denovo(structgroups, reference, structscore, norefseq)
            structgroups, ungrouped = group_to_reference(structgroups, reference, nonreference, structscore, norefseq)
            endungrouped = len(structgroups)
        #end while
        structgroups, reference = group_denovo(structgroups, ungrouped, structscore, norefseq)
        return structgroups


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
        stdout.flush()
        #hard limit of 500 sequences to align and fold for memory reasons
        if len(seqs) > 500:
            seqs = seqs[:500]
        #run BayesFold on sequences in the group
        #maxiters set to 5 because should have huge amount of sequences for some groups
        aln, struct = bayesfold(seqs, params={"-diags": True})
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

