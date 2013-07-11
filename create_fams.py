from selextrace.ctilib import score_local_rnaforester
from cogent.parse.fasta import MinimalFastaParser
from sys import argv
from os import walk
from time import time
from multiprocessing import Pool

#create_seqs.py /path/to/base/folder/ #forestercutoff #cpus 

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
        #remove empty and add remaining structs to group
        if "" in scores:
            scores.discard("")
        for struct in scores:
            if struct in fulldict:
                fulldict[currstruct].extend(fulldict[struct])
                fulldict.pop(struct)
    return fulldict

if __name__ == "__main__":
    basefolder = argv[1]
    if basefolder[-1] != "/":
        basefolder += "/"
    groups = {}
    for group in walk(basefolder).next()[1]:
        if group == "fasta_groups":
            continue
        filein = open(argv[1] + "/" + group + "/bayesfold-aln.fasta")
        for header, seq in MinimalFastaParser(filein):
            if "(" in seq:
                groups[seq] = [group]
        filein.close()
    secs = time()

    print len(groups), "starting groups"

    groups = group_by_forester(groups, 200, int(argv[3]))
    print len(groups), "ending groups"
    print (time() - secs)/60, "mins"
    fout = open(argv[1] + "families" + str(argv[2]) + ".txt", 'w')
    for fam in groups:
        fout.write("\t".join(groups[fam]) + "\n")
    fout.close()
