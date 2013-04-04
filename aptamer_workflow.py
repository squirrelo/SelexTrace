
from strip_primers import strip_primer
from sys import is_dir, mkdir, exit
from cogent.parse.fasta import MinimalFastaParser
#from qiime.pick_otus  import otu_picking_method_constructors
from subprocess import check_call, Popen, PIPE
from remove_duplicates import remove_duplicate, remove_duplicates
from tempfile import mktemp

#path to the Qiime scripts folder
PICKOTUS = "/Users/Ely/bin/Qiime/scripts/"
#threshold # of sequnces to be relevant

def write_fasta_list(lst, filename):
    '''writes MinimalFastaParser formatted list to fasta file filename'''
    fileout = open(filename, 'w')
    for seq in lst:
        fileout.write('>%s\n%s\n' % (seq[0], seq[1]))
    fileout.close()

def write_fasta_dict(dct, filename):
    '''writes MinimalFastaParser formatted list to fasta file filename'''
    fileout = open(filename, 'w')
    for header, seq in dct:
        fileout.write('>%s\n%s\n' % (header, seq))
    fileout.close()


def rem_N(seqs):
    '''Takes in a MinimalFastaParser formatted list of sequences and returns list 
    with N containing sequences removed'''
    rem = []
    for i, seq in enumerate(seqs):
        if "N" in seq[1].upper():
            rem.append(i)
        rem.sort(reverse=True)
    for i in rem:
        seqs.pop(i)
    return seqs


def rem_short(seqs, minlen):
    '''Takes in a MinimalFastaParser formatted list of sequences and returns list 
    with sequences shorter than minlen removed'''
    rem = []
    for i, seq in enumerate(seqs):
        if len(seq[1]) < minlen:
            rem.append(i)
        rem.sort(reverse=True)
    for i in rem:
        seqs.pop(i)
    return seqs


def run_uclust_pick_otus(input_seqs, output_dir, simmilarity = 0.90):
    '''Runs uclust pick_otus algorithm on input_seqs FASTA file using default 
    settings with simmmilarity simmilarity and write to output_dir'''
    if output_dir[:-1] != "/":
        output_dir += "/"

    otu_picker_constructor =\
     otu_picking_method_constructors['uclust']
    params = {
        'Similarity':simmilarity
        'enable_rev_strand_matching':False,
        'optimal':False,
        'exact':False,
        # suppress_sort=True when seqs are or will be pre-sorted
        'suppress_sort':False,
        'presort_by_abundance': True,
        'max_accepts':20,
        'max_rejects':500,
        'stepwords':20,
        'word_length':12,
        'new_cluster_identifier':None,
        'stable_sort':True,
        'save_uc_files':True,
        'output_dir':output_dir,
        'prefilter_identical_sequences':True
    }
    otu_picker = otu_picker_constructor(params)
    otu_picker(input_seqs, result_path=output_dir + "clusters.txt" \
    log_path=output_dir + "clusters.log", HALT_EXEC=False)


def find(lst, searchstring):
    '''returns first occurence of searchstring in list or -1 if not found
    searchstring can be a substring of a string in list'''
    try:
        hold = min([i for i, x in enumerate(lst) if searchstring in x])
        return hold
    except StopIteration:
        return -1


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Runs aptamer workflow on \
    FASTA files in directory. Each file is a round of selection.")
    parser.add_argument('-i', required=True, help="Folder containing fasta \
    files for each round of selection. Files must be in name format R1.fasta, \
    R2.fasta etc")
    parser.add_argument('-p', default="", required=True, help="3' primer \
    sequence to strip")
    parser.add_argument('-r', type=int, required=True, help="Number of rounds \
    in selection")
    parser.add_argument('-m', default = -1, type=int, help="minimum length of \
    sequence required to not be removed from pool (default 0.75 * longest seq)")
    parser.add_argument('-c', default = 1, type=int, "Number of CPUs to use \
    default 1)")

    args = parser.parse_args()

    if args.m != -1:
        if args.m < 0:
            print "ERROR: mim sequence length must be greater than 1!"
            exit(1)

    if args.c < 1:
            print "ERROR: Can't use less than 1 CPU!"
            exit(1)

    #add trailing slash if necessary to folderin
    folderin = args.i
    if folderin[:-1] != '/':
        folderin += "/"

    rounds = args.r

    print "===================="
    print "Folderin: " + args.i
    print "3' primer: " + args.p
    print "Number of rounds: " + str(args.r)
    print "Number of CPUs: " + stre(args.c)
    print "====================\n"

    print "==Cleaning input sequences=="
    #make directory to store cleaned sequences in
    if not is_dir(folderin + "cleaned_seqs"):
        mkdir(folderin + "cleaned_seqs")

    currfolder = folderin + "cleaned_seqs/" 
    #keep sequences for later use in roundseqs dict, is dict of dicts where 
    #each round is a key in the dict and holds a dict of unique sequences keyed
    #to headers for round
    roundseqs = {}
    for i in range(1, rounds + 1):  #cycle through all rounds fasta files
        print folderin + "R" + str(i) + currfile + ".fasta"
        #strip primers from sequences, print out not stripped to be safe
        kept, rem = strip_primer(currfile + ".fasta", args.p, 2)
        print str(len(kept)) + " sequences left, " + str(len(rem)) + \
        " sequences removed"
        write_fasta_list(rem, currfolder + "R" + str(i) + "-NotStripped.fasta")
        #remove all sequences with Ns and short sequences
        kept = rem_N(kept)
        kept = rem_short(kept, 75)
        print str(len(kept)) + " sequences left"
        write_fasta_list(kept, currfolder + "R" + str(i) + "-CleanStripped.fasta")

        #remove duplicate sequences from the fasta file and store for later
        kept = remove_duplicates(kept)
        roundseqs[i] = {}
        for seq in kept:
            roundseqs[i][seq[0]] = seq[1]
        write_fasta_list(kept, currfolder + "R" + str(i) + "-CleanStrippedUnique.fasta")

    print "==ENTERING ITERATIVE PHASE=="
    if not is_dir(folderin + "clusters"):
        mkdir(folderin + "clusters")
    if not is_dir(folderin + "locarna"):
        mkdir(folderin + "locarna")
    if not is_dir(folderin + "infernal"):
        mkdir(folderin + "infernal")

    for i in range(rounds + 1, 1, -1):  #start at highest round and go down to 1
        print "====ROUND "+ str(i) +"===="
        print "==Picking clusters=="
        #running uclust pickOTUs from qiime on sequences in roundseqs for round
        #NEED TO MAKE THIS QIIME INDEPENDENT OR BETTER THAN DIRECT CALL
        fn = mktemp(prefix="R" + str(i) + "clusterseqs_", suffix=".fasta")
        write_fasta_dict(roundseqs[i], fn)
        check_call(['.' + PICKOTUS + 'pick_otus.py', '-i', fn, \
        '-o', folderin + "clusters", '-s', "0.90"])
        print "==Sorting clusters by abundance of sequences=="
        #make folder to store round's cluster fastas
        if not is_dir(folderin + "clusters/R"+str(i)):
            mkdir(folderin + "clusters/R"+str(i))
        fn = fn[fn.rfind("/") + 1:fn.rfind(".")]
        otufile= open(folderin + "clusters/" + fn + "_otus.txt",'rU')
        #parse through OTU file and output fasta files for each significant OTU
        otucount = []
        for line in otufile:
            lineinfo = line.split()
            count = len(lineinfo) - 1
            for header in lineinfo[1:]:
                count += int(header.split("_")[1])
            if count > THRESHOLD:
                #sort from most abundant to least abundant by header info
                #and output fasta
                temparray = []
                for header in lineinfo[1:]:
                    temparray.append((header, int(header.split("_")[1])))
                temparray.sort(reverse=True, key = lambda pair: pair[1])
                fileout.open(folderin + "clusters/R"+str(i) + "/" + lineinfo[0] \
                + ".fasta", 'w')
                for item in temparray:
                    fileout.write('>%s\n%s\n' % (item[0], roundseqs[i][item[0]]))
                fileout.close()
                otucount.append(folderin + "clusters/R"+str(i) + "/" + lineinfo[0] \
                + ".fasta", count)
        #end for loop through OTU file
        #sort OTU list by sequence abundance, most to least
        otucount.sort(reverse=True, key = lambda pair: pair[1])
        
        #make sequence count list from OTU table
        print "==Folding most abundant sequences=="
        #Run locana-P on 50 most abundant sequences for each OTU
        print "==Running Infernal over all rounds=="
        #Use the 50 seqs + locarna-P structure to create CM model for infernal
        #run infernal on each pool unique sequences for these models
    print "==Compiling results=="
    #need to repeat above starting with last pool, moving to first, removing all used sequences as we go.
