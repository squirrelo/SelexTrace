#dupestest.py /path/to/OTU.fasta /path/to/folder/out CPUs
from sys import argv, exit
from os.path import exists
from os import mkdir
from cogent import LoadSeqs, RNA
from cogent.app.locarna import create_locarnap_alignment
from cogent.app.infernal_v11 import cmsearch_from_file, cmbuild_from_alignment
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from remove_duplicates import remove_duplicates
from numpy import log2

#locarna_test.py path/to/otuslist.txt

def find(lst, searchstring):
    '''returns first occurence of searchstring in list or -1 if not found
    searchstring can be a substring of a string in list'''
    try:
        hold = min([i for i, x in enumerate(lst) if searchstring in x])
        return hold
    except StopIteration:
        return -1

if __name__ == "__main__":
    #Get list of OTUs from file, populate dictionary
    #need to add -i for input, -o for out folder, -c for cpus, -r current round
    otus = []
    fn = open(argv[1], 'rU')
    for line in fn:
        lineinfo = line.strip().split()
        otus.append((lineinfo[0], lineinfo[1]))
    fn.close()
    for currotu in otus:
        otu = currotu[0]
        print "==" + otu + "=="
        print "Reading in 30 most abundant sequences"
        #assuming that the fasta has more than 30 sequences in it. Safe assumption
        #if this is a significant cluster
        seqs = [(header, seq) for header,seq in MinimalFastaParser(open(currotu[1], 'rU'))]
        seqs, headers = remove_duplicates(seqs)
        #blank headers to save memory
        headers = 0
        #headers come out in format Header_# so split to get # and sort by abundance
        seqs.sort(reverse=True, key=lambda count:int(count[0].split('_')[1]))
        #cut to 30 most abundant sequences
        seqs = seqs[:30]
        print "Running locarna-p on sequences"
        args = {'--cpus':'24'}
        aln, struct = create_locarnap_alignment(seqs, RNA, struct=True, params=args)
        #create output folder for OTU
        otufolder = "/Users/Ely/Desktop/Ely_selection/R7/lead_clusters/"
        if not exists(otufolder):
            mkdir(otufolder)
        otufolder += otu
        if not exists(otufolder):
            mkdir(otufolder)
        #print out alignment and structure in fasta and stockholm formats
        alnout = open(otufolder + "/locarnap-aln.fasta", 'w')
        alnout.write(aln.toFasta() + "\n>SS_struct\n" + struct + "\n")
        alnout.close()
        alnout = open(otufolder + "/locarnap-aln.sto", 'w')
        struct_dict = {'SS_cons':struct}
        alnout.write(stockholm_from_alignment(aln,GC_annotation=struct_dict))
        alnout.close()
        print struct
            #remove found sequences from the round files
     
