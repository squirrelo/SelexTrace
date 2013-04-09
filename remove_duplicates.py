#!/usr/bin/env python
"""Removes duplicate sequences from a fasta file. Prints unique sequences as
   csv and fasta files with counts
"""

__author__ = "Joshua Shorenstein"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Joshua Shorenstein"]
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Joshua Shorenstein"
__email__ = "joshua.shorenstein@colorado.edu"
__status__ = "Prototype"

from cogent.parse.fasta import MinimalFastaParser
from collections import defaultdict

def remove_duplicate(lst, seq):
    '''takes in MinimalFastaParser formatted sequences list and a sequence seq,
    returns a list with all occurences of seq removed'''
    return [x for x in lst if x[1] != seq]

def remove_duplicates(seqsin):
    '''Takes in minimalfastaparser formatted list, removes duplicate sequences 
    and returns a MinimalFastaParser formatted list of unique sequences (with 
    one sequence from each duplicate set)'''
    
    uniques = {}
    counts = defaultdict(int)
    #repseqs: dict of lists keyed to a sequence. Holds all headers of duplicate 
    #sequences with that sequence
    repseqs = defaultdict(list)
    for seq in seqsin:
        if seq[1] not in uniques:
            #append unique sequence to uniques list
            uniques[seq[1]] = seq[0]
        #first item in repseqs list will be header of representative sequence
        repseqs[seq[1]].append(seq[0])
        counts[seq[1]] += 1
    #append count to end of header
    uniquesret = []
    for unique in uniques:
        #get count for this seuqnce
        count = counts[unique]
        #append to header and add to return list
        uniquesret.append((uniques[unique] + "_" + str(count), unique))
    #sort by sequence count, highest to lowest
    uniquesret.sort(reverse=True, key=lambda item: int(item[0].split("_")[1]))
    return uniquesret, repseqs

if __name__ == "__main__":
    from sys import argv, exit
    folderout = argv[2]
    if folderout[:-1] != "/":
        folderout += "/"
    if len(argv) < 3:
        print "remove_duplicates.py /path/to/finein.fasta /path/to/folder/out"
        exit(1)
    seqs = [(header, seq) for header, seq in MinimalFastaParser(open(argv[1]))]
    uniques = remove_duplicates(seqs)
    csvout = open(folderout + "Unique_sequences.csv", 'w')
    csvout.write("Count,Header,Sequence")
    fastaout = open(folderout + "Unique_sequences.fasta", 'w')
    for seq in uniques:
        fastaout.write('>%s\n%s\n' % (seq[0], seq[1]))
        csvout.write('%s,%s,%s\n' % (seq[2],seq[0],seq[1]))
