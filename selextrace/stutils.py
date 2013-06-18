
from cogent.parse.fasta import MinimalFastaParser
from collections import defaultdict
from cogent import LoadSeqs
from cogent.app.uclust import Uclust, clusters_from_uc_file
from cogent.parse.fasta import MinimalFastaParser
from sys import argv
from os.path import exists
from subprocess import Popen, PIPE
from qiime.split_libraries import local_align_primer_seq


def write_fasta_list(lst, filename):
    '''writes MinimalFastaParser formatted list [(header,sequence)] to fasta file filename'''
    fileout = open(filename, 'w')
    for seq in lst:
        fileout.write('>%s\n%s\n' % (seq[0], seq[1]))
    fileout.close()


def write_fasta_dict(dct, filename):
    '''writes MinimalFastaParser formatted dict {header: sequence} to fasta file filename'''
    fileout = open(filename, 'w')
    for header, seq in dct:
        fileout.write('>%s\n%s\n' % (header, seq))
    fileout.close()


def strip_primer(seqfile, primer, maxmismatch=0, keep_primer=False):
    '''strips 3 prime primer from sequences in fasta file and returns MinimalFastaParser
    formatted arrays for stripped and not stripped sequences'''
    seqs = MinimalFastaParser(open(seqfile, 'rU'))
    nostripped = []
    stripped = []
    for seq in seqs:
        #code adapted from truncate_reverse_primers.py in qiime
        rev_primer_mm, rev_primer_index =\
            local_align_primer_seq(primer, seq[1])
        if rev_primer_mm > maxmismatch:
            nostripped.append(seq)
        else:
            seqnew = seq[1][:rev_primer_index]
            if keep_primer:
                seqnew += primer
            stripped.append((seq[0], seqnew))
    #end for
    return stripped, nostripped


def cluster_seqs(seqspath, simmilarity, folderout='/tmp', gapopen=None, gapext=None):
    if folderout[-1] != "/":
        folderout += "/"

    params = {
        '--usersort': False,
        '--id': float(simmilarity),
        '--maxaccepts': 20,
        '--maxrejects': 500,
        '--stepwords': 20,
        '--hsp': 0
    }
    if gapopen is not None:
        params['--gapopen'] = gapopen
    if gapext is not None:
        params['--gapext'] = gapext
    uclust = Uclust(params, WorkingDir='/tmp')
    input_data = {
        '--input': seqspath,
        '--uc': folderout + "clusters.uc",
        '--log': folderout + "clusters.log"
    }
    result = uclust(input_data)
    clusters, failures, new_seeds = clusters_from_uc_file(result['ClusterFile'])

    seqs = LoadSeqs(seqspath, aligned=False)
    convheader = {}
    clusterseqs = {}
    #create dictinary to convert shortened headers to full headers
    for header in seqs.getSeqNames():
        convheader[header.split()[0]] = header
    #match headers in each cluster to seqs to create cluster tuples list
    for num, cluster in enumerate(clusters):
        clusterseqs["cluster_" + str(num)] = []
        for header in clusters[cluster]:
            clusterseqs["cluster_" + str(num)].append((convheader[header], seqs.getSeq(convheader[header])))

    return clusterseqs


def remove_duplicate(lst, seq):
    '''Generic function for removing duplicates of a sequence in a list
    Takes MinimalFastaParser formatted tuple list'''
    return [x for x in lst if x[1] != seq]


def remove_duplicates(seqsin):
    '''Takes in minimalfastaparser formatted list, removes duplicate sequences
    and returns a MinimalFastaParser formatted list of unique sequences (with
    one sequence from each duplicate set)'''

    uniques = {}
    counts = defaultdict(int)
    #uniquesret: standard (header, seq) tuple list of unique sequeces
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


def get_shape(struct):
    '''Converts a dot-bracket notation to abstract shape notation'''
    p = Popen(["rnashapes", "-D", struct], stdout=PIPE)
    return p.communicate()[0].strip()