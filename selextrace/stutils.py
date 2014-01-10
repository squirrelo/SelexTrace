
from cogent import LoadSeqs
from cogent.app.uclust import Uclust, clusters_from_uc_file
from subprocess import Popen, PIPE
from qiime.split_libraries import local_align_primer_seq


def write_fasta_list(lst, filename):
    '''writes formatted list [(header,sequence)] to fasta file filename'''
    fileout = open(filename, 'w')
    for seq in lst:
        fileout.write('>%s\n%s\n' % (seq[0], seq[1]))
    fileout.close()


def write_fasta_dict(dct, filename):
    '''writes formatted dict {header: sequence} to fasta file filename'''
    fileout = open(filename, 'w')
    for header, seq in dct:
        fileout.write('>%s\n%s\n' % (header, seq))
    fileout.close()


def strip_primer(seqs, primer, maxmismatch=0, keep_primer=False):
    '''strips 3 prime primer from sequences in fasta file and returns MinimalFastaParser
    formatted arrays for stripped and not stripped sequences'''
    nostripped = []
    stripped = []
    pri = primer.upper()
    for head, seq in seqs:
        RNA = False
        seq = seq.upper()
        if 'U' in seq:
            seq = seq.replace('U', 'T')
            RNA = True
        #code adapted from truncate_reverse_primers.py in qiime
        rev_primer_mm, rev_primer_index =\
            local_align_primer_seq(pri, seq)
        if rev_primer_mm > maxmismatch:
            nostripped.append((head, seq))
            continue
        if keep_primer:
            seqnew = seq[:rev_primer_index + len(primer)]
        else:
            seqnew = seq[:rev_primer_index]
        if RNA:
            seqnew = seqnew.replace('T', 'U')
        stripped.append((head, seqnew))
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


def remove_duplicates(seqsin):
    '''Takes in LoadSeqs loadable sequences, removes duplicate sequences
    and returns a list of unique sequence tuples, formated (sequence, count)
    sorted most abundant to least abundant'''

    parsable_seqs = LoadSeqs(data=seqsin, aligned=False)
    uniques = {}
    for header, seq in parsable_seqs.items():
        seq = str(seq)
        if seq in uniques:
            uniques[seq] += 1
        else:
            uniques[seq] = 1
    uniques = uniques.items()
    uniques.sort(key=lambda x: x[1], reverse=True)
    return uniques

def get_shape(struct):
    '''Converts a dot-bracket notation to abstract shape notation'''
    p = Popen(["RNAshapes", "-D", struct], stdout=PIPE)
    return p.communicate()[0].strip()