from cogent import LoadSeqs
from cogent.app.uclust import Uclust, clusters_from_uc_file
from cogent.parse.fasta import MinimalFastaParser
from sys import argv
from os.path import exists
from os import mkdir

#clustertest.py /path/to/unique/fasta/in.fasta path/to/folder/out simmilarity
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
    if gapopen != None:
        params['--gapopen'] = gapopen
    if gapext != None:
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

#if __name__ == "__main__":
#    if argv[2][-1] != "/":
#        argv[2] += "/"
#    if not exists(argv[2]):
#        mkdir(argv[2])#

#    print "RESULTS: ", len(clusters)
#    headers = {}
#    for header, seq in MinimalFastaParser(open(argv[1])):
#        headers[header.split()[0]] = header
#    otusout = open(argv[2] + argv[3] + "_clusters.txt", 'w')
#    for group, cluster in enumerate(clusters):
#        otusout.write(str(group) + "\t")
#        #map headers back to orignal ones with counts
#        for header in clusters[cluster]:
#            otusout.write(headers[header] + "\t")
#        otusout.write("\n")
#    otusout.close()
#    #log = open(argv[2] + argv[3] + "_clusters.log", 'w')
#    #log.write('\n'.join(input_data) + '\n'.join(params) +
#    #    "clusters: " + str(len(clusters)) + "\n")
