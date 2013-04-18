from cogent.app.uclust import Uclust, clusters_from_uc_file
from cogent.parse.fasta import MinimalFastaParser
from sys import argv
from os.path import exists
from os import mkdir

#clustertest.py /path/to/unique/fasta/in.fasta path/to/folder/out simmilarity

if __name__ == "__main__":
    if argv[2][-1] != "/":
        argv[2] += "/"
    if not exists(argv[2]):
        mkdir(argv[2])

    params = {
        '--usersort': False,
        '--id': float(argv[3]),
        '--maxaccepts': 20,
        '--maxrejects': 500,
        '--stepwords': 20,
        '--w': 12,
        '--gapopen': '50.0',
        '--gapext': '50.0'
    }
    uclust = Uclust(params, WorkingDir='/tmp')
    input_data = {
        '--input': argv[1],
        '--uc': argv[2] + argv[3] + "_clusters.uc",
        '--log': argv[2] + argv[3] + "_clusters.log"
    }
    result = uclust(input_data)
    clusters, failures, new_seeds = clusters_from_uc_file(result['ClusterFile'])
    print "RESULTS: ", len(clusters)
    headers = {}
    for header, seq in MinimalFastaParser(open(argv[1])):
        headers[header.split()[0]] = header
    otusout = open(argv[2] + argv[3] + "_clusters.txt", 'w')
    for group, cluster in enumerate(clusters):
        otusout.write(str(group) + "\t")
        #map headers back to orignal ones with counts
        for header in clusters[cluster]:
            otusout.write(headers[header] + "\t")
        otusout.write("\n")
    otusout.close()
    #log = open(argv[2] + argv[3] + "_clusters.log", 'w')
    #log.write('\n'.join(input_data) + '\n'.join(params) +
    #    "clusters: " + str(len(clusters)) + "\n")
