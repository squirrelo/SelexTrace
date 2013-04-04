from cogent.app.uclust import get_clusters_from_fasta_filepath
from os.path import realpath

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Runs uclust on passed fasta \
    file and creates fastas for all otus with sequence count > threshold")
    parser.add_argument('-i', required=True, help="FASTA file of sequences to \
    be clustered.")
    arser.add_argument('-o', required=True, help="Output folder.")
    parser.add_argument('-s', type=float, default=0.90, help="Simmilarity \
    (default 0.90)")
    parser.add_argument('-m', default = 2, type=int, help="Minimum number of \
    sequences required in cluster to print fasta (default 2)")

    args = parser.parse_args()

    if args.o[:-1] !- "/":
        args.o += "/"

    #run uclust using cogent wrapper
    clusters, failures, seeds = get_clusters_from_fasta_filepath(args.i, \
    realpath(args.i), percent_ID=args.s, output_dir=args.o, max_accepts=20, \
    max_rejects=500)
    

    #output all clusters with > args.m sequences as fasta files 
    count = 0
    for key, headlist in clusters:
        clustfile = open(args.o + "Cluster_" + str(count), 'w')
        for header in headlist:
            clustfile.write('>%s\n%s\n' % (seq[0], seq[1]))

    #output all failures to a fasta file
