def pickotu_to_fasta(fasta, otulist, folderout, minseqs=0):
    '''converts list of sequence names from pick_otu to seperate fasta files
    for each OTU found.

    fasta: fasta file used to pick otulist
    otulist: output file of pick_otus.py
    folderout: output folder for the OTU fasta files
    minseqs: (optional) minimum sequences in an OTU required before making a fasta file'''

    from cogent.parse.fasta import MinimalFastaParser

    fastain = open(fasta, 'r')
    seqs = {}
    for header, seq in MinimalFastaParser(fastain):  # Clean seq names so only using the identifier before first space
        seqs[header] = seq
    fastain.close()

    otusin = open(otulist, 'rU')
    summary = []
    for line in otusin:
        splitnames = line.strip().split("\t")
        #get sequence count
        numseqs = 0
        for header in splitnames[1:]:
            numseqs += int(header.split("_")[1])
        if numseqs > minseqs:
            summary.append((splitnames[0], numseqs))
            otuout = open(''.join([folderout, "cluster_", splitnames[0], ".fasta"]), 'w')
            for i in range(1, len(splitnames)):
                otuout.write(''.join([">", splitnames[i], "\n", seqs[splitnames[i]], "\n"]))
            otuout.close()
    otusin.close()

    summary.sort(reverse=True, key=lambda seq: seq[1])
    pathfile = open(folderout + "cluster_paths.txt", 'w') #file with OTU name and fasta path
    summaryfile = open(folderout + "cluster_summary.txt", 'w')  #file containting OTU names and sequence counts for each OTU
    summaryfile.write(str(len(summary)) + " total clusters Min " + str(minseqs) + " sequences\n")
    for otu in summary:
        summaryfile.write(''.join(["cluster_", otu[0], "\t", str(otu[1]), "\n"]))
        pathfile.write(''.join(["cluster_", otu[0], " ", folderout, "cluster_", otu[0], ".fasta\n"]))
    summaryfile.close()


def help():
    print "Breaks OTUs into seperate fasta files of sequences for each OTU"
    print ""
    print "USE: pickotu_to_fasta.py /path/to/fastafile.fa /path/to/otufile.txt /path/to/output/folder [minseqs]"


if __name__ == "__main__":
    from os.path import isdir
    from os import makedirs
    from sys import argv, exit
    if len(argv) < 4:
        help()
        exit(0)

    if argv[3][:-1] != "/":
        argv[3] += "/"

    if not isdir(argv[3]):
        makedirs(argv[3])

    minseqs = 0
    if len(argv) > 4:
        minseqs = int(argv[4])
    print "minseqs:" + str(minseqs)
    pickotu_to_fasta(argv[1], argv[2], argv[3], minseqs)
