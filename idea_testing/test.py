from cluster_seqs import cluster_seqs
from sys import argv
from cogent import RNA, LoadSeqs
from cogent.app.muscle_v38 import align_unaligned_seqs
from weblogolib import *
from os.path import exists
# test.py /path/to/filein.fasta /path/to/folderout/ #simmilarity

if __name__ == "__main__":
    if argv[2][-1] != "/":
        argv[2] += "/"

    clusters = cluster_seqs(argv[1], argv[3], folderout=argv[2], gapopen='10.0/*TI', gapext='10.0')

    #remove tiny clusters
    topop = []
    countarray = []
    for cluster in clusters:
        totalseqs = 0
        for seq in clusters[cluster]:
            totalseqs += int(seq[0].split("_")[1])
        if totalseqs < 100:
            topop.append(cluster)
        else:
            countarray.append((cluster,totalseqs))
    for remove in topop:
        clusters.pop(remove)
    countarray.sort(reverse=True, key=lambda c: c[1])

    for c in countarray:
        print c[0] + "\t" + str(c[1])

    print str(len(clusters)) + " clusters"

    #print that shit to file
    clustfiles = []
    for cnum,cinfo in enumerate(countarray):
        if not exists(''.join([argv[2],"cluster", str(cnum), "_", str(cinfo[1]), ".txt"])):
            clusters[cinfo[0]].sort(reverse=True,key=lambda count: int(count[0].split('_')[1]))
            cout = open(''.join([argv[2],"cluster", str(cnum), "_", str(cinfo[1]), ".txt"]), 'w')
            for seq in clusters[cinfo[0]]:
                cout.write(">%s\n%s\n" % seq)
            cout.close()
        clustfiles.append(''.join([argv[2],"cluster", str(cnum), "_", str(cinfo[1]), ".txt"]))
    for cfile in clustfiles:
        print cfile
        seqs = LoadSeqs(cfile, moltype=RNA, format='fasta')
        aln = align_unaligned_seqs(seqs,RNA)
        cout = open(cfile[:cfile.rfind(".")] + "_aln.txt", 'w')
        cout.write(aln.toFasta())
        cout.close()
        fin = open(cfile[:cfile.rfind(".")] + "_aln.txt")
        seqs = read_seq_data(fin) 
        fin.close()
        data = LogoData.from_seqs(seqs)
        options = LogoOptions()
        options.title = "cluster logo"
        format = LogoFormat(data, options)
        fout = open(cfile[:cfile.rfind(".")] + ".eps", 'w')
        eps_formatter(data, format, fout)
        fout.close()