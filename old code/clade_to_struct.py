#clade_to_struct.py /path/to/clade/list.txt /path/to/dupes.fasta /path/to/tree.newick /path/to/folder/out 

if __name__ == "__main__":
    from sys import argv, exit
    from os import mkdir
    from os.path import isdir
    from cogent.parse.fasta import MinimalFastaParser
    from cogent import LoadSeqs, RNA, LoadTree
    from cogent.app.muscle_v38 import align_unaligned_seqs
    from subprocess import check_call
    
    #directory where the PPfold jar is kept
    PPFOLDDIR = "/Users/Ely/bin/"

	#cleaing up the output directory string and creating the folder if necessary 
    folderout = argv[4]
    if folderout[:-1] == "/":
        folderout = folderout[0:-1]
    if not isdir(folderout):
        mkdir(folderout)

	#getting the base name to use for all files as the folder name
    basenames = folderout
    if "/" in basenames:
        basenames = basenames[basenames.rfind("/") + 1:]

    print "Creating fasta of sequences in clade"
    fastain = open(argv[2], 'rU')
    seqs = {}
    #load in all the sequences in the fasta used to make the tree
    for name, seq in MinimalFastaParser(fastain):
        seqs[name.split()[0]] = seq
    fastain.close()
	
	#load in all haeders for sequences in the clade and match them to their sequence
    listin = open(argv[1], 'rU')
    fileout = open(folderout + "/" + basenames + "-seqs.fasta", 'w')
    listin.readline()
    rawseqs = []
    tips = []
    for line in listin:
        header = line.split()[0]
        fileout.write(''.join([">", header, "\n", seqs[header], "\n"]))
        rawseqs.append((header, seqs[header]))
        tips.append(header)
    fileout.close()

    print "Aligning seqs using muscle with -diags"
    seqs = LoadSeqs(data=rawseqs, moltype=RNA, aligned=False)
    aln = align_unaligned_seqs(seqs, RNA, {"-diags": True})
    fileout = open(folderout + "/" + basenames + "-seqsaligned.fasta", 'w')
    fileout.write(str(aln))
    fileout.close()

    print "Folding sequences"
    #get subtree of the clade being folded to pass to PPfold
    tr = LoadTree(argv[3])
    sub_tree = tr.getSubTree(tips, keep_root=True)
    filesubtree = open(folderout + "/" + basenames + "-subtreeDistances.nwk", 'w')
    filesubtree.write(sub_tree.getNewick(with_distances=True))
    filesubtree.close()
    filesubtree = open(folderout + "/" + basenames + "-subtree.nwk", 'w')
    filesubtree.write(sub_tree.getNewick(with_distances=False))
    #call PPfold with aligned sequences and subtree
    args = ["java", "-jar", PPFOLDDIR + "PPfold.jar", folderout + "/" + basenames + "-seqsaligned.fasta", "--outputd", folderout]
    check_call(args)
    
    print "Converting sequences to vienna"
    check_call(["ct2b.pl", folderout + basenames + "-seqsaligned.ct", ">",folderout + basenames + "-vienna.txt"]) 

    print "DONE"
