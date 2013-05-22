from os import walk
from os.path import exists
from sys import argv, exit
from cogent.parse.fasta import MinimalFastaParser

if __name__ == "__main__":
    if len(argv) < 4:
        print "remove_infernal_seqs_rounds.py /path/to/groups/folder /path/to/base/selection/folder #rounds"
        print "ex: remove_infernal_seqs_rounds.py Ely_Selection/R6/97percent-lead-groups Ely_Selection 6"
        exit(1)

    startdir = argv[1]
    if startdir[:-1] != "/":
        startdir += "/"

    rounddir = argv[2]
    if rounddir[:-1] != "/":
        rounddir += "/"

    #go through each round seperately
    for r in range(1, int(argv[3]) + 1):
        print "Round " + str(r)
        curround = "R" + str(r)
        #make a list of hits
        hits = set([])
        #iterate through all hroup folders and load headers of hits
        for dirname in walk(startdir).next()[1]:
            if exists(startdir + dirname + "/" + curround + "hits.txt"):
                hitsfile = open(startdir + dirname + "/" + curround + "hits.txt", 'rU')
                #get rid of header
                hitsfile.readline()
                for line in hitsfile:
                    hit = line.split(",")[0]
                    if hit not in hits:
                        hits.add(hit)

        #make the new unique sequences file
        print "Writing uniques file"
        if exists(rounddir + curround + "/" + curround + "-Unique-Remaining.fasta"):
            uniques = open(rounddir + curround + "/" + curround + "-Unique-Remaining.fasta", 'rU')
        else:
            uniques = open(rounddir + curround + "/" + curround + "-Unique.fasta", 'rU')
        uniques_seqs = uniques.readlines()
        uniques.close()
        newuniques = open(rounddir + curround + "/" + curround + "-Unique-Remaining.fasta", 'w')
        for header, seq in MinimalFastaParser(uniques_seqs):
            if header not in hits:
                newuniques.write(">" + header + "\n" + seq + "\n")
        newuniques.close()
