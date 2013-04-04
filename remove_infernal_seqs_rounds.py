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
    for r in range (1, int(argv[3]) + 1):
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
                    hit = line.split(",")[0].split("_")[0]
                    if hit not in hits:
                        hits.add(hit)
        

        #make the new unique sequences file
        print "Writing uniques file"
        if exists(rounddir + curround + "/" + curround + "-Unique-Remaining.fasta")
            uniques = open(rounddir + curround + "/" + curround + "-Unique-Remaining.fasta", 'rU')
        else:
            uniques = open(rounddir + curround + "/" + curround + "-Unique.fasta", 'rU')
        uniques_seqs = uniques.readlines()
        uniques.close()
        newuniques = open(rounddir + curround + "/" + curround + "-Unique-Remaining.fasta", 'w')
        for header, seq in MinimalFastaParser(uniques_seqs):
            if header.split("_")[0] not in hits:
                newuniques.write(">" + header + "\n" + seq + "\n")
        uniques = 0
        newuniques.close()
        #make new clean stripped file. More dificult because one hit to multiple seqs
        print "Writing CleanStripped file"
        #load current round's map of sequences to headers
        mapfile = open(rounddir + curround + "/" + curround + "-seqtoheaders.txt")
        seqtoheaders = {}
        for line in mapfile:
            lineinfo = line.split("\t")
            seqtoheaders[lineinfo[0]] = lineinfo[1].split(',')[:-1]
        mapfile.close()

        #load cleanstripped into a dictionary
        cleanstrip = {}
        for header, seq in MinimalFastaParser(open(rounddir + curround + "/" + curround + "-CleanStripped.fasta", 'rU')):
            cleanstrip[header] = seq
        #remove every copy of a sequence that maps to the unique one in mapfile
        for hit in hits:
            #remove the count from the header and get the sequence
            hitseq = cleanstrip[hit]
            for remove in seqtoheaders[hitseq]:
                cleanstrip.pop(remove)
        #now cleanstrip only contains sequences that havent been found already
        #remap headers and print out remaining sequs to file
        newcleanstripped = open(rounddir + curround + "/" + curround + "-CleanStripped-Remaining.fasta", 'w')
        for key in cleanstrip:
             newcleanstripped.write(">" + key + "\n" + cleanstrip[key] + "\n")
        newcleanstripped.close()

            