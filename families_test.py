from sys import argv
from os import walk
from numpy import zeros

#overlap_check.py basefolder #round

if __name__ == "__main__":
    basefolder = argv[1]
    if basefolder[:-1] != "/":
        basefolder += "/"

    rnd = "R" + argv[2]


    #create csv file for data
    #empty first cell of header

    #get all group names
    groups = [g for g in walk(basefolder).next()[1]]
    if "fasta_groups" in groups:
        groups.remove("fasta_groups")
    sizegroups = len(groups)
    groups.sort()
    #get total number fo groups in the round, add 1 for looping
    print str(sizegroups) + " groups"
    #compare all group sequences to other sequences
    countmatrix = zeros(shape=(sizegroups, sizegroups), dtype=float)
    for group, gfile in enumerate(groups):
        #load in headers of current group
        headers = set([])
        firstgroup = open(basefolder + gfile + "/" + rnd + "hits.txt")
        #get rid of header
        firstgroup.readline()
        firstgroup.readline()
        for line in firstgroup:
            headers.add(line.split()[0])
        firstgroup.close()
        groupsize = len(headers)
        #go through all other groups and compare sequence headers
        for secgroup in range(group, sizegroups):
            count = 0
            if gfile != groups[secgroup]:  # only do comparison if needed
                compfile = open(basefolder + groups[secgroup] + "/" + rnd + "hits.txt")
                #get rid of header
                compfile.readline()
                for line in compfile:
                    if line.split()[0] in headers:
                        count += 1
                compfile.close()
            countmatrix[group][secgroup] = float(count)/groupsize

    #group families
    famed = set([])
    fams = []
    for pos, currgroup in enumerate(groups):
        currfam = -1
        #find which family it belongs to, if any
        for famnum, fam in enumerate(fams):
            if groups[pos] in fam:
                currfam = famnum
                break
        if currfam == -1:
            fams.append(set([]))

        #count over row and check if any are over 50%
        for rowpos in range(pos+1, sizegroups):
            if countmatrix[pos][rowpos] > 0.5:
                    fams[currfam].add(groups[rowpos])
                    famed.add(groups[rowpos])
        #count over the row and check if any are over 50%
        for colpos in range(pos):
            if countmatrix[colpos][pos] > 0.5:
                fams[currfam].add(groups[colpos])
                famed.add(groups[colpos])

    #print out families to file
    fout = open(basefolder + "families.txt", 'w')
    for num, fam in enumerate(fams):
        fout.write(str(num+1) + "\t")
        for group in fam:
            fout.write(group + "\t")
        fout.write("\n")
    fout.close()
