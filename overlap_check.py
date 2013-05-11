from sys import argv
from os import walk
from os.path import exists
from os import walk
from numpy import zeros

#overlap_check.py basefolder #round

if __name__ == "__main__":
    basefolder = argv[1]
    if basefolder[:-1] != "/":
        basefolder += "/"

    rnd = "R" + argv[2]


    #create csv file for data
    overlapfile = open(basefolder + rnd + "overlap.csv", 'w')
    #empty first cell of header
    overlapfile.write(",")

    #get all group names
    groups = [g for g in walk(basefolder).next()[1]]
    if "fasta_groups" in groups:
        groups.remove("fasta_groups")
    sizegroups = len(groups)
    #get total number fo groups in the round, add 1 for looping
    print str(sizegroups-1) + " groups"
    #write out all group names for column
    for x in range(sizegroups):
        overlapfile.write(groups[x] + ",")
    overlapfile.write("\n")
    #compare all group sequences to other sequences
    countmatrix = zeros(shape=(sizegroups, sizegroups), dtype=int)
    for group, gfile in enumerate(groups):
        #write out groupname to csv file as row indicator
        overlapfile.write(gfile + ",")

        #compare seq counts over all groups
        for x in range(group):
            overlapfile.write("-,")
        #load in headers of current group
        headers = set([])
        firstgroup = open(basefolder + gfile + "/" + rnd + "hits.txt")
        #get rid of header
        firstgroup.readline()
        firstgroup.readline()
        for line in firstgroup:
            headers.add(line.split()[0])
        firstgroup.close()
        dupefile = open(basefolder + gfile + "/" + rnd + "dupes.txt", 'w')
        #go through all other groups and compare sequence headers
        for secgroup in range(group, sizegroups):
            count = 0
            if gfile != groups[secgroup]:  # only do comparison if needed
                compfile = open(basefolder + groups[secgroup] + "/" + rnd + "hits.txt")
                dupefile.write(groups[secgroup] + "\n")
                #get rid of header
                compfile.readline()
                for line in compfile:
                    if line.split()[0] in headers:
                        dupefile.write(line.split()[0] + "\n")
                        count += 1
                compfile.close()
            else:  # we are comparing group to itself so only need seq count
                count = len(headers)
            countmatrix[group][secgroup] = count
            overlapfile.write(str(count) + ",")
        overlapfile.write("\n")
    overlapfile.write("\n")

    #now create percentages matrix
    for group in range(sizegroups):
        overlapfile.write("group" + str(group+1) + ",")
        if countmatrix[group][group] != -1:
             #write out the dashes needed to line up the matrix in file
            for x in range(group):
                overlapfile.write("-,")

            groupcount = float(countmatrix[group][group])
            for secgroup in range(group, sizegroups):
                if countmatrix[group][secgroup] != -1:
                    overlapfile.write(str(countmatrix[group][secgroup]/groupcount*100) + ",")
                else:
                    overlapfile.write("-,")
        else:
            for secgroup in range(sizegroups):
                overlapfile.write("-,")
        overlapfile.write("\n")
