from sys import argv
from os import walk
from os.path import exists
from numpy import zeros

#overlap_check.py basefolder #round

if __name__ == "__main__":
    basefolder = argv[1]
    if basefolder[:-1] != "/":
        basefolder += "/"

    rnd = "R" + argv[2]

    #get total number fo groups in the round, add 1 for looping
    sizegroups = len(walk(basefolder).next()[1]) + 1
    print str(sizegroups-1) + " groups"

    #create csv file for data
    overlapfile = open(basefolder + "overlap.csv", 'w')
    #empty first cell of header
    overlapfile.write(",")
    #write out all group names for column
    for x in range(1, sizegroups):
        overlapfile.write("group" + str(x) + ",")
    overlapfile.write("\n")

    countmatrix = zeros(shape=(sizegroups-1, sizegroups-1), dtype=int)
    #compare all group sequences to other sequences
    #starting at 1 for folder name ease
    for group in range(1, sizegroups):
        #write out groupname to csv file as row indicator
        overlapfile.write("group" + str(group) + ",")

        #check if group was run in infernal
        #if so, compare seq counts over all groups
        if exists(basefolder + "group_" + str(group) + "/" + rnd + "hits.txt"):
            #write out the dashes needed to line up the matrix in file
            for x in range(1, group):
                overlapfile.write("-,")
            #load in headers of current group
            headers = set([])
            firstgroup = open(basefolder + "group_" + str(group) + "/" + rnd + "hits.txt")
            #get rid of header
            firstgroup.readline()
            firstgroup.readline()
            for line in firstgroup:
                headers.add(line.split()[0])
            firstgroup.close()
            dupefile = open(basefolder + "group_" + str(group) + "/" + rnd + "dupes.txt", 'w')
            #go through all other groups and compare sequence headers
            for secgroup in range(group, sizegroups):
                if exists(basefolder + "group_" + str(secgroup) + "/" + rnd + "hits.txt"):
                    count = 0
                    if group != secgroup:  # only do comparison if needed
                        compfile = open(basefolder + "group_" + str(secgroup) + "/" + rnd + "hits.txt")
                        dupefile.write("Group " + str(secgroup) + "\n")
                        #get rid of header
                        compfile.readline()
                        for line in compfile:
                            if line.split()[0] in headers:
                                dupefile.write(line.split()[0] + "\n")
                                count += 1
                        compfile.close()
                    else:  # we are comparing group to itself so only need seq count
                        count = len(headers)
                    countmatrix[group-1][secgroup-1] = count
                    overlapfile.write(str(count) + ",")
                else:
                    overlapfile.write("-,")
                    countmatrix[group-1][secgroup-1] = -1
        else:
            for secgroup in range(1, sizegroups):
                overlapfile.write("-,")
                countmatrix[group-1][secgroup-1] = -1
        overlapfile.write("\n")
    overlapfile.write("\n")

    #now create percentages matrix
    #dealing with matrix starting at 0: change sizegroups accordingly
    sizegroups -= 1
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
