from cogent.parse.infernal_v11 import CmsearchParser
from sys import argv

#count.py /path/to/folder/with/hits #cutoff#

if __name__ == "__main__":
    minscore = float(argv[2])
    for i in range(1,8):
        count = 0
        filein = open(argv[1] + "/R" + str(i) + "hits.txt")
        filein.readline()
        for line in filein:
            lineinfo = line.split(",")
            if float(lineinfo[1]) > minscore:
                count += 1
            else:
                break
        print "R" + str(i) + ": " + str(count)
        filein.close()
