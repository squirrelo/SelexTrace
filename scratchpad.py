from cogent.struct.rna2d import ViennaStructure
from subprocess import Popen, PIPE
from StringIO import StringIO
from sys import argv

if __name__ == "__main__":    
   basefolder = "/Users/Ely/Desktop/Ely_selection/R7/97percent-lead_clusters-forester_gtLEN/group_"
   for group in range(1,5):
        print "GROUP " + str(group)
        #load in headers
        headers = []
        for line in open(basefolder + str(group) + "/R7hits.txt").readlines()[1:]:
            headers.append(line.split()[0])

        dupefile = open(basefolder + str(group) + "/R7dupes.txt", 'w')

        for secgroup in range(1,5):
            if secgroup == group:
                continue
            dupefile.write("Group " + str(secgroup) + "\n")
            print "Group " + str(secgroup)
            compfile = open(basefolder + str(secgroup) + "/R7hits.txt")
            compfile.readline()
            count = 0
            for line in compfile:
                if line.split()[0] in headers:
                    dupefile.write(line.split()[0] + "\n")
                    count += 1
            print str(count) + " duplicates"
            compfile.close()
