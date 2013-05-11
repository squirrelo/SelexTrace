from sys import argv
from os import walk

if __name__ == "__main__":
    basefolder = argv[1]
    if basefolder[-1] != "/":
        basefolder += "/"
    fileout = open(argv[2], 'w')
    fileout.write("group\tseq count\tunique count\n")
    for folder in walk(basefolder).next()[1]:
        if folder == "fasta_groups":
            continue
        info = open(basefolder + folder + "/log.txt")
        lines = info.readlines()
        info.close()
        count = int(lines[1].split()[0])
        unique = int(lines[2].split()[0])
        fileout.write("%s\t%s\t%s\n" % (folder, count, unique))
    fileout.close()


