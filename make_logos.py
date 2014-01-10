clfrom sys import argv
from os import walk
from subprocess import Popen, PIPE

if __name__ == "__main__":
    basefolder = argv[1]
    if basefolder[-1] != "/":
        basefolder += "/"
    for folder in walk(basefolder).next()[1]:
        if folder == "fasta_groups":
            continue
        #feed sequences to weblogo
        p = Popen(["weblogo", "-F", "pdf", "-A", "rna", "-c", "classic",
            "-t", folder, "-f", basefolder + folder + "/bayesfold-aln.sto",
            "-o", basefolder + folder + "/logo.pdf", "-U", "probability"])
        p.communicate()
