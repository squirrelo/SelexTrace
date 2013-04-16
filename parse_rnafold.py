from remove_duplicates import remove_duplicates
from subprocess import Popen, PIPE
from sys import argv
from time import time
from cogent.app.vienna_package import RNAfold
from cogent.parse.rnaforester import cluster_parser


def parse_RNAfold(lines):
    ret = {}
    hold = []
    first = True
    for line in lines:
        line = line.strip()
        if line[0] == ">" and not first:
            structinfo = hold[2].split()
            ret[hold[0][1:]] = (hold[1], structinfo[0], 
                float(structinfo[-1].strip("(").strip(")")))
            hold = [line]
        else:
            hold.append(line)
            first = False

    return ret


def run_RNAfold(filein, temperature=37):
    r = RNAfold(InputHandler='_input_as_path', WorkingDir="/tmp", HALT_EXEC=True)
    r.Parameters['-T'].on(temperature)
    res = r(filein)
    lines = res['StdOut'].readlines()
    res.cleanUp()
    return parse_RNAfold(lines.split("\n"))

def score_rnaforester(struct1, struct2):
    '''returns relative score of two structures'''
    #return gigantically negative number if no structure for one struct
    if "(" not in struct1 or "(" not in struct2:
        raise ValueError(struct1 + "\n" + struct2 + "\nNo pairing in given structures!")
    p = Popen(["RNAforester", "--score", "-r"], stdin=PIPE, stdout=PIPE)
    p.stdin.write(''.join([struct1, "\n", struct2, "\n&"]))
    return float(p.communicate()[0].split("\n")[-2])


def group_by_forester(filein, foresterscore):
    '''run rnaforester clustering and return list of cluster groups'''
    p = Popen(["RNAforester", "--score", "-m", "-r", "-mt=" + str(foresterscore), "-f", filein], stdout=PIPE)
    lines = p.communicate()[0].split("\n")
    structclust = []
    for group in cluster_parser(lines):
        numclust = int(group[2].split(':')[1])
        cluster = []
        for n in range(numclust):
            cluster.append(group[(4+n)].split()[0])
        structclust.append(cluster)
    return structclustorgin

if __name__ == '__main__':
    start = time()
    filein = argv[1].strip(" '\"")
    #initial = run_RNAfold(filein, 25)
    initial = parse_RNAfold(open(filein, 'U'))
    structs = open("structs.fasta", 'w')
    for header in initial:
        structs.write(">" + header + "\n" + initial[header][1] + "\n")
    structs.close()
    structgroups = group_by_forester("structs.fasta", 0.6)
    print "STRUCTGROUPS:", len(structgroups)
    structs = open("struct_clusters.fasta", 'w')
    for groupnum in range(len(structgroups)):
        print len(structgroups[groupnum])
        structs.write('\t'.join(structgroups[groupnum]) + "\n")
    structs.close()
    print "TIME:", (time()-start)/60, "m"
