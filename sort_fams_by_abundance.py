from sys import argv

#sort_fams_by_abundance.py  /path/to/group_sizes.txt /path/to/families.txt

if __name__ == "__main__":
    abundanceorder = open(argv[1])
    #get rid of header
    abundanceorder.readline()
    aborder = []
    for line in abundanceorder:
        aborder.append(line.split("\t")[0])
    abundanceorder.close()
    famsin = open(argv[2])
    #create list of families, each having list of groups in fam
    fams = []
    for line in famsin:
        currfam = line.strip().split("\t")
        fams.append(currfam)
    #sort each family by abundance in pool round
    for fam in fams:
        fam.sort(key=(aborder+fam).index)
    #sort all families by which contains most abundant sequences
    groupssorted = set([])
    sortedfams = []
    for group in aborder:
        if group in groupssorted:
            continue
        for fam in fams:
            if group in fam:
                sortedfams.append(fam)
                groupssorted |= set(fam)
                break
    #print out families now that they are sorted
    outpath = argv[1][:argv[1].rfind("/")]
    outpath += "/families_sorted.txt"
    outfile = open(outpath, 'w')
    for fam in sortedfams:
        outfile.write('\t'.join(fam) + "\n")
    outfile.close()
