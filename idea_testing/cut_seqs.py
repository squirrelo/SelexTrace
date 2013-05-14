from sys import argv
from cogent.parse.fasta import MinimalFastaParser

if __name__ == "__main__":
    out = open(argv[2], 'w')
    for header, seq in MinimalFastaParser(open(argv[1])):
        if seq[:3] == 'GAC':
            if len(seq) < 89:
                print "WARNING: " + header + " shorter than 89nt"
                continue
            seq = seq[:89]
            seq = seq[20:]
            out.write(">%s\n%s\n" % (header,seq))
        elif seq[:3] == 'ACT':
            if len(seq) < 88:
                print "WARNING: " + header + " shorter than 88nt"
                continue
            seq = seq[:88]
            seq= seq[19:]
            out.write(">%s\n%s\n" % (header,seq))
        elif seq[:3] == 'CTT':
            if len(seq) < 87:
                print "WARNING: " + header + " shorter than 87nt"
                continue
            seq = seq[:87]
            seq= seq[18:]
            out.write(">%s\n%s\n" % (header,seq))
        else:
            print "MISMATCH"

