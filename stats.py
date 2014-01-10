from cogent import LoadSeqs, RNA
from cogent.parse.stockholm import StockholmParser
from sys import argv
from numpy import zeros

#stats.py /path/to/file.sto /path/to/folder/out/

if __name__ == "__main__":
	if argv[2][-1] != "/":
		argv[2] += "/"
	fin = open(argv[1])
	aln = LoadSeqs(data=StockholmParser(fin).next(), moltype=RNA)
	fin.close()
	consensus = aln.majorityConsensus()
	#counts    A:0 U:1 G:2 C:3 -:4
	counts = zeros(5, dtype=int)
	countsout = open(argv[2] + "counts.txt", 'w')
	countsout.write('\t'.join(['pos','maj','A','T','G','C','-','\n']))
	#count all nucs that do not conform to consensus for each position
	for pos, nucs in enumerate(aln.iterPositions()):
		majnuc = consensus[pos]
		for nuc in nucs:
			if nuc != majnuc:
				if nuc == 'A':
					counts[0] += 1
				elif nuc == 'U':
					counts[1] += 1
				elif nuc == 'G':
					counts[2] += 1
				elif nuc == 'C':
					counts[3] += 1
				elif nuc == '-':
					counts[4] += 1
		countsout.write(str(pos+1) + "\t" + majnuc)
		for count in counts:
			countsout.write('\t' + str(count))
		countsout.write('\n')
		counts = zeros(5, dtype=int)




