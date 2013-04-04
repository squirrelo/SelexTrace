if __name__ == "__main__":
	from cogent.app.infernal import cmsearch_from_file
	from cogent import LoadSeqs, RNA
	from math import log
	from os import mkdir
	from os.path import isdir
	
	for i in range(1,4):
		print "making folder structural/Construct" + str(i)
		if not isdir("structural/Construct" + str(i)):
			mkdir("structural/Construct" + str(i))
		for j in range(1,8):
			print "Loading R"+ str(j) + "_CleanedStrippedUnique.fas"
			seqs = LoadSeqs("R" + str(j) + "_CleanedStrippedUnique.fas", \
			moltype=RNA, aligned=False, format='fasta', label_to_name=lambda x: x.split()[0])
			print "Running Infernal for structural/Construct" + str(i) + "/R" + str(j) + "_tabfile.txt"
			cmsearch_result = cmsearch_from_file("structural/Construct" + str(i) + ".cm", seqs, RNA, \
			params={"--tabfile": "structural/Construct" + str(i) + "/R" + str(j) + "_tabfile.txt"})
			print "result: " + str(cmsearch_result)
			
			#cutoff=log(len(seqs),2) ??
	