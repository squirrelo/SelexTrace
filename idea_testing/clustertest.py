from qiime.pick_otus import UclustOtuPicker
from sys import argv


if __name__ == "__main__":
    if argv[2][-1] != "/":
        argv[2] += "/"

    uclust = UclustOtuPicker({'Similarity': float(argv[3]), 'output_dir': argv[2],
        'gapopen': 10.0, 'gapext': 1.0})
    results = uclust(argv[1], 
        log_path=argv[2] + argv[3] + "_clusters_log.txt")
    if results != None:
        print "RESULTS: ", len(results)
        otusout = open(argv[2]+ argv[3] + "_clusters.txt", 'w')
        for result in results:
            otusout.write(str(result) + "\t" + "\t".join(results[result]) + "\n")
        otusout.close()
