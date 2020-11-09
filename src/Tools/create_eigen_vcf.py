import csv
import gzip
import sys


## uncomment commented code lines if non-coding information is needed


chroms = str(sys.argv[1]).split(",")
outfile_coding = sys.argv[2]
#outfile_noncoding = sys.argv[3]

coding = csv.writer(open(outfile_coding, "w"), delimiter="\t")
#noncoding = csv.writer(open(outfile_noncoding, "w"), delimiter="\t")
coding.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
#noncoding.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])

for chrom in chroms:
    with gzip.open(chrom, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            else:
                split_line = line.split("\t")
                if split_line[4] == "c":
                    coding.writerow([split_line[0], split_line[1], ".", split_line[2], split_line[3], ".", ".", str("EIGEN_PHRED=" + split_line[6])])
                #elif split_line[4] == "n":
                #    noncoding.writerow([split_line[0], split_line[1], ".", split_line[2], split_line[3], ".", ".", str("EIGEN_PHRED=" + split_line[6])])
