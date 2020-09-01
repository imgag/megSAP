import csv
import gzip
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

reader = csv.DictReader(gzip.open(infile, "rt"), dialect="excel-tab", fieldnames=["#CHROM", "POS", "REF", "ALT", "FATHMM_XF"])
writer = csv.DictWriter(open(outfile, "w"), fieldnames=["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], dialect="excel-tab")

writer.writeheader()

for row in reader:
    writer.writerow({"#CHROM": row["#CHROM"], "POS": row["POS"], "ID": ".", "REF": row["REF"], "ALT": row["ALT"], "QUAL": ".", "FILTER": ".", "INFO": str("FATHMM_XF=" + str(row["FATHMM_XF"]))})
