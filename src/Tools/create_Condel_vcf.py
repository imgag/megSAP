import gzip
import sys
import csv


infile = sys.argv[1]
outfile = sys.argv[2]

reader = csv.DictReader(gzip.open(infile, "rt"), dialect="excel-tab")
writer = csv.DictWriter(open(outfile, "w"), fieldnames=["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], dialect="excel-tab")

writer.writeheader()

for row in reader:
    writer.writerow({"#CHROM": row["CHR"], "POS": row["START"], "ID": ".", "REF": row["REF"], "ALT": row["ALT"], "QUAL": ".", "FILTER": ".", "INFO": str("CONDEL=" + str(row["CONDEL"]))})
