import csv
import gzip
import sys

chroms = str(sys.argv[1]).split(",")
outfile = sys.argv[2]

writer = csv.DictWriter(open(outfile, "w"), fieldnames=["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], dialect="excel-tab")

writer.writeheader()

for chrom in chroms:
    reader = csv.DictReader(open(chrom, "rt"), dialect="excel")
    for row in reader:
        mutation = str(row["Mutation"]).split(",")
        writer.writerow({"#CHROM": mutation[1], "POS": mutation[2], "ID": ".", "REF": mutation[3], "ALT": mutation[4], "QUAL": ".", "FILTER": ".", "INFO": str("MutationAssessor=" + str(row["FI score"]))})
