import csv
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

reader = csv.DictReader(open(infile, "rt"), dialect="excel-tab")
writer = csv.writer(open(outfile, "w"), dialect="excel-tab")

for row in reader:
    if row["CHROM"] == 23:
        chrom = "X"
    elif row["CHROM"] == 24:
        chrom = "Y"
    else:
        chrom = row["CHROM"]
    start = int(row["POS"]) - 1
    end = int(row["POS"])
    score = float(row["ABB"])
    writer.writerow([chrom, start, end, score])
