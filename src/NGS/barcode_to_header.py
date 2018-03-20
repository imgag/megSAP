#!/usr/bin/env python
import argparse
import gzip
from itertools import izip

# from Bio import bgzf

parser = argparse.ArgumentParser(description='Getting barcodes in fastq.gz file and labels')
parser.add_argument('-i', '--infile', type=str, help='Input fastq(.gz)', required=True)
parser.add_argument('-bc1', '--bc1', type=str, help='Barcode 1 fastq(.gz) file', default='')
parser.add_argument('-bc2', '--bc2', type=str, help='Barcode 2 fastq(.gz) file', default='')
parser.add_argument('-o', '--outfile', type=str, help='Prefix of the output file', required=True)

args = parser.parse_args()

FILE = args.infile
bc_file1 = args.bc1
bc_file2 = args.bc2
out = args.outfile

# out_seq = out + ".fastq"
# OUT_SEQ = open(out_seq,'w')
out_seq = out
OUT_SEQ = gzip.open(out_seq, 'w')

if (len(bc_file1) == 0 and len(bc_file2) == 0):
    print("Error:\nNo fastq files with barcodes provided")
    exit()

L = 0

tmp = open(FILE, 'r')
magic_number = tmp.read(2)
tmp.close()

# Two barcode provided
if (len(bc_file1) > 0 and len(bc_file2) > 0):
    with open(FILE) if magic_number != '\x1f\x8b' else gzip.open(FILE) as f1:
        with open(bc_file1) if magic_number != '\x1f\x8b' else gzip.open(bc_file1) as f2, open(
                bc_file2) if magic_number != '\x1f\x8b' else gzip.open(bc_file2) as f3:
            for (c1, c2, c3) in izip(f1, f2, f3):
                line1 = c1.rstrip('\n')
                line2 = c2.rstrip('\n')
                line3 = c3.rstrip('\n')

                L = L + 1

                INFO1 = line1.split(" ")
                INFO2 = line2.split(" ")
                INFO3 = line3.split(" ")

                if line1.startswith('@') and line2.startswith('@') and line3.startswith('@') and L == 1:
                    read = line1
                    READ = line1.split(" ")
                    READ_NAME = READ[0]
                    READ_NAME2 = line2.split(" ")[0]
                    READ_NAME3 = line3.split(" ")[0]

                if READ_NAME != READ_NAME2 and READ_NAME != READ_NAME3:
                    exit(
                        "Read names differ: " + line1 + ", " + line2 + ", " + line3 + ". Check order of reads in input files.")

                elif (L == 2):
                    Seq = line1

                    # We ignore last base of the barcode due to low reliability
                    BARCODE1 = line2
                    BARCODE1_l = len(BARCODE1)

                    # We ignore last base of the barcode due to low reliability
                    BARCODE2 = line3
                    BARCODE2_l = len(BARCODE2)

                elif (L == 3):
                    description = line1

                elif (L == 4):
                    L = 0
                    QUAL = line1

                    ### READS FILE
                    # 1. Read name "@NAME:bc1_length,bc2_length:bc1_seq,bc2_seq
                    LINE_1 = READ_NAME + ":" + str(BARCODE1_l) + "," + str(
                        BARCODE2_l) + ":" + BARCODE1 + "," + BARCODE2 + "\n"

                    # 2. Sequence without the barcodes
                    LINE_2_r = Seq + "\n"

                    # 3. Description of read
                    LINE_3_r = description + "\n"

                    # 4. Qualities
                    LINE_4_r = QUAL + "\n"

                    # SAVING READ INFO in fastq file
                    OUT_SEQ.write(LINE_1)
                    OUT_SEQ.write(LINE_2_r)
                    OUT_SEQ.write(LINE_3_r)
                    OUT_SEQ.write(LINE_4_r)

                else:
                    print
                    "error"
                    print
                    line1
                    print
                    line2
    f1.close()
    f2.close()
    f3.close()

# Only barcode 1 provided    
elif (len(bc_file1) > 0 and len(bc_file2) == 0):
    with open(FILE) if magic_number != '\x1f\x8b' else gzip.open(FILE, 'r') as f1:
        with open(bc_file1) if magic_number != '\x1f\x8b' else gzip.open(bc_file1, 'r') as f2:
            for (c1, c2) in izip(f1, f2):
                line1 = c1.rstrip('\n')
                line2 = c2.rstrip('\n')

                L = L + 1

                INFO1 = line1.split(" ")
                INFO2 = line2.split(" ")

                if line1.startswith('@') and line2.startswith('@') and L == 1:
                    READ = line1.split(" ")
                    READ_NAME = READ[0]
                    READ_NAME2 = line2.split(" ")[0]
                    if READ_NAME != READ_NAME2:
                        exit("Read names differ: " + line1 + ", " + line2 + ". Check order of reads in input files.")
                    read = line1

                elif (L == 2):
                    Seq = line1

                    # We ignore last base of the barcode due to low reliability
                    BARCODE1 = line2
                    BARCODE1_l = len(BARCODE1)

                    # We ignore last base of the barcode due to low reliability
                    BARCODE2 = "."
                    BARCODE2_l = 0

                elif (L == 3):
                    description = line1

                elif (L == 4):
                    L = 0
                    QUAL = line1

                    ### READS FILE
                    # 1. Read name "@NAME:bc1_length,bc2_length:bc1_seq,bc2_seq
                    LINE_1 = READ_NAME + ":" + str(BARCODE1_l) + "," + str(
                        BARCODE2_l) + ":" + BARCODE1 + "," + BARCODE2 + "\n"

                    # 2. Sequence without the barcodes
                    LINE_2_r = Seq + "\n"

                    # 3. Description of read
                    LINE_3_r = description + "\n"

                    # 4. Qualities
                    LINE_4_r = QUAL + "\n"

                    # SAVING READ INFO in fastq file
                    OUT_SEQ.write(LINE_1)
                    OUT_SEQ.write(LINE_2_r)
                    OUT_SEQ.write(LINE_3_r)
                    OUT_SEQ.write(LINE_4_r)

                else:
                    print
                    "error"
                    print
                    line1
                    print
                    line2
    f1.close()
    f2.close()

# Only barcode 2 provided
if (len(bc_file1) == 0 and len(bc_file2) > 0):
    with open(FILE) if magic_number != '\x1f\x8b' else gzip.open(FILE) as f1:
        with open(bc_file2) if magic_number != '\x1f\x8b' else gzip.open(bc_file2) as f3:
            for (c1, c3) in izip(f1, f3):
                line1 = c1.rstrip('\n')
                line3 = c3.rstrip('\n')

                L = L + 1

                INFO1 = line1.split(" ")
                INFO3 = line3.split(" ")

                if line1.startswith('@') and line3.startswith('@') and L == 1:
                    read = line1
                    READ = line1.split(" ")
                    READ_NAME = READ[0]
                    READ_NAME3 = line3.split(" ")[0]
                    if READ_NAME != READ_NAME3:
                        exit("Read names differ: " + line1 + ", " + line3 + ". Check order of reads in input files.")

                elif (L == 2):
                    Seq = line1

                    # We ignore last base of the barcode due to low reliability
                    BARCODE1 = "."
                    BARCODE1_l = 0

                    # We ignore last base of the barcode due to low reliability
                    BARCODE2 = line3
                    BARCODE2_l = len(BARCODE2)

                elif (L == 3):
                    description = line1

                elif (L == 4):
                    L = 0
                    QUAL = line1

                    ### READS FILE
                    # 1. Read name "@NAME:bc1_length,bc2_length:bc1_seq,bc2_seq
                    LINE_1 = READ_NAME + ":" + str(BARCODE1_l) + "," + str(
                        BARCODE2_l) + ":" + BARCODE1 + "," + BARCODE2 + "\n"

                    # 2. Sequence without the barcodes
                    LINE_2_r = Seq + "\n"

                    # 3. Description of read
                    LINE_3_r = description + "\n"

                    # 4. Qualities
                    LINE_4_r = QUAL + "\n"

                    # SAVING READ INFO in fastq file
                    OUT_SEQ.write(LINE_1)
                    OUT_SEQ.write(LINE_2_r)
                    OUT_SEQ.write(LINE_3_r)
                    OUT_SEQ.write(LINE_4_r)

                else:
                    print
                    "error"
                    print
                    line1
                    print
                    line2
    f1.close()
    f3.close()

OUT_SEQ.close()
