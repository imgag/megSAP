import pysam
import argparse
from collections import defaultdict
import re
import copy

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def Rcount(infile, MIPS_DICT):
	## Getting pair info
	for read1, read2 in read_pair_generator(infile):
		
		# Coordinates read 1
		CHROM1 = infile.getrname(read1.tid)
		CHROM1 = str(re.sub('^chr', '', CHROM1))
		
		Start1 = read1.reference_start + 1 # To solve the 0 based coordinate
		End1 = read1.reference_end
		
		# Coordinates read 2
		CHROM2 = infile.getrname(read2.tid)
		CHROM2 = str(re.sub('^chr', '', CHROM2))
		
		Start2 = read2.reference_start + 1 # To solve the 0 based coordinate
		End2 = read2.reference_end
		
		# Min Start and max end
		Start = str(min(Start1, Start2))
		End = str(max(End1, End2))
		
		if (CHROM1 == CHROM2):
			Code = ';'.join([CHROM1, Start, End])
			
			if not Code in MIPS_DICT:
				MIPS_DICT['NA;NA;NA'][1] += 1
			else:
				MIPS_DICT[Code][1] += 1
			
	return MIPS_DICT

### Read parameters
parser = argparse.ArgumentParser(description='Filter bam file based on different values like mismatches, indels, paired reads and mapping quality. The input bam file has to be sorted by read name.')
parser.add_argument('--infile1', required=True, dest='infile1', help='Input BAM file before deduplication.')
parser.add_argument('--infile2', required=True, dest='infile2', help='Input BAM file after deduplication.')
parser.add_argument('--mips', required=True, dest='mips', help='MIPs design')
parser.add_argument('--outfile', required=True, dest='outfile', help='Output BAM file.')


args = ''
try:
	args = parser.parse_args()
except IOError as io:
	print io
	sys.exit('Error reading parameters.')


### Input BAM
try:
	infile1 = pysam.AlignmentFile(args.infile1, 'rb')
	infile2 = pysam.AlignmentFile(args.infile2, 'rb')
except:
	exit("Cannot open input file.")


## Saving mip probes in a dictionary
MIPS_DICT = {}
MIPS_DICT['NA;NA;NA'] = ['Others', 0]

mips=open(args.mips,'r')
for line in mips:
	line = line.rstrip('\n')
	
	if not line.startswith('>'):
		probes = line.split("\t")
		
		# Mips ID
		ID = probes[-1] # If we want mip name, we have to change this column by the last columnq
		CHROM = probes[2]
		CHROM = str(re.sub('^chr', '', CHROM))
		
		# We assume always the same input format (same column numbers)
		START = probes[11]
		END = probes[12]
		
		CODE = ';'.join([CHROM, START, END])
		
		
		if not CODE in MIPS_DICT: # If you want to speed up a bit this mini tool, you can create a dictionary of dictionary spliting by chromosomes
			MIPS_DICT[CODE] = [ID, 0]

# One dictionary for before and other for after
MIPS_DICT_before_dedup = copy.deepcopy(MIPS_DICT)
MIPS_DICT_dedup = copy.deepcopy(MIPS_DICT)		

## Running functions in both bam files	
BEFORE = Rcount(infile1, MIPS_DICT_before_dedup)
AFTER = Rcount(infile2, MIPS_DICT_dedup)

## Printing counts
outfile = open(args.outfile,'w')
HEADER = ['Mip_name', 'Chromosome', 'Start', 'End', 'Fragments_before_deduplication', 'Fragments_after_deduplication']
HEADER = '\t'.join(HEADER)+'\n'
outfile.write(HEADER)

for mip in MIPS_DICT:
	CHROM, START, END = mip.split(";")
	
	# Counting before deduplication
	ID_bef, Count_bef = BEFORE[mip]
	
	# Counting after deduplication
	ID_aft, Count_aft = AFTER[mip]

	LINE = [str(ID_bef), str(CHROM), str(START), str(END), str(Count_bef), str(Count_aft)]
	LINE = '\t'.join(LINE)+'\n'
	outfile.write(LINE)

### Finish
infile1.close()
infile2.close()
outfile.close()
exit(0)



### pysam documentation
# http://pysam.readthedocs.org/en/latest/
