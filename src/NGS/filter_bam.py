import pysam
import argparse
import timeit

start = timeit.default_timer()

# Quality filter flag. 1 if it passes the filter and 0 if not.
def QC_read(read):
	read_out=''
	
	if (read.is_unmapped or read.mapq < args.minMQ or read.cigarstring == None or read.opt("NM") > args.maxMM or (read.cigarstring.count("I") + read.cigarstring.count("D")) > args.maxGAP or not read.is_paired or read.mate_is_unmapped):
		read_out = 0
	else:
		read_out = 1
		
	return (read_out)



### Read parameters
parser = argparse.ArgumentParser(description='Filter bam file based on different values like mismatches, indels, paired reads and mapping quality. The input bam file has to be sorted by read name.')
parser.add_argument('--infile', required=True, dest='infile', help='Input BAM file.')
parser.add_argument('--outfile', required=True, dest='outfile', help='Output BAM file.')
parser.add_argument('--minMQ', required=False, dest='minMQ', type=int, default=30, help='Minimum required mapping quality of aligned read. Default = 30')
parser.add_argument('--maxMM', required=False, dest='maxMM', type=int, default=4, help='Maximum number of mismatches of aligned read . Default = 4')
parser.add_argument('--maxGAP', required=False, dest='maxGAP', type=int, default=1, help='Maximum number of GAPs (indels) . Default = 1')


args = ''
try:
	args = parser.parse_args()
except IOError as io:
	print io
	sys.exit('Error reading parameters.')


### Input BAM
try:
	infile = pysam.Samfile( args.infile, "rb" )
except:
	exit("Cannot open input file.")


### Output BAM
try:
	outfile = pysam.Samfile(args.outfile, mode="wb", template = infile)
except:
	exit("Cannot open output file.")


READNAME = ""
PAIRED_COUNT = 0

### Parse BAM file.
while 1:

	# Get read
	try:
		read = infile.next()
	except StopIteration:
		break
	
	# Grouping and filtering
	if (read.qname != READNAME):
		if (PAIRED_COUNT == 1):
			print "ERROR", READ1.qname
		PAIRED_COUNT = 1
		READNAME = read.qname
		READ1 = read
		
	elif (read.qname == READNAME and PAIRED_COUNT > 0):
		PAIRED_COUNT = 2
		READ2 = read
		
		READ1_out = QC_read(READ1)
		#print READ1.qname, READ1_out
		# If the first read in the paired is already low quality, we do not consider none of them
		if (READ1_out == 0):
			
			
			print "PAIR_REMOVED:", read.qname
			#print read.tlen
			#print READ1.is_unmapped, READ1.mapq,READ1.cigarstring, str(READ1.is_paired)
			#print READ2.is_unmapped, READ2.mapq,READ2.cigarstring, str(READ2.is_paired)
			#print "\n"
			
			READ1_out = ""
			READ2_out = ""
			
			PAIRED_COUNT = 0
			READ1 = ""
			READ2 = ""			
			
		else:
			
			READ2_out = QC_read(READ2)
			#print READ2.qname, READ2_out
			
			if (READ2_out == 1):
				outfile.write(READ1)
				outfile.write(READ2)
			else:
				print "PAIR_REMOVED:", read.qname
			
			READ1_out = ""
			READ2_out = ""
			
			READ1 = ""
			READ2 = ""

			PAIRED_COUNT = 0


print 'TIME'
stop = timeit.default_timer()
print stop - start

### Finish
infile.close()
outfile.close()
exit(0)



### pysam documentation
# http://pysam.readthedocs.org/en/latest/
