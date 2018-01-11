#!/usr/bin/env python

import argparse
import pysam
from difflib import SequenceMatcher
import numpy


## EXAMPLE COMMAND LINE
# python correct_bam_barcodes_21.11.2016.py --infile PATH/TO/test.bam --outfile PATH/TO/corrected.test.bam --barcodes BOTH

def similarity(a, b):
    ratio = SequenceMatcher(None, a, b).ratio()
    return ratio

def keywithmaxval(d):
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]

def extract_barcode(READ,BC_TYPE): # Individual read. It returns the barcode of the read
	NAME = str(READ.qname)
	BC = NAME.split(':')[-1]
	BC1 = BC.split(',')[0]
	BC2 = BC.split(',')[1]
	
	# Grouping by barcodes
	if (BC_TYPE == "BEGINNING"):
		bc = BC1

	elif (BC_TYPE == "END"):
		bc = BC2
	    
	elif (BC_TYPE == "BOTH"):
		bc = BC1 + BC2
	#print bc
	return(bc)

def extract_bc_groups(LIST,BC_TYPE): # Input, list of reads that start are end are equal
	BARCODE_DICT={}
	#BARCODE = ""
	BC={}
		
	for i in LIST: # i represents each read of the list
		bc = extract_barcode(i,BC_TYPE) # Extract the barcode
		
		if bc in BC:
			BC[bc].append(i)
	   
		else:
			BC[bc] = list()
			BC[bc].append(i)
	
		
	while len(BC) > 0:
		MAX = "".join(([k for k in BC.keys() if BC[k]==max(BC.values(),key=len)])) # For getting the key (Barcode) with more number of reads
		
		SIM = [BC.keys()[x] for x in range(0,len(BC.keys())) if similarity(BC.keys()[x], MAX) >= ((len(MAX)-Errors)/float(len(MAX)))] # Getting similar barcodes
		
		BARCODE_DICT[MAX] = list() # Creating key in the hash based on our barcode where reads of this group will be saved
	
		for i in SIM:
			BARCODE_DICT[MAX].extend(BC[i]) # Grouping reads based on similarity of the barcodes
			del BC[i] # Removing barcodes already considered
	
	
	return BARCODE_DICT # Dictionary with Barcode as a key and reads as values
    

def ReverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])
    
# Reduce mapping quality to < 20
def reduce_mapq(read):
    if read.mapq >= 20:
	read.mapq = 19
	NEW_READ = read
    else:
	NEW_READ = read
	
    return (NEW_READ)

def ascii2dec(ASCII):    
    qual = [ord(i)-33 for i in ASCII]    
    return qual


def most_common_base(LIST, QUALITIES, minBQ,STEP):
    QUALITIES_NATIVE = get_native_qual(QUALITIES,STEP)
    
    HQ = [x for x in range(0,len(QUALITIES_NATIVE)) if QUALITIES_NATIVE[x] >= minBQ]
    
    result_list = [LIST[i] for i in HQ]    
    result_qual = [QUALITIES[i] for i in HQ] 
    
    BASE = max(set(result_list), key=result_list.count)
    NUM_DIFF_BASES = len(set(result_list))
    
    DIFF_COUNT = len([x for x in range(0,len(result_list)) if result_list[x] != BASE])
    BASE_COUNT = len([x for x in range(0,len(result_list)) if result_list[x] == BASE])
    
    return BASE,NUM_DIFF_BASES,BASE_COUNT,DIFF_COUNT

def most_common_base_low_qual(LIST):
    BASE = max(set(LIST), key=LIST.count)
    NUM_DIFF_BASES = len(set(LIST))
    
    return BASE,NUM_DIFF_BASES


def get_qualities(BASE,BASES,QUALITIES):
	QUAL = [(ord(QUALITIES[x])-33) for x in range(0,len(BASES)) if BASES[x] == BASE]
	#LIST = [BASES[x] for x in range(0,len(BASES)) if BASES[x] == BASE]
	return QUAL
'''
def get_native_qual(QUALITIES):
    NATIVE_QUAL=[]
    for i in QUALITIES:
        if i > 50:
            Q = (i-50)*10
        else:
            Q = i   
        
        NATIVE_QUAL.append(Q)
    return(NATIVE_QUAL)
'''
def get_native_qual(QUALITIES,STEP):
    NATIVE_QUAL=[]
    #print QUALITIES
    if (STEP == 1 or STEP == 0):
        for i in QUALITIES:
            if i > 50:
                Q = (i-50)*10
            else:
                Q = i   
            
            NATIVE_QUAL.append(Q)
            
    elif (STEP == 2):
        NATIVE_QUAL = QUALITIES
    
    return(NATIVE_QUAL)


def one_duplicate_qual(read, STEP, minBQ):
    #READ = reduce_mapq(read)
    READ = read
    qualities = READ.qual
    QUALITIES = ascii2dec(qualities)
    QUAL=[]
    QUALITIES_NATIVE = get_native_qual(QUALITIES,STEP)
    if (STEP == 1):
        for i in range(0,len(QUALITIES_NATIVE)):
            if (QUALITIES_NATIVE[i] >= minBQ):
                QUAL.append(str(chr(10 + 33)))
            else:
                QUAL.append("!") # 0 quality = chr(0 + 33)           
    if (STEP == 2):
        for i in range(0,len(QUALITIES_NATIVE)):
		if (QUALITIES_NATIVE[i] >= 10):
			MEDIAN = int(round(QUALITIES_NATIVE[i]/10,0))
			# Fixing max value
			if (MEDIAN > 5):
				MEDIAN = 5				
			QUALi = 10 + MEDIAN
		else:
			QUALi = 0
		
		QUAL.append(str(chr(QUALi + 33)))

    if (STEP == 0):
        for i in range(0,len(QUALITIES_NATIVE)):
		if (QUALITIES_NATIVE[i] >= minBQ):
			MEDIAN = int(round(QUALITIES_NATIVE[i]/10,0))
			# Fixing max value
			if (MEDIAN > 5):
				MEDIAN = 5				
			QUALi = 10 + MEDIAN
		else:
			QUALi = 0
		
		QUAL.append(str(chr(QUALi + 33)))
		
    QUAL = ''.join(QUAL)
    #print STEP
    #print QUAL
    #print QUALITIES_NATIVE
    
    READ.qual = QUAL	   
    return(READ)


'''
def one_duplicate_qual(qualities, STEP, minBQ):
    QUALITIES = ascii2dec(qualities)
    QUAL=[]
    QUALITIES_NATIVE = get_native_qual(QUALITIES,STEP)
    if (STEP == 1):
        for i in range(0,len(QUALITIES_NATIVE)):
            if (QUALITIES_NATIVE[i] >= minBQ):
                QUAL.append(str(chr(10 + 33)))
            else:
                QUAL.append("!") # 0 quality = chr(0 + 33)           
    if (STEP == 2):
        for i in range(0,len(QUALITIES_NATIVE)):
		if (QUALITIES_NATIVE[i] >= 10):
			MEDIAN = int(round(QUALITIES_NATIVE[i]/10,0))
			# Fixing max value
			if (MEDIAN > 5):
				MEDIAN = 5				
			QUALi = 10 + MEDIAN
		else:
			QUALi = 0
		
		QUAL.append(str(chr(QUALi + 33)))
		
    QUAL = ''.join(QUAL)
    #print STEP
    #print QUAL
    #print QUALITIES_NATIVE
    	   
    return(QUAL)
'''

def overlap_quality_base_check(QUALITIES,minBQ):
    return (len([x for x in range(0,len(QUALITIES)) if QUALITIES[x] >= 50+(minBQ/10)]))

def low_quality_base_check(QUALITIES,minBQ,STEP):
	QUALITIES_NATIVE = get_native_qual(QUALITIES,STEP)
	return (len(QUALITIES_NATIVE) >  len([x for x in range(0,len(QUALITIES_NATIVE)) if QUALITIES_NATIVE[x] < minBQ]))
    
def consensus_quality(QUALITIES,minBQ,ERRORS,STEP):
    QUALITIES_NATIVE = get_native_qual(QUALITIES,STEP)
    
    # How many bases with high good quality
    COPIES = len([x for x in range(0,len(QUALITIES_NATIVE)) if QUALITIES_NATIVE[x] >= minBQ])
    
    # How many bases with high quality have overlap correction
    OVERLAPPED = overlap_quality_base_check(QUALITIES,minBQ)
    	
    # Consider more or equal than 5 as best agreement
    if (COPIES >= 5):
        NEW_COPIES = 5
    else:
        NEW_COPIES = COPIES    
    
    
    if (ERRORS > 5):
	ERRORS = 5

	# Doble check
    MAX_QUAL = max(QUALITIES_NATIVE)
    if (MAX_QUAL < minBQ):
        NEW_QUAL = 0
    else:
	# Base qualities are correlated with how many reads with good quality are showing the same allele and how many of them had errors
	# COPIES * 10 + ERRORS
        if (STEP == 1):
		NEW_QUAL = NEW_COPIES * 10 + int(ERRORS)
	
        if (STEP == 2 or STEP == 0):
		if ((ERRORS >= 1 and float(ERRORS)/(COPIES+ERRORS) > 0.25) or ERRORS >= 3):
			NEW_QUAL = 0
		else:	
			MEDIAN = int(round((numpy.mean(QUALITIES_NATIVE))/10,0))
			
			# Fixing max value
			if (MEDIAN > 5):
				MEDIAN = 5
				
			NEW_QUAL = NEW_COPIES * 10 + MEDIAN
		
	#NEW_QUAL = NEW_COPIES * 10 + OVERLAPPED
    
    #print len(QUALITIES), len(QUALITIES_NATIVE), NEW_QUAL, NEW_COPIES, OVERLAPPED, MAX_QUAL, max(QUALITIES)
    
    return (NEW_QUAL)
    
# Function to correct reads grouped by barcodes
def GET_FINAL_READ(reads,minBQ,STEP):
    
    CONSENSUS_SEQ = list()
    CONSENSUS_QUAL = list()
    
    # Getting the amount 
    #print len(reads)    
    
    if (len(reads) <= 1):
	for i in reads:
	    # Encoding quality
	    #SINGLE_QUAL = i.qual
	    #RECODE_QUAL = one_duplicate_qual(SINGLE_QUAL, STEP, minBQ)
	    #i.qual = RECODE_QUAL
	    #CONSENSUS_READ = i
	    CONSENSUS_READ = one_duplicate_qual(i, STEP, minBQ)
    else:
	SEQ = {}
	QUAL = {}
	MAPQ = list()
	
	for i in reads:
	    
	    LEN = i.rlen
	    #print LEN
	    seq = i.seq
	    qual = i.qual
	    MAPQ.append(i.mapq)
	    
	    for b in range(0,LEN):
		
		BASE = seq[b]
		BASE_QUAL = qual[b]
		
		

		if b in SEQ:
		    SEQ[b].append(BASE)
		    QUAL[b].append(BASE_QUAL)
		    
		else:
		    SEQ[b] = list()
		    QUAL[b] = list()
		    SEQ[b].append(BASE)
		    QUAL[b].append(BASE_QUAL)		
		
		#if ((ord(BASE_QUAL)-33) >= minQUAL):
		#    if b in SEQ:
		#	SEQ[b].append(BASE)
		#	QUAL[b].append(BASE_QUAL)
		#    
		#    else:
		#	SEQ[b] = list()
		#	QUAL[b] = list()
		#	SEQ[b].append(BASE)
		#	QUAL[b].append(BASE_QUAL)
	
	
	for position in SEQ:
		CONSENSUS_BASE = ""
		NUM_DIFF_BASES = ""
		CONSENSUS_BASE_COUNT = ""
		DIFF_COUNT = ""

		Q = ascii2dec(QUAL[position])
		if (low_quality_base_check(Q,minBQ,STEP)):
			BASE = most_common_base(SEQ[position],Q,minBQ,STEP)
		
			CONSENSUS_BASE = BASE[0]
			NUM_DIFF_BASES = BASE[1]
			CONSENSUS_BASE_COUNT = BASE[2]
			DIFF_COUNT = BASE[3]
	
			CONSENSUS_SEQ.append(CONSENSUS_BASE)
	    
			QUALITIES = get_qualities(CONSENSUS_BASE,SEQ[position],QUAL[position])
	    
			if (NUM_DIFF_BASES < 2 and CONSENSUS_BASE_COUNT > DIFF_COUNT):
		
				CONSENSUS_QUALITY_num = consensus_quality(QUALITIES,minBQ,DIFF_COUNT,STEP)
				CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
				CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)
		
			elif (NUM_DIFF_BASES >= 2 or CONSENSUS_BASE_COUNT <= DIFF_COUNT ):
			#else:
				CONSENSUS_QUALITY_num = 0
				CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
				CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)
	    
			else:
				print "Error"
				print SEQ[position],BASE
		
		else:
			CONSENSUS_BASE = most_common_base_low_qual(SEQ[position])[0]
			CONSENSUS_SEQ.append(CONSENSUS_BASE)

			CONSENSUS_QUALITY_num = 0
			CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
			CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)

	# We take the info from the last read in the group
	CONSENSUS_READ = i
	
	# Mapping quality == mean of reads' mapq
	CONSENSUS_READ.mapq = int(round(float(sum(MAPQ))/len(MAPQ)))
	
	# Consensus seq per position
	CONSENSUS_SEQ = ''.join(CONSENSUS_SEQ)
	CONSENSUS_READ.seq = CONSENSUS_SEQ
	
	# Base qualities are calcultated as the mean of the base qualities of each read.
	# In case there were more than one SNP in the position, its base quality was 0.
	CONSENSUS_QUAL = ''.join(CONSENSUS_QUAL)
	CONSENSUS_READ.qual = CONSENSUS_QUAL

    return (CONSENSUS_READ)


### Read parameters
parser = argparse.ArgumentParser(description='Correcting bamfiles using barcodes info')
parser.add_argument('--infile', required=True, dest='infile', help='Input BAM file.')
parser.add_argument('--outfile', required=True, dest='outfile', help='Output BAM file.')
parser.add_argument('--barcodes', required=False, dest='barcodes', type=str,choices=['BEGINNING', 'END', 'BOTH'], default='BOTH', help='Barcodes to use. BEGGINING = Barcode 1; END = Barcode 2; BOTH = Barcode 1 and 2.')
parser.add_argument('--minBQ', required=False, dest='minBQ', type=int, default=30, help='Minimum base quality to consider. Default = 30')
parser.add_argument('--BCerror', required=False, dest='BCerror', type=int, default=0, help='Maximum sequencing errors alowed in barcode sequence. Default = 0')
parser.add_argument('--step', required=False, dest='step', type=int, default=0, choices=[0, 1, 2], help='Genebits step. 0: Unique barcode correction; 1: Subfamily correction; 2: Family correction')


args = ''
try:
	args = parser.parse_args()
except IOError as io:
	print io
	sys.exit('Error reading parameters.')


### Input BAM
try:
	samfile = pysam.Samfile( args.infile, "rb" )
except:
	exit("Cannot open input file.")

### Output BAM
try:
	outfile = pysam.Samfile(args.outfile, mode="wb", template = samfile)
except:
	exit("Cannot open output file.")



logfile1 = args.outfile + ".log1"
LOGFILE1=open(logfile1,'w')

logfile2 = args.outfile + ".log2"
LOGFILE2=open(logfile2,'w')

minBQ = args.minBQ
Errors=float(args.BCerror)
STEP=args.step

if (args.barcodes == "BOTH"):
	Errors = Errors*2

pos=0
end=0
BARCODE = ""
READ = {}
ENDS = []
for read in samfile.fetch():
	length = str(read.rlen)
	ref_end = str(read.aend)
	ref_start = str(read.pos)
	
	ref_length = ref_end + "," + length
	
	# Getting the barcodes
	NAME = str(read.qname)
        BC = NAME.split(':')[-1]
        BC1 = BC.split(',')[0]
        BC2 = BC.split(',')[1]
	

	if (ref_start == pos):
		
		if ref_length in READ:
			READ[ref_length].append(read)
	   
		else:
			READ[ref_length] = list()
			READ[ref_length].append(read)

	
	else:
		#print len(READ)
		if (len(READ) > 0):
			#SAME_END_READS = []
			#print ENDS
			for i in READ: # Grouping reads with the same start and end

				log1 = [str(ref_start), str(i), str(len(READ[i]))]
				LOG1 = "\t".join(log1)+"\n"
				LOGFILE1.write(LOG1)
				
				DICTIONARY = extract_bc_groups(READ[i],args.barcodes)
				
				for barcode in DICTIONARY:
					#print barcode, str(pos), str(i), str(len(DICTIONARY[barcode]))

					log2 = [barcode, str(ref_start), str(i), str(len(DICTIONARY[barcode]))]
					LOG2 = "\t".join(log2)+"\n"
					LOGFILE2.write(LOG2)
				
					# Printing consensus reads to a new bam file
					NEW_READ = GET_FINAL_READ(DICTIONARY[barcode],minBQ,STEP)
					#NEW_READ = DICTIONARY[barcode][0]
					#print STEP, NEW_READ.qual
					outfile.write(NEW_READ)
			
			READ = {}
			#ENDS = []
			
			pos = ref_start
			end = ref_length
			
			READ[ref_length] = list()
			READ[ref_length].append(read)
			#ENDS.append(ref_end)
			
		else:
			READ = {}
			#ENDS = []
			
			pos = ref_start
			end = ref_end
			
			
			
			READ[ref_length] = list()
			READ[ref_length].append(read)



#### We need to print the last groups of reads
if (len(READ) > 0):
	#SAME_END_READS = []
	#print ENDS
	for i in READ: # Grouping reads with the same start and end

		log1 = [str(ref_start), str(i), str(len(READ[i]))]
		LOG1 = "\t".join(log1)+"\n"
		LOGFILE1.write(LOG1)
				
		DICTIONARY = extract_bc_groups(READ[i],args.barcodes)
				
		for barcode in DICTIONARY:
			#print barcode, str(pos), str(i), str(len(DICTIONARY[barcode]))

			log2 = [barcode, str(ref_start), str(i), str(len(DICTIONARY[barcode]))]
			LOG2 = "\t".join(log2)+"\n"
			LOGFILE2.write(LOG2)
				
			# Printing consensus reads to a new bam file
			NEW_READ = GET_FINAL_READ(DICTIONARY[barcode],minBQ,STEP)
			#NEW_READ = DICTIONARY[barcode][0]
			#print STEP, NEW_READ.qual
			outfile.write(NEW_READ)
		
		
		
samfile.close()
LOGFILE1.close()
LOGFILE2.close()
outfile.close()
			

