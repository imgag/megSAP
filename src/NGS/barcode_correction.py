#!/usr/bin/env python

import argparse
import pysam
import numpy
import timeit

start = timeit.default_timer()

## EXAMPLE COMMAND LINE
# python correct_bam_barcodes_21.11.2016.py --infile PATH/TO/test.bam --outfile PATH/TO/corrected.test.bam --barcodes BOTH

# Function to get the amount of differences
def similarity(a, b):
    return sum(x!=y for x,y in zip(a,b))

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
	##print bc
	return(bc)

def extract_bc_groups(BC, Errors): # Input, list of reads that start are end are equal
	BARCODE_DICT={}
	
	##print BC.keys()
	
	while len(BC) > 0:
		
		# The reads are stored in dict BC. Their keys are Barcode;Start_read_2.
                # We get the more common KEY and split it in Barcode and Start_read_2 (POS2).
                MAX = "".join(([k for k in BC.keys() if BC[k]==max(BC.values(),key=len)]))
                Barcode = MAX.split(';')[0]
                POS2 = MAX.split(';')[1]
		
                # We get those indexes that Barcodes are similar to the most common one (allowing Errors) and with the same Start_read_2 
		SIM = [BC.keys()[x] for x in range(0,len(BC.keys())) if (similarity(BC.keys()[x].split(';')[0], MAX) <= Errors and BC.keys()[x].split(';')[1] == POS2)]
		
		# We create a new key of the most common Barcode and we add later all the similar reads.
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
    QUALITIES_NATIVE = QUALITIES
    
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

# Change all base qualities of the read to 0 ("!")
def error_read_qual(read):
    READ = read
    qualities = READ.qual
    LEN = len(qualities)
    QUAL = ''.join(["!" for i in xrange(LEN)])
    READ.qual = QUAL	   
    return(READ)

def error_read_seq(read):
    READ = read
    seq = READ.seq
    LEN = len(seq)
    SEQ = ''.join(["N" for i in xrange(LEN)])
    READ.seq = SEQ	   
    return(READ)

def one_duplicate_qual(read, STEP, minBQ):
    #READ = reduce_mapq(read)
    READ = read
    
    # No Variant Quality score recalibration. Just error correction. Not necessary anything else in this 	
    if (STEP == 3):
        return (READ)
    
    qualities = READ.qual
    QUALITIES = ascii2dec(qualities)
    QUAL=[]
    QUALITIES_NATIVE = QUALITIES
    
    # Variant quality score recalibration for general barcode strategy
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

    # Variant quality score recalibration for sub-family barcode
    elif (STEP == 1):
        for i in range(0,len(QUALITIES_NATIVE)):
            if (QUALITIES_NATIVE[i] >= minBQ):
                QUAL.append(str(chr(10 + 33)))
            else:
                QUAL.append("!") # 0 quality = chr(0 + 33)
    
    # Variant quality score recalibration for family barcode
    elif (STEP == 2):
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
    
    READ.qual = QUAL	   
    return(READ)

def low_quality_base_check(QUALITIES,minBQ,STEP):
	QUALITIES_NATIVE = QUALITIES
	return (len(QUALITIES_NATIVE) >  len([x for x in range(0,len(QUALITIES_NATIVE)) if QUALITIES_NATIVE[x] < minBQ]))
    
def consensus_quality(QUALITIES,minBQ,ERRORS,STEP):
    QUALITIES_NATIVE = QUALITIES
    
    # How many bases with high good quality
    COPIES = len([x for x in range(0,len(QUALITIES_NATIVE)) if QUALITIES_NATIVE[x] >= minBQ])
        	
    # Consider more or equal than 5 as best agreement
    if (COPIES >= 5):
        NEW_COPIES = 5
    else:
        NEW_COPIES = COPIES    
    
    
    if (ERRORS > 5):
	ERRORS = 5
    
	# Double check
    MAX_QUAL = max(QUALITIES_NATIVE)
    if (MAX_QUAL < minBQ):
        NEW_QUAL = 0
    else:
	# Base qualities are correlated with how many reads with good quality are showing the same allele and how many of them had errors
	# COPIES * 10 + ERRORS
        if (STEP == 1):
		if ((ERRORS >= 1 and float(ERRORS)/(COPIES+ERRORS) > 0.25) or ERRORS >= 3):
		    NEW_QUAL = 0
		else:
		    NEW_QUAL = NEW_COPIES * 10 + int(ERRORS)
		    
        elif (STEP == 2 or STEP == 0):
		if ((ERRORS >= 1 and float(ERRORS)/(COPIES+ERRORS) > 0.25) or ERRORS >= 3):
			NEW_QUAL = 0
		else:	
			MEDIAN = int(round((numpy.mean(QUALITIES_NATIVE))/10,0))
			
			# Fixing max value
			if (MEDIAN > 5):
				MEDIAN = 5
				
			NEW_QUAL = NEW_COPIES * 10 + MEDIAN
		
	elif (STEP == 3):
		if ((ERRORS >= 1 and float(ERRORS)/(COPIES+ERRORS) > 0.25) or ERRORS >= 3):
			NEW_QUAL = 0
		else:
		    # We take the max base quality as the concensus one
		    MAX = max(QUALITIES_NATIVE)
		    NEW_QUAL = MAX
    
    return (NEW_QUAL)
    
# Function to correct reads grouped by barcodes
def GET_FINAL_READ(reads,minBQ,STEP,SET_N):
    
    CONSENSUS_SEQ = list()
    CONSENSUS_QUAL = list()
    
    # Getting the amount 
    #print len(reads)    
    
    # Objects to save info for lexicoorder of the reads
    LIST_READNAMES = list()
    DICT_READS = {}
    READNAME = ''
    
    if (len(reads) <= 1):
		for i in reads:
		    
		    # Encoding quality
		    CONSENSUS_READ = one_duplicate_qual(i, STEP, minBQ)
		
		
		    # We add info about the amount of duplicates per family group
		    count = len(reads)
		    color = ['230,242,255','179,215,255','128,187,255','77,160,255','26,133,255']
		    current_color = '0,0,0'
		    if (count > 5):
			current_color = color[4]
		    else:
			current_color = color[count-1]
		    CONSENSUS_READ.tags += [('DP', count)]
		    CONSENSUS_READ.tags += [('YC', current_color)]

		    # Info about barcode groups
		    LOG_INFO = (CONSENSUS_READ.qname, str(CONSENSUS_READ.pos), str(len(reads))) 
		    LOG_INFO = "\t".join(LOG_INFO)+"\n"
		    
    else:
		#SEQ = collections.OrderedDict()
		#QUAL = collections.OrderedDict()
		SEQ = {}
		QUAL = {}
		MAPQ = list()
		
		# Variable count for compare indels between duplicates
		# The first read is taken as reference
		##print reads[0].qname, reads[0].cigarstring
		if (reads[0].cigarstring == None):
		    Ind_count = 0
		else:
		    Ind_count = reads[0].cigarstring.count("I") + reads[0].cigarstring.count("D")
		
		for i in reads:
						
			#print i.pos, i.aend, i.rlen, len(i.seq), i.qname, i.cigarstring, i.qstart, i.query
			#print i.seq, i.qname, i.qual, i.cigarstring

			# In case that the amount of indels differ  between duplicates, we take the the first read in the lexico order and we label all them base qualities to 0
			if (i.cigarstring == None or (i.cigarstring.count("I") + i.cigarstring.count("D")) != Ind_count):
				
			    for j in reads:
				# Adding reads to a hash to later, sort them in lexicographical order
				READNAME = j.qname
			        DICT_READS[READNAME] = j
			        LIST_READNAMES.append(READNAME)
			    
			    # We take the info from the last read in the group
			    SORTED_READNAMES = sorted(LIST_READNAMES)
			    
			    # When problems with Different indels between duplicates appear, we get as consensus read the first lexico read, and change or their base qualities to 0
			    if (not SET_N):
				CONSENSUS_READ = error_read_qual(DICT_READS[SORTED_READNAMES[0]])
			    else:
				CONSENSUS_READ = error_read_seq(DICT_READS[SORTED_READNAMES[0]])
				
			    ###print "ERROR2", CONSENSUS_READ.qname, len(CONSENSUS_SEQ), DICT_READS[SORTED_READNAMES[0]].cigarstring, 
			    
			    # We add info about the amount of duplicates per family group
			    count = len(reads)
			    color = ['230,242,255','179,215,255','128,187,255','77,160,255','26,133,255']
			    current_color = '0,0,0'
			    if (count > 5):
				current_color = color[4]
			    else:
				current_color = color[count-1]
			    CONSENSUS_READ.tags += [('DP', count)]
			    CONSENSUS_READ.tags += [('YC', current_color)]
			    
			    LOG_INFO = (CONSENSUS_READ.qname, str(CONSENSUS_READ.pos), str(len(reads))) 
			    LOG_INFO = "\t".join(LOG_INFO)+"\n"
			    return (CONSENSUS_READ, LOG_INFO)
			
			
			# Adding reads to a hash to later, sort them in lexicographical order
			READNAME = i.qname
			DICT_READS[READNAME] = i
			LIST_READNAMES.append(READNAME)
			
			LEN = i.rlen
			##print LEN
			seq = i.seq
			qual = i.qual
			MAPQ.append(i.mapq)
			
			
			SOFT = i.pos - i.qstart
			##print "SOFT", i.pos, i.qstart
			
			for b in range(0,LEN):
			    BASE = seq[b]
			    BASE_QUAL = qual[b]
			    
			    B = SOFT + b
			    
			    if B in SEQ:
				SEQ[B].append(BASE)
				QUAL[B].append(BASE_QUAL)
				    
			    else:
				SEQ[B] = list()
				QUAL[B] = list()
				SEQ[B].append(BASE)
				QUAL[B].append(BASE_QUAL)		

		for position in sorted(SEQ):
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
					
				QUALITIES = get_qualities(CONSENSUS_BASE,SEQ[position],QUAL[position])
			
				if (NUM_DIFF_BASES < 3 and CONSENSUS_BASE_COUNT > DIFF_COUNT):
					CONSENSUS_QUALITY_num = consensus_quality(QUALITIES,minBQ,DIFF_COUNT,STEP)

					if (SET_N and CONSENSUS_QUALITY_num==0):
						CONSENSUS_QUALITY_num = QUALITIES[0]
						CONSENSUS_BASE = "N"
										
					CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
					CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)
					CONSENSUS_SEQ.append(CONSENSUS_BASE)
			
				elif (NUM_DIFF_BASES >= 3 or CONSENSUS_BASE_COUNT <= DIFF_COUNT ):
					CONSENSUS_QUALITY_num = 0

					if (SET_N):
						CONSENSUS_QUALITY_num = QUALITIES[0]
						CONSENSUS_BASE = "N"
						
					CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
					CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)
					CONSENSUS_SEQ.append(CONSENSUS_BASE)
			
				else:
					print "Error"
					#print SEQ[position], BASE
			
			else:
				CONSENSUS_BASE = most_common_base_low_qual(SEQ[position])[0]
				CONSENSUS_QUALITY_num = 0

				if (SET_N):
					QUALITIES = get_qualities(CONSENSUS_BASE,SEQ[position],QUAL[position])
					CONSENSUS_QUALITY_num = QUALITIES[0]
					CONSENSUS_BASE = "N"
				
				CONSENSUS_SEQ.append(CONSENSUS_BASE)
				CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
				CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)
		    
		# We take the info from the last read in the group
		SORTED_READNAMES = sorted(LIST_READNAMES)
		
		##print LIST_READNAMES
		READ_COUNT = 0

		while (READ_COUNT < len(SORTED_READNAMES) and len(CONSENSUS_SEQ) != DICT_READS[SORTED_READNAMES[READ_COUNT]].rlen):
		    ##print "NO", len(CONSENSUS_SEQ), DICT_READS[SORTED_READNAMES[READ_COUNT]].rlen, DICT_READS[SORTED_READNAMES[0]].cigarstring, DICT_READS[SORTED_READNAMES[READ_COUNT]].rlen, DICT_READS[SORTED_READNAMES[READ_COUNT]].cigarstring
		    READ_COUNT = READ_COUNT + 1
		
		if (READ_COUNT < len(SORTED_READNAMES) and len(CONSENSUS_SEQ) == DICT_READS[SORTED_READNAMES[READ_COUNT]].rlen):
		    CONSENSUS_READ = DICT_READS[SORTED_READNAMES[READ_COUNT]]
		    CONSENSUS_READ.qname = DICT_READS[SORTED_READNAMES[0]].qname
		    ##print ''.join(CONSENSUS_SEQ), ''.join(CONSENSUS_QUAL)
		    ##print "YES", CONSENSUS_READ.qname, len(CONSENSUS_SEQ), DICT_READS[SORTED_READNAMES[0]].cigarstring, DICT_READS[SORTED_READNAMES[READ_COUNT]].rlen, DICT_READS[SORTED_READNAMES[READ_COUNT]].cigarstring

		else:
		    # When problems with Different indels between duplicates appear, we get as consensus read the first lexico read, and change or their base qualities to 0
		    if (not SET_N):
			CONSENSUS_READ = error_read_qual(DICT_READS[SORTED_READNAMES[0]])
		    else:
			CONSENSUS_READ = error_read_seq(DICT_READS[SORTED_READNAMES[0]])
		    
		    ##print "ERROR", CONSENSUS_READ.qname, len(CONSENSUS_SEQ), DICT_READS[SORTED_READNAMES[0]].cigarstring, DICT_READS[SORTED_READNAMES[READ_COUNT]].rlen, DICT_READS[SORTED_READNAMES[READ_COUNT]].cigarstring
		    
		    # We add info about the amount of duplicates per family group
		    count = len(reads)
		    color = ['230,242,255','179,215,255','128,187,255','77,160,255','26,133,255']
		    current_color = '0,0,0'
		    if (count > 5):
			current_color = color[4]
		    else:
			current_color = color[count-1]
		    CONSENSUS_READ.tags += [('DP', count)]
		    CONSENSUS_READ.tags += [('YC', current_color)]
		
		    LOG_INFO = (CONSENSUS_READ.qname, str(CONSENSUS_READ.pos), str(len(reads))) 
		    LOG_INFO = "\t".join(LOG_INFO)+"\n"
		    return (CONSENSUS_READ, LOG_INFO)

		
		# Mapping quality == mean of reads' mapq
		CONSENSUS_READ.mapq = int(round(float(sum(MAPQ))/len(MAPQ)))
	
		# Consensus seq per position
		CONSENSUS_SEQ = ''.join(CONSENSUS_SEQ)
		CONSENSUS_READ.seq = CONSENSUS_SEQ
		
		# Base qualities are calcultated as the mean of the base qualities of each read.
		# In case there were more than one SNP in the position, its base quality was 0.
		CONSENSUS_QUAL = ''.join(CONSENSUS_QUAL)
		CONSENSUS_READ.qual = CONSENSUS_QUAL
		
		# We add info about the amount of duplicates per family group
		count = len(reads)
		color = ['230,242,255','179,215,255','128,187,255','77,160,255','26,133,255']
		current_color = '0,0,0'
		if (count > 5):
			current_color = color[4]
		else:
			current_color = color[count-1]
		CONSENSUS_READ.tags += [('DP', count)]
		CONSENSUS_READ.tags += [('YC', current_color)]
		    		
		# Info about barcode groups
		LOG_INFO = (CONSENSUS_READ.qname, str(CONSENSUS_READ.pos), str(len(reads))) 
		LOG_INFO = "\t".join(LOG_INFO)+"\n"
		#print "\n"
		
    return (CONSENSUS_READ, LOG_INFO)


### Read parameters
parser = argparse.ArgumentParser(description='Correcting bamfiles using barcodes info')
parser.add_argument('--infile', required=True, dest='infile', help='Input BAM file.')
parser.add_argument('--outfile', required=True, dest='outfile', help='Output BAM file.')
parser.add_argument('--barcodes', required=False, dest='barcodes', type=str,choices=['BEGINNING', 'END', 'BOTH'], default='BOTH', help='Barcodes to use. BEGGINING = Barcode 1; END = Barcode 2; BOTH = Barcode 1 and 2. Default = BOTH')
parser.add_argument('--minBQ', required=False, dest='minBQ', type=int, default=30, help='Minimum base quality to be considered. Default = 30')
parser.add_argument('--BCerror', required=False, dest='BCerror', type=int, default=0, help='Maximum number of sequencing errors allowed in barcode sequence. Default = 0')
parser.add_argument('--step', required=False, dest='step', type=int, default=0, choices=[0, 1, 2, 3], help='Protocol step. 0: Unique barcode correction; 1: Subfamily correction; 2: Family correction; 3: Not full Base quality recalibration, just error correction. Default = 0')
parser.add_argument('--n', required=False, dest='n', action='store_true', help='Use Ns instead of reducing base quality.')


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

logfile2 = args.outfile + ".log"
LOGFILE2=open(logfile2,'w')

minBQ = args.minBQ
Errors=(args.BCerror)
STEP=args.step
SET_N = args.n


pos=0
end=0
BARCODE = ""
READ = {}
POSITIONS_DICT = {}
ENDS = []
for read in samfile.fetch():
	length = str(read.rlen)
	ref_end = str(read.aend)
	ref_start = str(read.pos)
	
	# Ref start of the paired read
	#ref_length = str(read.next_reference_start)
	ref_length = str(read.tlen)
	
	# Getting the barcodes
	NAME = str(read.qname)
	bc = extract_barcode(read,args.barcodes) # Extract the barcode
	
	CODE = bc + ";" + ref_length 

	if (ref_start == pos):
		
		# To store the CODEs for each ref_length
		if ref_length in POSITIONS_DICT:
			POSITIONS_DICT[ref_length].append(CODE)
		else:
			POSITIONS_DICT[ref_length] = list()
			POSITIONS_DICT[ref_length].append(CODE)
			
		# To create a Dictionary to store all reads with the code as key (create above)
		if CODE in READ:
			READ[CODE].append(read)
		else:
			READ[CODE] = list()
			READ[CODE].append(read)

	
	else:

		if (len(READ) > 0 and Errors > 0):
		    ##print POSITIONS_DICT.keys(), len(READ), len(POSITIONS_DICT)
		    for pos2 in POSITIONS_DICT:
			##print pos2, len(POSITIONS_DICT[pos2])

			NEW_READS =  dict((key,value) for key, value in READ.iteritems() if key in POSITIONS_DICT[pos2])
			
			# When we allow errors in the Barcodes, we re-group them by similarity (Errors specified in parameter)
			DICTIONARY = extract_bc_groups(NEW_READS,Errors)
			#print len(DICTIONARY)
			
			for barcode in DICTIONARY:
			    # Printing consensus reads to a new bam file
			    NEW_READ,LOG2 = GET_FINAL_READ(DICTIONARY[barcode],minBQ,STEP,SET_N)
			    #NEW_READ = DICTIONARY[barcode][0]
			    ##print STEP, NEW_READ.qual
			    
			    LOGFILE2.write(LOG2)
			    outfile.write(NEW_READ)
			
		    READ = {}
		    POSITIONS_DICT = {}
			
		    pos = ref_start
		    end = ref_length
			
		    READ[CODE] = list()
		    READ[CODE].append(read)
			
		    POSITIONS_DICT[ref_length] = list()
		    POSITIONS_DICT[ref_length].append(CODE)
			
		elif (len(READ) > 0 and Errors == 0):
		    
		    DICTIONARY = READ
		    ##print len(DICTIONARY)
		    for barcode in DICTIONARY:
			# printing consensus reads to a new bam file
			NEW_READ,LOG2 = GET_FINAL_READ(DICTIONARY[barcode],minBQ,STEP,SET_N)
			#NEW_READ = DICTIONARY[barcode][0]
			##print STEP, NEW_READ.qual
			
			LOGFILE2.write(LOG2)
			outfile.write(NEW_READ)
		
		    READ = {}
		    POSITIONS_DICT = {}
		    
		    pos = ref_start
		    end = ref_length
		    
		    READ[CODE] = list()
		    READ[CODE].append(read)
		    
		    POSITIONS_DICT[ref_length] = list()
		    POSITIONS_DICT[ref_length].append(CODE)
		    	
		else:
			READ = {}
			POSITIONS_DICT = {}
			
			pos = ref_start
			end = ref_end
						
			READ[CODE] = list()
			READ[CODE].append(read)
			
			POSITIONS_DICT[ref_length] = list()
			POSITIONS_DICT[ref_length].append(CODE)




#### We need to print the last groups of reads
if (len(READ) > 0 and Errors > 0):
    for pos2 in POSITIONS_DICT:
	NEW_READS =  dict((key,value) for key, value in READ.iteritems() if key in POSITIONS_DICT[pos2])
			
	# When we allow errors in the Barcodes, we re-group them by similarity (Errors specified in parameter)
	DICTIONARY = extract_bc_groups(NEW_READS,Errors)
	
	for barcode in DICTIONARY:
	    
	    # Printing consensus reads to a new bam file
	    NEW_READ,LOG2 = GET_FINAL_READ(DICTIONARY[barcode],minBQ,STEP,SET_N)
	    #NEW_READ = DICTIONARY[barcode][0]
	    #print STEP, NEW_READ.qual
	    
	    LOGFILE2.write(LOG2)
	    outfile.write(NEW_READ)

elif (len(READ) > 0 and Errors == 0):
	DICTIONARY = READ

	for barcode in DICTIONARY:
	    # Printing consensus reads to a new bam file
	    NEW_READ,LOG2 = GET_FINAL_READ(DICTIONARY[barcode],minBQ,STEP,SET_N)
	
	    LOGFILE2.write(LOG2)
	    outfile.write(NEW_READ)
	
		
samfile.close()
LOGFILE2.close()
outfile.close()
			
stop = timeit.default_timer()
print 'TIME'
print stop - start
