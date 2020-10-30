#!/usr/bin/env python

import argparse
import pysam
import numpy
import timeit
import re
import networkx as nx

start = timeit.default_timer()


## EXAMPLE COMMAND LINE
# python correct_bam_barcodes_21.11.2016.py --infile PATH/TO/test.bam --outfile PATH/TO/corrected.test.bam --barcodes BOTH

# Function to get the amount of differences
def similarity(a, b):
    return sum(x != y for x, y in zip(a, b))


def keywithmaxval(d):
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


def update_tag(TAG, VALUE):
    return [VALUE if x[0] == VALUE[0] else x for x in TAG]


def get_edge(CODE, LIST, Errors):
    EDGE = [[CODE, LIST[x]] for x in range(0, len(LIST)) if (similarity(LIST[x], CODE) <= Errors)]
    return EDGE


def extract_barcode(READ, BC_TYPE):  # Individual read. It returns the barcode of the read
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
    return (bc)


def extract_bc_groups(BC, BC_NETWORK):  # Input, list of reads that start are end are equal
    BARCODE_DICT = {}

    sortedKeyList = sorted(BC.keys(), key=lambda s: len(BC.get(s)), reverse=True)

    while len(BC) > 0:

        # The reads are stored in dict BC. Their keys are the Barcodes.
        # We get the most frequent KEY (Barcode)
        MAX = sortedKeyList[0]

        # We create a new key of the most common Barcode and we add later all the similar reads.
        BARCODE_DICT[
            MAX] = list()  # Creating key in the hash based on our barcode where reads of this group will be saved

        SIM = list(BC_NETWORK.adj[MAX])

        for i in SIM:
            BARCODE_DICT[MAX].extend(BC[i])  # Grouping reads based on similarity of the barcodes
            del BC[i]  # Removing barcodes already considered
            BC_NETWORK.remove_node(i)
            sortedKeyList.remove(i)

    return BARCODE_DICT  # Dictionary with Barcode as a key and reads as values


def ReverseComplement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
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
    qual = [ord(i) - 33 for i in ASCII]
    return qual


def most_common_base(LIST, QUALITIES, minBQ):
    QUALITIES_NATIVE = QUALITIES

    HQ = [x for x in range(0, len(QUALITIES_NATIVE)) if QUALITIES_NATIVE[x] >= minBQ]

    result_list = [LIST[i] for i in HQ]
    result_qual = [QUALITIES[i] for i in HQ]

    BASE = max(sorted(set(result_list)), key=result_list.count)
    NUM_DIFF_BASES = len(set(result_list))

    DIFF_COUNT = len([x for x in range(0, len(result_list)) if result_list[x] != BASE])
    BASE_COUNT = len([x for x in range(0, len(result_list)) if result_list[x] == BASE])

    return BASE, NUM_DIFF_BASES, BASE_COUNT, DIFF_COUNT


def most_common_base_low_qual(LIST):
    BASE = max(sorted(set(LIST)), key=LIST.count)
    NUM_DIFF_BASES = len(set(LIST))

    return BASE, NUM_DIFF_BASES


def recalculate_NM(READ):
    #   # There is a problem when the cigar length does not fit to the MD tag
    #
    #   print READ.reference_name, READ.pos
    #   print READ.cigarstring, get_md_reference_length(read.get_tag("MD"))

    #
    #   refSeq = READ.get_reference_sequence()
    #   readSeq = READ.query_alignment_sequence
    #
    #
    #   if ('I' not in READ.cigarstring and 'I' not in READ.cigarstring and len(refSeq) == len(readSeq)):
    #       New_NM = sum(x.upper()!=y.upper() for x,y in zip(refSeq,readSeq))
    #   else:
    #   New_NM = READ.opt("NM")
    New_NM = READ.opt("NM")
    return New_NM


def get_qualities(BASE, BASES, QUALITIES):
    QUAL = [(ord(QUALITIES[x]) - 33) for x in range(0, len(BASES)) if BASES[x] == BASE]
    # LIST = [BASES[x] for x in range(0,len(BASES)) if BASES[x] == BASE]
    return QUAL


# Change all base qualities of the read to 0 ("!")
def error_read_qual(read):
    READ = read
    qualities = READ.qual
    LEN = len(qualities)
    QUAL = ''.join(["!" for i in xrange(LEN)])
    READ.qual = QUAL
    return (READ)


# Change base errors to Ns
def error_read_seq(read):
    READ = read
    seq = READ.seq
    LEN = len(seq)
    SEQ = ''.join(["N" for i in xrange(LEN)])
    READ.seq = SEQ
    return (READ)


# Check if there is at least one high quality base
def low_quality_base_check(QUALITIES, minBQ):
    return (max(QUALITIES) >= minBQ)


# Get consensus read qualities
def consensus_quality(QUALITIES, minBQ, ERRORS, STEP):
    QUALITIES_NATIVE = QUALITIES

    # How many bases with high good quality
    COPIES = len([x for x in range(0, len(QUALITIES_NATIVE)) if QUALITIES_NATIVE[x] >= minBQ])

    MAX_QUAL = max(QUALITIES_NATIVE)
    if (MAX_QUAL < minBQ):
        NEW_QUAL = MAX_QUAL
    else:
        # We label errors with base quality 0
        # For consensus bases, we take the max as the new quality base
        if (STEP == 1 or STEP == 2):
            if ((ERRORS >= 1 and float(ERRORS) / (COPIES + ERRORS) > 0.25) or ERRORS >= 3):
                NEW_QUAL = 0
            else:
                # We take the max base quality as the concensus one
                MAX1 = max(QUALITIES_NATIVE)
                MAX = max(30, MAX1)
                NEW_QUAL = MAX

    return (NEW_QUAL)


# Function to correct reads grouped by barcodes
def GET_FINAL_READ(reads, minBQ, STEP, SET_N):
    CONSENSUS_SEQ = list()
    CONSENSUS_QUAL = list()

    # Getting the amount
    # print(len(reads))

    # Objects to save info for lexicoorder of the reads
    LIST_READNAMES = list()
    DICT_READS = {}
    READNAME = ''
    FLAG = ''
    DICT_FLAG = {}
    DP1_TAG = list()

    if (len(reads) <= 1):
        for i in reads:

            # Encoding quality
            CONSENSUS_READ = i

            # We add info about the amount of duplicates per family group
            count = len(reads)
            color = ['230,242,255', '179,215,255', '128,187,255', '77,160,255', '26,133,255']
            current_color = '0,0,0'
            if (count > 5):
                current_color = color[4]
            else:
                current_color = color[count - 1]

            # adding barcodes to tag in bam file
            if (STEP == 1):
                CONSENSUS_READ.tags += [('DP', count)]
                CONSENSUS_READ.tags += [('YC', current_color)]
            elif (STEP == 2):
                CONSENSUS_READ.tags += [('DF', count)]

            # Info about barcode groups
            LOG_INFO = (CONSENSUS_READ.qname, str(CONSENSUS_READ.pos), str(len(reads)))
            LOG_INFO = "\t".join(LOG_INFO) + "\n"

    else:
        REF_LENGTH = [(k.reference_length, k.rlen, k.cigarstring) for k in reads]

        MAX_INFO = max(sorted(set(REF_LENGTH), reverse=True), key=REF_LENGTH.count)

        # Reference length
        MAX_REF_LENGTH = MAX_INFO[0]

        # Read length
        MAX_READ_LENGTH = MAX_INFO[1]

        # Max cigarstring
        MAX_CIGAR_LENGTH = MAX_INFO[2]

        SEQ = {}
        QUAL = {}
        MAPQ = list()
        LQ_read = 0

        # Variable count for compare indels between duplicates
        # The first read is taken as reference

        # if (reads[0].cigarstring == None):
        #    Ind_count = 0
        # else:
        #    Ind_count = reads[0].cigarstring.count("I") + reads[0].cigarstring.count("D")

        for i in reads:

            # In case that the amount of indels differ  between duplicates, we take the the first read in the lexico order and we label all them base qualities to 0
            if ((MAX_REF_LENGTH != i.reference_length) or (MAX_READ_LENGTH != i.rlen) or (MAX_CIGAR_LENGTH != i.cigarstring)):  # and (i.cigarstring == None or "I" in i.cigarstring or "D" in i.cigarstring)):
                # if ((MAX_REF_LENGTH != i.reference_length or MAX_READ_LENGTH != i.rlen) and (i.cigarstring == None or "I" in i.cigarstring or "D" in i.cigarstring)):
                # print 'BULSHIT READ'
                # print MAX_REF_LENGTH, len(i.get_reference_sequence()), i.rlen
                READNAME = i.qname
                FLAG = i.flag
                DICT_READS[READNAME] = i
                DICT_FLAG[READNAME] = FLAG
                LIST_READNAMES.append(READNAME)
                LQ_read = LQ_read + 1
                continue

            # else:
            #    print 'PASS READ'
            #    print MAX_REF_LENGTH, len(i.get_reference_sequence()), i.rlen

            # Saving the amount of duplicates from first round of correction
            LAST_READ = i

            if (STEP == 2):
                DP1_TAG.append(i.opt("DP"))

            # Adding reads to a hash to later, sort them in lexicographical order
            READNAME = i.qname
            FLAG = i.flag
            DICT_READS[READNAME] = i
            DICT_FLAG[READNAME] = FLAG
            LIST_READNAMES.append(READNAME)

            LEN = i.rlen
            ##print LEN
            seq = i.seq
            qual = i.qual
            MAPQ.append(i.mapq)

            SOFT = i.pos - i.qstart

            for b in range(0, LEN):
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
            if (low_quality_base_check(Q, minBQ)):
                BASE = most_common_base(SEQ[position], Q, minBQ)

                CONSENSUS_BASE = BASE[0]
                NUM_DIFF_BASES = BASE[1]
                CONSENSUS_BASE_COUNT = BASE[2]
                DIFF_COUNT = BASE[3]

                QUALITIES = get_qualities(CONSENSUS_BASE, SEQ[position], QUAL[position])

                if (NUM_DIFF_BASES < 3 and CONSENSUS_BASE_COUNT > DIFF_COUNT):
                    CONSENSUS_QUALITY_num = consensus_quality(QUALITIES, minBQ, DIFF_COUNT, STEP)

                    if (SET_N and CONSENSUS_QUALITY_num == 0):
                        CONSENSUS_QUALITY_num = QUALITIES[0]
                        CONSENSUS_BASE = "N"

                    CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
                    CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)
                    CONSENSUS_SEQ.append(CONSENSUS_BASE)

                elif (NUM_DIFF_BASES >= 3 or CONSENSUS_BASE_COUNT <= DIFF_COUNT):
                    CONSENSUS_QUALITY_num = 0

                    if (SET_N):
                        CONSENSUS_QUALITY_num = QUALITIES[0]
                        CONSENSUS_BASE = "N"

                    CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
                    CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)
                    CONSENSUS_SEQ.append(CONSENSUS_BASE)

                else:
                    print("Error")
                    # print SEQ[position], BASE

            else:
                CONSENSUS_BASE = most_common_base_low_qual(SEQ[position])[0]
                CONSENSUS_QUALITY_num = 0

                if (SET_N):
                    QUALITIES = get_qualities(CONSENSUS_BASE, SEQ[position], QUAL[position])
                    CONSENSUS_QUALITY_num = QUALITIES[0]
                    CONSENSUS_BASE = "N"

                CONSENSUS_SEQ.append(CONSENSUS_BASE)
                CONSENSUS_QUALITY_ascii = chr(CONSENSUS_QUALITY_num + 33)
                CONSENSUS_QUAL.append(CONSENSUS_QUALITY_ascii)

        # We take the info from the last read in the group
        SORTED_READNAMES = sorted(LIST_READNAMES)

        ##print LIST_READNAMES
        READ_COUNT = 0

        # We take as template the last HQ read, but we change the read name and the flag
        CONSENSUS_READ = LAST_READ
        CONSENSUS_READ.qname = SORTED_READNAMES[0]

        # Mapping quality == mean of reads' mapq
        CONSENSUS_READ.mapq = int(round(float(sum(MAPQ)) / len(MAPQ)))
        CONSENSUS_READ.flag = DICT_FLAG[SORTED_READNAMES[0]]

        # Consensus seq per position
        CONSENSUS_SEQ = ''.join(CONSENSUS_SEQ)
        CONSENSUS_READ.seq = CONSENSUS_SEQ

        # Base qualities are calcultated as the mean of the base qualities of each read.
        # In case there were more than one SNP in the position, its base quality was 0.
        CONSENSUS_QUAL = ''.join(CONSENSUS_QUAL)
        CONSENSUS_READ.qual = CONSENSUS_QUAL

        # We add info about the amount of duplicates per family group
        count = len(reads) - LQ_read
        color = ['230,242,255', '179,215,255', '128,187,255', '77,160,255', '26,133,255']
        current_color = '0,0,0'

        if (count > 5):
            current_color = color[4]
        else:
            current_color = color[count - 1]

        # Different correction type
        if (STEP == 1):
            # New_NM = recalculate_NM(CONSENSUS_READ)

            # Add DP tag
            CONSENSUS_READ.tags += [('DP', count)]

            # Add color
            CONSENSUS_READ.tags += [('YC', current_color)]

            ## Update NM tag
            # CONSENSUS_READ.tags = update_tag(CONSENSUS_READ.tags,('NM', New_NM))

        elif (STEP == 2):
            DP1_value = int(round(numpy.median(DP1_TAG)))
            # New_NM = recalculate_NM(CONSENSUS_READ)

            # Update DP tag
            CONSENSUS_READ.tags = update_tag(CONSENSUS_READ.tags, ('DP', DP1_value))

            # Add DF tag
            CONSENSUS_READ.tags += [('DF', count)]

            # Update Color tag
            CONSENSUS_READ.tags = update_tag(CONSENSUS_READ.tags, ('YC', current_color))

            ## Update NM tag
            # CONSENSUS_READ.tags = update_tag(CONSENSUS_READ.tags,('NM', New_NM))

        # Info about barcode groups
        LOG_INFO = (CONSENSUS_READ.qname, str(CONSENSUS_READ.pos), str(len(reads)))
        LOG_INFO = "\t".join(LOG_INFO) + "\n"

    return (CONSENSUS_READ, LOG_INFO)


### Read parameters
parser = argparse.ArgumentParser(description='Correcting bamfiles using barcodes info')
parser.add_argument('--infile', required=True, dest='infile', help='Input BAM file.')
parser.add_argument('--outfile', required=True, dest='outfile', help='Output BAM file.')
parser.add_argument('--barcodes', required=False, dest='barcodes', type=str, choices=['BEGINNING', 'END', 'BOTH'],
                    default='BOTH',
                    help='Barcodes to use. BEGGINING = Barcode 1; END = Barcode 2; BOTH = Barcode 1 and 2. Default = BOTH')
parser.add_argument('--minBQ', required=False, dest='minBQ', type=int, default=10,
                    help='Minimum base quality to be considered. Default = 30')
parser.add_argument('--BCerror', required=False, dest='BCerror', type=int, default=0,
                    help='Maximum number of sequencing errors allowed in barcode sequence. Default = 0')
parser.add_argument('--step', required=False, dest='step', type=int, default=1, choices=[1, 2],
                    help='Protocol step. 1: Unique barcode correction; 2: Family correction. Default = 1')
parser.add_argument('--n', required=False, dest='n', action='store_true',
                    help='Use Ns instead of reducing base quality.')

args = ''
try:
    args = parser.parse_args()
except IOError as io:
    print(io)
    sys.exit('Error reading parameters.')

### Input BAM
try:
    samfile = pysam.Samfile(args.infile, "rb")
except:
    exit("Cannot open input file.")

### Output BAM
try:
    outfile = pysam.Samfile(args.outfile, mode="wb", template=samfile)
except:
    exit("Cannot open output file.")

logfile2 = args.outfile + ".log"
LOGFILE2 = open(logfile2, 'w')

minBQ = args.minBQ
Errors = (args.BCerror)
STEP = args.step
SET_N = args.n

pos = 0
end = 0
BARCODE = ""
POSITIONS_DICT = {}
ENDS = []
UNIQUE_BARCODES = {}
EDGE = ''
for read in samfile.fetch():
    if not read.is_secondary:
        length = str(read.rlen)
        ref_end = str(read.aend)
        ref_start = str(read.pos)

        # Both are required. Start of next read, and tlen shows the sign of othe read (- or +), which helps to separate pair reads when they map to the same coordinates
        ref_length = str(read.next_reference_start) + ',' + str(read.tlen)

        # Getting the barcodes
        NAME = str(read.qname)
        bc = extract_barcode(read, args.barcodes)  # Extract the barcode

        # CODE = bc + ";" + ref_length
        CODE = bc

        if (ref_start == pos):

            # To store the CODEs for each ref_length
            if ref_length in POSITIONS_DICT:
                if CODE in POSITIONS_DICT[ref_length]:
                    POSITIONS_DICT[ref_length][CODE].append(read)
                else:
                    POSITIONS_DICT[ref_length][CODE] = list()
                    POSITIONS_DICT[ref_length][CODE].append(read)
            else:
                POSITIONS_DICT[ref_length] = {}
                POSITIONS_DICT[ref_length][CODE] = list()
                POSITIONS_DICT[ref_length][CODE].append(read)

            # Allowing errors
            if (Errors > 0):
                if ref_length in UNIQUE_BARCODES:
                    if CODE in list(UNIQUE_BARCODES[ref_length].nodes()):
                        UNIQUE_BARCODES[ref_length].add_node(CODE)
                    else:
                        UNIQUE_BARCODES[ref_length].add_node(CODE)
                        EDGE = get_edge(CODE, list(UNIQUE_BARCODES[ref_length].nodes()), Errors)
                        UNIQUE_BARCODES[ref_length].add_edges_from(EDGE)
                else:
                    UNIQUE_BARCODES[ref_length] = nx.Graph()
                    UNIQUE_BARCODES[ref_length].add_node(CODE)
                    EDGE = get_edge(CODE, list(UNIQUE_BARCODES[ref_length].nodes()), Errors)
                    UNIQUE_BARCODES[ref_length].add_edges_from(EDGE)

        else:

            if (len(POSITIONS_DICT) > 0 and Errors > 0):
                for pos2 in POSITIONS_DICT:
                    # When we allow errors in the Barcodes, we re-group them by similarity (Errors specified in parameter)
                    DICTIONARY = extract_bc_groups(POSITIONS_DICT[pos2], UNIQUE_BARCODES[pos2])

                    for barcode in DICTIONARY:
                        # Printing consensus reads to a new bam file
                        NEW_READ, LOG2 = GET_FINAL_READ(list(DICTIONARY[barcode]), minBQ, STEP, SET_N)
                        LOGFILE2.write(LOG2)
                        outfile.write(NEW_READ)

                POSITIONS_DICT = {}
                UNIQUE_BARCODES = {}
                EDGE = ''

                pos = ref_start
                end = ref_length

                if ref_length in POSITIONS_DICT:
                    if CODE in POSITIONS_DICT[ref_length]:
                        POSITIONS_DICT[ref_length][CODE].append(read)
                    else:
                        POSITIONS_DICT[ref_length][CODE] = list()
                        POSITIONS_DICT[ref_length][CODE].append(read)
                else:
                    POSITIONS_DICT[ref_length] = {}
                    POSITIONS_DICT[ref_length][CODE] = list()
                    POSITIONS_DICT[ref_length][CODE].append(read)

                # Allowing errors
                if (Errors > 0):
                    if ref_length in UNIQUE_BARCODES:
                        if CODE in list(UNIQUE_BARCODES[ref_length].nodes()):
                            UNIQUE_BARCODES[ref_length].add_node(CODE)
                        else:
                            UNIQUE_BARCODES[ref_length].add_node(CODE)
                            EDGE = get_edge(CODE, list(UNIQUE_BARCODES[ref_length].nodes()), Errors)
                            UNIQUE_BARCODES[ref_length].add_edges_from(EDGE)
                    else:
                        UNIQUE_BARCODES[ref_length] = nx.Graph()
                        UNIQUE_BARCODES[ref_length].add_node(CODE)
                        EDGE = get_edge(CODE, list(UNIQUE_BARCODES[ref_length].nodes()), Errors)
                        UNIQUE_BARCODES[ref_length].add_edges_from(EDGE)

            elif (len(POSITIONS_DICT) > 0 and Errors == 0):

                DICTIONARY = POSITIONS_DICT

                for pos2 in DICTIONARY:
                    # printing consensus reads to a new bam file
                    for barcode in DICTIONARY[pos2]:
                        NEW_READ, LOG2 = GET_FINAL_READ(list(DICTIONARY[pos2][barcode]), minBQ, STEP, SET_N)
                        # NEW_READ = DICTIONARY[barcode][0]
                        ##print STEP, NEW_READ.qual

                        LOGFILE2.write(LOG2)
                        outfile.write(NEW_READ)

                POSITIONS_DICT = {}

                pos = ref_start
                end = ref_length

                if ref_length in POSITIONS_DICT:
                    if CODE in POSITIONS_DICT[ref_length]:
                        POSITIONS_DICT[ref_length][CODE].append(read)
                    else:
                        POSITIONS_DICT[ref_length][CODE] = list()
                        POSITIONS_DICT[ref_length][CODE].append(read)
                else:
                    POSITIONS_DICT[ref_length] = {}
                    POSITIONS_DICT[ref_length][CODE] = list()
                    POSITIONS_DICT[ref_length][CODE].append(read)

            else:
                POSITIONS_DICT = {}
                UNIQUE_BARCODES = {}
                EDGE = ''

                pos = ref_start
                end = ref_end

                if ref_length in POSITIONS_DICT:
                    if CODE in POSITIONS_DICT[ref_length]:
                        POSITIONS_DICT[ref_length][CODE].append(read)
                    else:
                        POSITIONS_DICT[ref_length][CODE] = list()
                        POSITIONS_DICT[ref_length][CODE].append(read)
                else:
                    POSITIONS_DICT[ref_length] = {}
                    POSITIONS_DICT[ref_length][CODE] = list()
                    POSITIONS_DICT[ref_length][CODE].append(read)

                # Allowing errors
                if (Errors > 0):
                    if ref_length in UNIQUE_BARCODES:
                        if CODE in list(UNIQUE_BARCODES[ref_length].nodes()):
                            UNIQUE_BARCODES[ref_length].add_node(CODE)
                        else:
                            UNIQUE_BARCODES[ref_length].add_node(CODE)
                            EDGE = get_edge(CODE, list(UNIQUE_BARCODES[ref_length].nodes()), Errors)
                            UNIQUE_BARCODES[ref_length].add_edges_from(EDGE)
                    else:
                        UNIQUE_BARCODES[ref_length] = nx.Graph()
                        UNIQUE_BARCODES[ref_length].add_node(CODE)
                        EDGE = get_edge(CODE, list(UNIQUE_BARCODES[ref_length].nodes()), Errors)
                        # print EDGE
                        UNIQUE_BARCODES[ref_length].add_edges_from(EDGE)

#### We need to print the last groups of reads
if (len(POSITIONS_DICT) > 0 and Errors > 0):
    for pos2 in POSITIONS_DICT:

        # When we allow errors in the Barcodes, we re-group them by similarity (Errors specified in parameter)
        DICTIONARY = extract_bc_groups(POSITIONS_DICT[pos2], UNIQUE_BARCODES[pos2])

        for barcode in DICTIONARY:
            # Printing consensus reads to a new bam file
            NEW_READ, LOG2 = GET_FINAL_READ(list(DICTIONARY[barcode]), minBQ, STEP, SET_N)

            LOGFILE2.write(LOG2)
            outfile.write(NEW_READ)

elif (len(POSITIONS_DICT) > 0 and Errors == 0):

    DICTIONARY = POSITIONS_DICT
    for pos2 in DICTIONARY:

        # printing consensus reads to a new bam file
        for barcode in DICTIONARY[pos2]:
            NEW_READ, LOG2 = GET_FINAL_READ(list(DICTIONARY[pos2][barcode]), minBQ, STEP, SET_N)

            LOGFILE2.write(LOG2)
            outfile.write(NEW_READ)

samfile.close()
LOGFILE2.close()
outfile.close()

stop = timeit.default_timer()
print('TIME')
print(stop - start)
