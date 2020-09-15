#!/usr/bin/env python3

from keras import backend as K
import tensorflow as tf

config = tf.ConfigProto(intra_op_parallelism_threads=1, 
                        inter_op_parallelism_threads=5,
                        allow_soft_placement=True,
                        device_count = {'CPU': 5})
session = tf.Session(config=config)
K.set_session(session)

# Import
import sys
import argparse
from os import stat

from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
from mmsplice.utils import max_varEff, writeVCF

exon_sep = ['[', ']']

def writeTempVCF(vcf_in, vcf_out, dict):

    info_header = "##INFO=<ID=mmsplice,Number=.,Type=String,Description=\"mmsplice splice variant effect: delta_logit_psi(logit scale of variant\'s effect on the exon inclusion, positive score shows higher exon inclusion, negative higher exclusion rate - a score beyond 2, -2 can be considered strong); pathogenicity(Potential pathogenic effect of the variant)\">\n"

    if vcf_in.endswith(".vcf.gz"):
        import gzip
        in_file = gzip.open(vcf_in, 'rt')
    elif vcf_in.endswith(".vcf"):      
        in_file = open(vcf_in, 'r')
    else:
        sys.exit('Wrong file format. Support only \'vcf\' and \'vcf.gz\'')

    out = open(vcf_out, "w")

    for line in in_file.readlines():
        if line.startswith("##"):
            out.write(line)
        elif line.startswith("#CHROM"):
            header = line.split('\t')
            if(not header[7].startswith("INFO")):
                sys.exit(f"Wrong file format. 8th column of vcf file input is no info column: \'{header[7]}\'.")
            out.write(info_header)
            out.write(line)
        else:
            vcf_fields = line.split('\t')
            out.write('\t'.join(vcf_fields[0:7]))

            #generate an ID of CHROM:POS:REF>ALT
            ID = f"{vcf_fields[0]}:{vcf_fields[1]}:{vcf_fields[3]}>{vcf_fields[4]}"
            used_pred = ""

            if ID in dict:
                used_pred = dict[ID]

            if used_pred:
                new_info = vcf_fields[7] + ";" + "mmsplice=" + used_pred
            else:
                new_info = vcf_fields[7]

            if(len(vcf_fields) > 8):
                line_rest = "\t" + new_info + "\t" + "\t".join(vcf_fields[8:len(vcf_fields)])
            else:
                line_rest = "\t" + new_info

            out.write(line_rest)    

def writeMMSpliceToVcf(vcf_in, vcf_out, gtf, fasta):

    # dataloader to load variants from vcf
    dl = SplicingVCFDataloader(gtf, fasta, vcf_in, tissue_specific=False)

    # Specify model
    model = MMSplice()

    # Or predict and return as df
    predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True, progress = True)

    #generate hash
    dict = {}
    for row in predictions.itertuples():
        id = row.ID
        string = exon_sep[0] + "exon:" + row.exons + "," + "delta_logit_psi:" + str(round(row.delta_logit_psi, 4)) + "," + "pathogenicity:" + str(round(row.pathogenicity, 4)) + exon_sep[1]
        if id in dict:
            dict[id] = dict[id] + string
        else:
            dict[id] = string    

    writeTempVCF(vcf_in, vcf_out, dict)

def checkIfEmpty(f, vcf_out):
    if not any(not line.startswith(b"#") for line in f):
        with open(vcf_out, "wb") as out:
            f.seek(0)
            for line in f:
                out.write(line)
        sys.exit()

def main():

    parser = argparse.ArgumentParser(description='Run MMSplice on Input VCF file and write Output VCF file with MMSplice predictions in Info.')
    parser.add_argument('--vcf_in', required=True, dest='vcf_in', help='Input  VCF file.')
    parser.add_argument('--vcf_out', required=True, dest='vcf_out', help='Output VCF file.')
    parser.add_argument('--gtf', required=True, dest='gtf', help='GTF annotation file.')
    parser.add_argument('--fasta', required=True, dest='fasta', help='reference Fasta file.')

    args = ''
    try:
	    args = parser.parse_args()
    except IOError as io:
	    print(io)
	    sys.exit('Error reading parameters.')

    #check if VCF empty (just output input file)
    # we do not want an INFO field mmsplice / mmsplice crashes if no contig line is given
    import gzip
    if args.vcf_in.endswith(".vcf.gz"):
        with gzip.open(args.vcf_in, 'rb') as f:
            checkIfEmpty(f, args.vcf_out)
    elif args.vcf_in.endswith(".vcf"):      
        with open(args.vcf_in, 'rb') as f:
            checkIfEmpty(f, args.vcf_out)
    else:
        sys.exit('Wrong file format. Support only \'vcf\' and \'vcf.gz\'')

    #call tool
    writeMMSpliceToVcf(args.vcf_in, args.vcf_out, args.gtf, args.fasta)

if __name__ == "__main__":
    main()