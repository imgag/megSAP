#!/mnt/share/opt/Python3/bin/python3

# Import
import sys
import argparse
import re 

parser = argparse.ArgumentParser(description='Run MMSplice on Input VCF file and write Output VCF file with MMSplice predictions in Info.')
parser.add_argument('--vcf_in', required=True, dest='vcf_in', help='Input  VCF file.')
parser.add_argument('--vcf_lowAF', required=True, dest='vcf_lowAF', help='VCF file storing low AF variants.')
parser.add_argument('--vcf_out', required=True, dest='vcf_out', help='Output VCF file.')
parser.add_argument('--gtf', required=True, dest='gtf', help='GTF annotation file.')
parser.add_argument('--fasta', required=True, dest='fasta', help='reference Fasta file.')
parser.add_argument('--threads', required=True, dest='threads', help='Number of threads for MMSplice predictions.')

args = ''
try:
    args = parser.parse_args()
except IOError as io:
    print(io)
    sys.exit('Error reading parameters.')

from keras import backend as K
import tensorflow as tf
config = tf.ConfigProto(intra_op_parallelism_threads=1, 
                        inter_op_parallelism_threads=int(args.threads),
                        allow_soft_placement=True,
                        device_count = {'CPU': int(args.threads)})
session = tf.Session(config=config)
K.set_session(session)

from os import stat
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
from mmsplice.utils import max_varEff

exon_sep = ['[', ']']

def lowFrequencyVariants(vcf_in, vcf_lowAF):

    #boolean to check if any variant is left to score with mmsplice
    empty_vcf_file = True

    if vcf_in.endswith(".vcf.gz"):
        import gzip
        in_file = gzip.open(vcf_in, 'rb')
    elif vcf_in.endswith(".vcf"):      
        in_file = open(vcf_in, 'rb')
    else:
        sys.exit('Wrong file format. Support only \'vcf\' and \'vcf.gz\'')

    out = open(vcf_lowAF, "wb")

    AF_id = None
    gnomAD_AF_id = None
    for line in in_file.readlines():
        if line.startswith(b'#'):
            #store all header lines
            out.write(line)
            if(line.startswith(b'##INFO=<ID=CSQ,')):
                pattern = re.compile(b'.*Description=\"(.*)\".*')
                description = re.match(pattern, line)
                entries = line.split(b'|')
                count = 0
                for entry in entries:
                    if(entry.strip()==b'AF'):
                        AF_id = count
                    elif(entry.strip()==b'gnomAD_AF'):
                        gnomAD_AF_id = count
                    count+=1

        else:
            #else store only lines with variants AF<=1%
            vcf_line = line.split(b'\t')
            if(len(vcf_line) < 8):
                sys.exit("Wrong vcf file format of vcf file " + vcf_in + ": Missing INFO column.")
            info = vcf_line[7]
            info = info.split(b';')

            af = 0
            gnomad_af = 0
            gnomad_af_genome = 0
            for info_col in info:
                
                #parse CSQ entry
                if(info_col.startswith(b'CSQ=')):

                    info_col=info_col.split(b'=')
                    csq_annotations=info_col[1].split(b'|')

                    #1000 Genomes AF
                    try:
                        af = float((csq_annotations[AF_id]).decode('ascii'))
                    except:
                        pass
                    #gnomAD exome AF
                    try:
                        gnomad_af = float((csq_annotations[gnomAD_AF_id]).decode('ascii'))
                    except:
                        pass
                #gnomAD genome AF
                elif(info_col.startswith(b'gnomADg_AF=')):
                    try:
                        info_col=info_col.split(b'=')
                        gnomad_af_genome = float((info_col[1]).decode('ascii'))
                    except:
                        pass
            if(af<=0.01 and gnomad_af<=0.01 and gnomad_af_genome<=0.01):
                out.write(line)
                empty_vcf_file = False
    out.close()
    in_file.close()
    return empty_vcf_file

def writeTempVCF(vcf_in, vcf_out, dict):

    info_header = b'##INFO=<ID=mmsplice,Number=.,Type=String,Description=\"mmsplice splice variant effect: delta_logit_psi(logit scale of variant\'s effect on the exon inclusion, positive score shows higher exon inclusion, negative higher exclusion rate - a score beyond 2, -2 can be considered strong); pathogenicity(Potential pathogenic effect of the variant)\">\n'

    if vcf_in.endswith(".vcf.gz"):
        import gzip
        in_file = gzip.open(vcf_in, 'rb')
    elif vcf_in.endswith(".vcf"):      
        in_file = open(vcf_in, 'rb')
    else:
        sys.exit('Wrong file format. Support only \'vcf\' and \'vcf.gz\'')

    out = open(vcf_out, "wb")

    for line in in_file.readlines():
        if line.startswith(b'##'):
            out.write(line)
        elif line.startswith(b'#CHROM'):
            header = line.split(b'\t')
            if(not header[7].startswith(b'INFO')):
                sys.exit("Wrong file format. 8th column of vcf file input is no info column: \'" + header[7].decode('ascii') + "\'.")
            out.write(info_header)
            out.write(line)
        else:
            vcf_fields = line.split(b'\t')
            out.write(b'\t'.join(vcf_fields[0:7]))

            #generate an ID of CHROM:POS:REF>ALT
            ID = vcf_fields[0] + b":" + vcf_fields[1] + b":" + vcf_fields[3] + b">" + vcf_fields[4]
            ID_ascii = ID.decode('ascii')
            used_pred = b""

            if ID_ascii in dict:
                used_pred = dict[ID_ascii]
                used_pred = used_pred.encode('ascii')

            if used_pred:
                new_info = vcf_fields[7] + b";" + b"mmsplice=" + used_pred
            else:
                new_info = vcf_fields[7]

            if(len(vcf_fields) > 8):
                line_rest = b"\t" + new_info + b"\t" + b"\t".join(vcf_fields[8:len(vcf_fields)])
            else:
                line_rest = b"\t" + new_info

            out.write(line_rest)
    out.close()
    in_file.close()

def writeMMSpliceToVcf(vcf_in, vcf_lowAF, vcf_out, gtf, fasta):

    #If VCF file is empty after removel of high AF variants, skip mmsplice annotation
    empty_vcf_file = lowFrequencyVariants(vcf_in, vcf_lowAF)
    if(empty_vcf_file):
        with open(vcf_in, "rb") as in_file:
            with open(vcf_out, "wb") as out_file:
                for line in in_file:
                    out_file.write(line)
        sys.exit()
           
    # dataloader to load variants from vcf
    dl = SplicingVCFDataloader(gtf, fasta, vcf_lowAF, tissue_specific=False)

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
    writeMMSpliceToVcf(args.vcf_in, args.vcf_lowAF, args.vcf_out, args.gtf, args.fasta)

if __name__ == "__main__":
    main()