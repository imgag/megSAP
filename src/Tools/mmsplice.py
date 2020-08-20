#!/usr/bin/env python3

# Import
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table, writeVCF
from mmsplice.utils import max_varEff

exon_sep = ['[', ']']

def writeTempVCF(vcf_in, vcf_out, predictions):

    from cyvcf2 import Writer, VCF
    vcf =  VCF(vcf_in)
    vcf.add_info_to_header({
            'ID': 'mmsplice',
            'Description': 'mmsplice splice variant effect',
            'Type': 'Character',
            'Number': '.'
    })

    w = Writer(vcf_out, vcf)
    for var in vcf:

        #generate an ID of CHROM:POS:REF>ALT
        alt = ','.join(var.ALT)
        ID = f"{var.CHROM}:{var.POS}:{var.REF}>{alt}"

        #get list of indices containing mmsplice info
        indices =  predictions.ID[predictions.ID == ID].index.tolist()

        #extract useful information of mmsplice
        used_pred = ""
        for index in indices:
            pred = predictions.iloc[index]
            used_pred += exon_sep[0] + "exon:" + pred.exons + "," + "delta_logit_psi:" + str(pred.delta_logit_psi) + "," + "pathogenicity:" + str(pred.pathogenicity) + exon_sep[1]

        #add information to vcf file

        if used_pred is not None:
            var.INFO['mmsplice'] = used_pred
        w.write_record(var)

def writeMMSpliceToVcf(vcf_in, vcf_out):

    gtf = '/mnt/users/ahstoht1/Splicing/Homo_sapiens_GRCh37_75/Homo_sapiens_GRCh37_75_uniq_exon.gtf.gz'
    fasta = '/mnt/share/data/genomes/GRCh37.fa'

    # dataloader to load variants from vcf
    dl = SplicingVCFDataloader(gtf, fasta, vcf_in, tissue_specific=False)

    # Specify model
    model = MMSplice()

    # Or predict and return as df
    predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True, progress = True)

    # Summerize with maximum effect size
    writeTempVCF(vcf_in, vcf_out, predictions)

def main():

    vcf_in = sys.argv[1]
    vcf_out = sys.argv[2]

    writeMMSpliceToVcf(vcf_in, vcf_out)

if __name__ == "__main__":
    main()