import sys
import os
import gzip
from optparse import OptionParser
import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_seq


def main(argv):
    if len(argv) == 1:
        argv.append("--help")

    usage = "usage: %prog -i counts.tsv -b batch.tsv -o batch_corrected.tsv"
    desc = "Batch correction of gene or exon expression. Input file must contain a matrix where the first column is the index (gene, exon) and any other column represents a sample. The batch file needs a column RNA and COHORT."
    
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("-i", "--in", action="store", dest="input", type="string", help="Input gene expression matrix. First column must be the gene/exon and any other column  the counts per sample.")
    parser.add_option("-b", "--batch", action="store", dest="batch", type="string", help="TSV file with batch information. Must contain a column #RNA and a column COHORT.")
    parser.add_option("-o", "--out", action="store", dest="out", type="string", help="gzipped output file with batch corrected counts.")
    #optinal
    parser.add_option("-c", "--covariants", action="store", dest="covar", type="string", default="", help="TSV file with covariant information. Must contain a column RNA and a column COVARIANT.")

    (options, args) = parser.parse_args()

    # Output directory
    odir = os.path.dirname(os.path.abspath(options.out))

    # Read expression matrix
    df_expression = pd.read_csv(options.input, delimiter="\t", index_col=0)
    
    # Read batch and extract batch information
    batch = pd.read_csv(options.batch, delimiter="\t", index_col=0)
    batch = batch.reindex(df_expression.columns)
    batch_list = batch["BATCH"].to_list()
    
    print(batch_list)
    
    covar_list = None
    if options.covar != "":
        covar = pd.read_csv(options.covar, delimiter="\t", index_col=0)
        covar = covar.reindex(df_expression.columns)
        covar_list = covar["COVARIANT"].to_list()
    
        # if there is only one covariant group - ignore it
        if (len(np.unique(covar_list)) == 1):
            covar_list = None
    
    print(covar_list)
    
    # Batch correction
    batch_corrected = pycombat_seq(df_expression, batch_list, covar_mod=covar_list)
    
    #write gzipped:
    with gzip.open(options.out, "wb") as out:
        batch_corrected.to_csv(out, sep="\t")
    

if __name__ == "__main__":
    main(sys.argv)