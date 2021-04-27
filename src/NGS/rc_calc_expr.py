#!/usr/bin/env python3

import pandas as pd
import numpy as np
import click

def read_counts_tpm_megsap(fp, sample_id):
    """Read TPM column from megSAP count file."""
    df = pd.read_csv(fp, sep='\t', usecols=['#gene_id','tpm'], index_col='#gene_id')
    df.index.name = 'gene_id'
    df.columns = [sample_id]
    return df

@click.group()
def cli():
    pass

@cli.command(help="Use expression data from cohort samples as background.")
@click.option('--counts', type=click.File('r'), help="Sample gene counts for which to calculate relative expression.")
@click.option('--counts_out', type=click.File('w'), help="Output sample gene counts with relative expression calculated.")
@click.option('--prefix', default="cohort_", help="Column name prefix for new columns.", show_default=True)
@click.option('--samples', type=click.File('r'), required=True, help="Sample table with sample ID and file names.")
@click.option('--cohort', type=click.File('w'), help="Output file for gene-sample expression data.")
@click.option('--stats', type=click.File('w'), help="Output file for per-gene statistics.")
def cohort(samples, counts, counts_out, prefix, cohort, stats):
    # read list of file names
    dat = pd.read_csv(samples, sep='\t', header=None, names=['psample', 'path'])

    # genes-samples matrix
    dfs = [read_counts_tpm_megsap(row.path, row.psample) for _, row in dat.iterrows()]
    d_tpm = dfs[0].join(dfs[1:], how='outer')

    if cohort:
        d_tpm.to_csv(cohort, sep='\t')

    # statistics per gene
    df_stats = pd.concat([
        d_tpm.mean(axis=1),
        np.log2(d_tpm+1).mean(axis=1),
        np.log2(d_tpm+1).std(axis=1) ],
        axis=1,
        keys=["mean", "meanlog2", "sdlog2"]
        )

    if stats:
        df_stats.to_csv(stats, sep='\t')
    
    if counts:
        counts_tbl =  pd.read_csv(counts, sep="\t", index_col=0)
        counts_tbl.drop([col for col in counts_tbl.columns if col.startswith(prefix)], axis=1, inplace=True)
        # calculate ratios per gene for sample
        expr = counts_tbl['tpm'].to_frame().join(df_stats)
        expr['log2tpm'] = np.log2(expr['tpm']+1)
        expr['log2fc'] = expr['log2tpm'] - expr['meanlog2']
        expr['zscore'] = (expr['log2tpm'] - expr['meanlog2'])/expr['sdlog2']

        expr.drop('tpm', axis=1, inplace=True)
        add_prefix = lambda s: prefix + str(s)
        expr.rename(columns=add_prefix, inplace=True)

        # annotated count file
        counts_annot = counts_tbl.join(expr, on="#gene_id", how="left")

        if counts_out:
            counts_annot.to_csv(counts_out, sep='\t', na_rep='n/a')

@cli.command(help="Use summarized HPA expression data as reference.")
@click.option('--counts', type=click.File('r'), help="Sample gene counts for which to calculate relative expression.")
@click.option('--counts_out', type=click.File('w'), help="Output sample gene counts with relative expression calculated.")
@click.option('--prefix', default="hpa_", help="Column name prefix for new columns.", show_default=True)
@click.option('--hpa', type=click.File('r'), required=True,
    default="/mnt/storage1/share/data/dbs/gene_expression/rna_tissue_hpa.tsv",
    help="HPA tissue expression data.")
@click.option('--tissue', required=True, help="HPA reference tissue")
def hpa(counts, counts_out, prefix, hpa, tissue):
    # read HPA data
    dat = pd.read_csv(hpa, sep='\t', index_col=0)

    # genes-samples matrix
    d_tpm = dat[dat['Tissue'] == tissue][['TPM']]
    d_tpm.rename(columns={'TPM':'tissue_tpm'}, inplace=True)
    d_tpm['meanlog2'] = np.log2(d_tpm['tissue_tpm']+1)
    
    if counts:
        counts_tbl =  pd.read_csv(counts, sep="\t", index_col=0)
        counts_tbl.drop([col for col in counts_tbl.columns if col.startswith(prefix)], axis=1, inplace=True)
        # calculate ratios per gene for sample
        expr = counts_tbl['tpm'].to_frame().join(d_tpm[['tissue_tpm','meanlog2']])
        expr['log2tpm'] = np.log2(expr['tpm']+1)
        expr['log2fc'] = expr['log2tpm'] - expr['meanlog2']

        expr.drop('tpm', axis=1, inplace=True)
        add_prefix = lambda s: prefix + str(s)
        expr.rename(columns=add_prefix, inplace=True)

        # annotated count file
        counts_annot = counts_tbl.join(expr, on="#gene_id", how="left")

        if counts_out:
            counts_annot.to_csv(counts_out, sep='\t', na_rep='n/a')

if __name__ == '__main__':
  cli()
