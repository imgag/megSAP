#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import norm
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
@click.option('--prefix', default="", help="Column name prefix for new columns.", show_default=True)
@click.option('--samples', type=click.File('r'), required=True, help="Sample table with sample ID and file names.")
@click.option('--cohort', type=click.File('w'), help="Output file for gene-sample expression data.")
@click.option('--stats', type=click.File('w'), help="Output file for per-gene statistics.")
@click.option('--corr', type=click.File('w'), help="Output file for sample--cohort correlation.")
def cohort(samples, counts, counts_out, prefix, cohort, stats, corr):
    # read list of file names
    dat = pd.read_csv(samples, sep='\t', header=None, names=['psample', 'path'])

    # genes-samples matrix
    dfs = [read_counts_tpm_megsap(row.path, row.psample) for _, row in dat.iterrows()]
    d_tpm = dfs[0].join(dfs[1:], how='outer')

    if cohort:
        cohort.write('#')
        d_tpm.to_csv(cohort, sep='\t')

    # statistics per gene
    df_stats = pd.concat([
        d_tpm.mean(axis=1),
        np.log2(d_tpm+1).mean(axis=1),
        np.log2(d_tpm+1).std(axis=1) ],
        axis=1,
        keys=["cohort_mean", "cohort_meanlog2", "cohort_sdlog2"]
        )

    if stats:
        stats.write('#')
        df_stats.to_csv(stats, sep='\t')
    
    if counts:
        counts_tbl =  pd.read_csv(counts, sep="\t", index_col=0)
        # calculate ratios per gene for sample
        expr = counts_tbl[['tpm']].join(df_stats)
        expr['log2tpm'] = np.log2(expr['tpm']+1)
        expr["log2fc"] = np.log2(expr["tpm"] + 1) - np.log2(expr["cohort_mean"] + 1)
        expr['zscore'] = (np.log2(expr["tpm"] + 1) - expr['cohort_meanlog2'])/expr['cohort_sdlog2']
        expr['pval'] = expr['zscore'].apply(lambda z: 2*norm.cdf(-abs(z)))

        # annotated count file
        if counts_out:
            counts_annot = counts_tbl.join(expr.drop('tpm', axis=1).rename(columns=lambda s: prefix + str(s)),
                how="left")
            counts_annot.to_csv(counts_out, sep='\t', na_rep='n/a')
        
        if corr:
            value = expr[['log2tpm', 'cohort_meanlog2']].query('log2tpm>0 | cohort_meanlog2>0').corr(method='spearman').iloc[0,1]
            corr.write(f'{value}\n')

@cli.command(help="Use summarized HPA expression data as reference.")
@click.option('--counts', type=click.File('r'), help="Sample gene counts for which to calculate relative expression.")
@click.option('--counts_out', type=click.File('w'), help="Output sample gene counts with relative expression calculated.")
@click.option('--prefix', default="hpa_", help="Column name prefix for new columns.", show_default=True)
@click.option('--hpa', type=click.File('r'), required=True,
    default="/mnt/storage2/megSAP/data/dbs/gene_expression/rna_tissue_consensus_v22.tsv",
    help="HPA tissue expression data.")
@click.option('--tissue', required=True, help="HPA reference tissue")
@click.option('--corr', type=click.File('w'), help="Output file for sample--cohort correlation.")
def hpa(counts, counts_out, prefix, hpa, tissue, corr):
    # read HPA data
    dat = pd.read_csv(hpa, sep='\t', index_col=0)

    # genes-samples matrix
    ref = dat[dat['Tissue'] == tissue][['TPM']]
    ref.rename(columns={'TPM':'tissue_tpm'}, inplace=True)
    ref['tissue_log2tpm'] = np.log2(ref['tissue_tpm']+1)
    
    if counts:
        counts_tbl =  pd.read_csv(counts, sep="\t", index_col=0)
        # calculate ratios per gene for sample
        expr = counts_tbl[['tpm']].join(ref[['tissue_tpm','tissue_log2tpm']])
        expr['sample_log2tpm'] = np.log2(expr['tpm']+1)
        expr['log2fc'] = expr['sample_log2tpm'] - expr['tissue_log2tpm']

        # annotated count file
        if counts_out:
            counts_annot = counts_tbl.join(expr.drop('tpm', axis=1).rename(columns=lambda s: prefix + str(s)),
                how="left")
            counts_annot.to_csv(counts_out, sep='\t', na_rep='n/a')
        
        if corr:
            value = expr[['sample_log2tpm', 'tissue_log2tpm']].query('sample_log2tpm>0 | tissue_log2tpm>0').corr(method='spearman').iloc[0,1]
            corr.write(f'{value}\n')

if __name__ == '__main__':
  cli()
