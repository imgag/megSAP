import pandas as pd
import numpy as np
import click
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@click.command(help="Plot cohort expression data for selected genes as boxplots.")
@click.option(
    "--cohort",
    type=click.File("r"),
    help="Matrix file with gene-sample TPM expression data.",
)
@click.option(
    "--annotation",
    type=click.File("r"),
    help="Annotation table with gene IDs, gene symbols and reference expression.",
)
@click.option(
    "--sample", default=None, help="Sample name of highlighted sample (red cross)."
)
@click.option(
    "--reference/--no-reference",
    default=False,
    help="Show reference expression (black square).",
    show_default=True,
)
@click.option(
    "--genelist",
    type=click.File("r"),
    help="File with list of genes to plot, one gene symbol per line.",
)
@click.option(
    "--log2/--no-log2",
    default=True,
    help="Use log2(TPM+1) transformation.",
    show_default=True,
)
@click.option("--plot", help="Plot output file.")
@click.option(
    "--title", default="Gene Expression", help="Plot title.", show_default=True
)
@click.option(
    "--jitter/--no-jitter",
    default=True,
    help="Show individual data points.",
    show_default=True,
)
@click.option(
    "--violin/--no-violin",
    default=False,
    help="Use violin plot instead boxplot.",
    show_default=True,
)
@click.option(
    "--outlier/--no-outlier",
    default=False,
    help="Show outliers of boxplot as black circles. Only has an effect if violin plot is not active.",
    show_default=True,
)
def plot(
    cohort, annotation, sample, reference, genelist, log2, plot, title, jitter, violin, outlier
):
    sample_expr = pd.read_csv(annotation, sep="\t", index_col=0)
    gene_id_names = sample_expr[["gene_name"]]

    cohort = pd.read_csv(cohort, sep="\t", index_col=0)

    genes = [symbol.strip() for symbol in genelist.readlines()]

    reference = reference and ("hpa_tissue_tpm" in sample_expr.columns)
    if reference:
        ref_expr = sample_expr[["gene_name", "hpa_tissue_tpm"]]
        ref_expr = ref_expr[ref_expr["gene_name"].isin(genes)]
        ref_expr = ref_expr.set_index("gene_name")
        ref_expr = ref_expr.loc[genes]

    dat = cohort.merge(gene_id_names, left_index=True, right_index=True).reset_index(
        drop=True
    )
    dat_selected = dat[dat["gene_name"].isin(genes)].copy()
    dat_genecol = dat_selected.set_index("gene_name").transpose()
    # column order
    dat_genecol = dat_genecol[genes]

    if log2:
        dat_genecol = dat_genecol.apply(lambda x: np.log2(x + 1))
        if reference:
            ref_expr = ref_expr.apply(lambda x: np.log2(x + 1))
        ylabel = "log2(TPM+1)"
    else:
        ylabel = "TPM"

    # plot boxplot

    plt.figure(figsize=(7, 5))
    fig, ax1 = plt.subplots()
    if not violin:
        dat_genecol.boxplot(grid=False, ax=ax1, showfliers=outlier)
    else:
        ax1.violinplot(
            dataset=dat_genecol.transpose(),
            widths=0.8,
            showmedians=True,
            showextrema=False,
        )
        labels = dat_genecol.columns
        ax1.set_xticks(np.arange(1, len(labels) + 1))
        ax1.set_xticklabels(labels)

    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
    ax1.set(title=title, ylabel=ylabel)

    # add highlighted sample
    if sample:
        for i in range(0, len(dat_genecol.columns)):
            ax1.plot(
                i + 1,
                dat_genecol.loc[sample].iloc[i],
                marker="x",
                color="red",
                markersize=12,
            )
    # add reference expression
    if reference:
        for i in range(0, len(dat_genecol.columns)):
            ax1.plot(
                i + 1,
                ref_expr.iloc[i],
                marker="s",
                color="black",
                markerfacecolor="none",
                markersize=6,
            )

    # add all data points as jitter points
    if jitter:
        for i in range(0, len(dat_genecol.columns)):
            pts = dat_genecol.iloc[:, i]
            x = np.random.normal(i + 1, 0.04, size=len(pts))
            ax1.scatter(x, pts, color="grey", s=10, alpha=0.5)

    # save figure
    plt.savefig(plot, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    np.random.seed(0)
    plot()
