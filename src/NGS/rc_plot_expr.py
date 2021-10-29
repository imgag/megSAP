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
    help="Annotation table with gene IDs and gene symbols.",
)
@click.option(
    "--sample", default=None, help="Sample name of highlighted sample (red cross)."
)
@click.option(
    "--sample2", default=None, help="Sample name of highlighted sample (block square)."
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
def plot(cohort, annotation, sample, sample2, genelist, log2, plot, title, jitter):
    gene_id_names = pd.read_csv(annotation, sep="\t", index_col=0)
    gene_id_names = gene_id_names[["gene_name"]]

    cohort = pd.read_csv(cohort, sep="\t", index_col=0)

    genes = [symbol.strip() for symbol in genelist.readlines()]

    dat = cohort.merge(gene_id_names, left_index=True, right_index=True).reset_index(
        drop=True
    )
    dat_selected = dat[dat["gene_name"].isin(genes)].copy()
    dat_genecol = dat_selected.set_index("gene_name").transpose()
    # column order
    dat_genecol = dat_genecol[genes]

    if log2:
        dat_genecol = dat_genecol.apply(lambda x: np.log2(x + 1))
        ylabel = "log2(TPM+1)"
    else:
        ylabel = "TPM"

    # plot boxplot
    ax1 = dat_genecol.boxplot(grid=False, figsize=(7, 5))
    ax1.set(title=title, ylabel=ylabel)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

    # add highlighted sample
    if sample:
        for i in range(0, len(dat_genecol.columns)):
            ax1.plot(
                i + 1,
                dat_genecol.loc[sample][i],
                marker="x",
                color="red",
                markersize=12,
            )
    # add highlighted sample (sample2)
    if sample2:
        for i in range(0, len(dat_genecol.columns)):
            ax1.plot(
                i + 1,
                dat_genecol.loc[sample2][i],
                marker="s",
                color="black",
                markerfacecolor="none",
                markersize=10,
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
