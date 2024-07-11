"""
Small python script to extract read counts from illumina NovaSeq X demultiplexing and display them in a table and bar chart.

"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Script to extract read counts and display them in a table and bar chart")
    parser.add_argument("--demux_stats", type=str, help="CSV file containing the demux stats (Demultiplex_Stats.csv).",
                        default="Analysis/1/Data/Demux/Demultiplex_Stats.csv", required=False)
    parser.add_argument("--unknown_barcodes", type=str, help="CSV file containing the top unknown barcodes (Top_Unknown_Barcodes.csv).",
                        default="Analysis/1/Data/Demux/Top_Unknown_Barcodes.csv", required=False)
    parser.add_argument("--table", type=str, help="HTML file containing the counts as table.", default="DemuxStats.html", required=False)
    parser.add_argument("--plot", type=str, help="PNG file containing the counts as bar plot.", default="DemuxStats.png", required=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    demux_stats = pd.read_csv(args.demux_stats)
    demux_stats["Sample"] = demux_stats["SampleID"]
    tmp_index = demux_stats[["Lane", "SampleID", "Sample", "Index"]].set_index(["Lane", "SampleID"]).unstack(level=0)[("Index", 1)]
    tmp_sample = demux_stats[["Lane", "SampleID", "Sample", "Index"]].set_index(["Lane", "SampleID"]).unstack(level=0)[("Sample", 1)]
    demux_stats = demux_stats[["Lane", "Sample", "# Reads"]].set_index(["Lane", "Sample"]).unstack(level=0)
    demux_stats.columns = ["# Reads (lane " + str(x[1]) + ")" for x in demux_stats.columns]
    demux_stats["# Reads (total)"] = demux_stats.sum(axis=1)
    demux_stats.insert(0, "Index", tmp_index)
    demux_stats.insert(0, "Sample", tmp_sample)

    unknown_barcodes = pd.read_csv(args.unknown_barcodes)
    unknown_barcodes["Index"] = unknown_barcodes["index"] + "-" + unknown_barcodes["index2"]
    unknown_barcodes = unknown_barcodes[["Lane", "Index", "# Reads"]].set_index(["Lane", "Index"]).unstack(level=0).fillna(0).astype(int)
    unknown_barcodes.columns = ["# Reads (lane " + str(x[1]) + ")" for x in unknown_barcodes.columns]
    unknown_barcodes["# Reads (total)"] = unknown_barcodes.sum(axis=1)
    unknown_barcodes = unknown_barcodes.sort_values("# Reads (total)", ascending=False).reset_index()
    unknown_barcodes.insert(0, "Sample", "")

    combined_table = pd.concat([demux_stats, unknown_barcodes])

    # combined_table.drop("# Reads (total)", axis=1).drop("Undetermined", axis=0).plot(kind="bar", stacked=True)

    # format table
    def format_read_count(x):
        return f'{x:,}'

    formatter = {}
    # add formatter for each read count column
    for col_name in combined_table.columns:
        if col_name.startswith("# Reads"):
            formatter[col_name] = format_read_count

    raw_html = combined_table.to_html(index=False, formatters=formatter, justify='center', classes='table table-striped text-center').replace('<td>', '<td align="right">')

    with open(args.table, 'w') as output_table:
        output_table.write(raw_html)

    print("Table generated!")

    def get_excluded_barcodes(row):
        idx1_len = len(row.split('-')[0])
        idx2_len = len(row.split('-')[1])

        if ("G" * idx1_len) in row or ("G" * idx2_len) in row:
            return True
        return False

    # cleanup data for plotting
    demux_stats = demux_stats.set_index("Sample").drop("# Reads (total)", axis=1).drop("Index", axis=1).drop("Undetermined", axis=0)
    unknown_barcodes = unknown_barcodes.iloc[:50, ].drop("# Reads (total)", axis=1).set_index("Index")
    unknown_barcodes["excluded"] = unknown_barcodes.index.to_series().apply(get_excluded_barcodes)

    included_barcodes = unknown_barcodes.copy()
    included_barcodes[included_barcodes["excluded"]] = np.nan
    excluded_barcodes = unknown_barcodes.copy()
    excluded_barcodes[excluded_barcodes["excluded"] == False] = np.nan

    lanes = demux_stats.columns
    colors_samples = sns.color_palette("crest", len(lanes))
    colors_barcodes = sns.color_palette("flare", len(lanes))
    colors_barcodes_exclude = sns.color_palette('Greys', len(lanes) + 4)[2:-2]

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(25, 12), sharey=True)
    fig.suptitle("Demultiplexing Statistics", fontsize=25)

    ax1 = demux_stats.plot(kind="bar", stacked=True, color=colors_samples, subplots=False, ax=axes[0], legend=False)
    ax1.tick_params(axis='x', labelsize=10)
    plt.xlabel("Samples")

    ax2 = included_barcodes.plot(kind="bar", stacked=True, color=colors_barcodes, subplots=False, ax=axes[1], legend=False)
    ax2 = excluded_barcodes.drop("excluded", axis=1).plot(kind="bar", stacked=True, color=colors_barcodes_exclude, subplots=False, ax=axes[1])
    ax2.tick_params(axis='x', labelsize=8)
    plt.xlabel("Unknown barcodes")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.savefig(args.plot, dpi=300, format="png")

    print("Plot generated!")


if __name__ == '__main__':
    main()

