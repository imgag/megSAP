"""
creates a methylation plot over a given cohort
"""

import sys
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import os.path


def parse_args():
    parser = argparse.ArgumentParser(description="Script to create a methylation plot for a given sample and cohort.")
    parser.add_argument("--plot_data", type=str, help="Path to the cohort data file.", required=True)
    parser.add_argument("--out", type=str, help="output file path of plot.", required=True)
    parser.add_argument("--chr", type=str, help="Chromosome of the plot region.", required=True)
    parser.add_argument("--start", type=int, help="Start base of the plot region.", required=True)
    parser.add_argument("--end", type=int, help="End base of the plot region.", required=True)
    # optional
    parser.add_argument("--site", type=str, help="Name of the methylation site.", required=False, default="")
    parser.add_argument("--sample_name", type=str, help="Name of the sample.", required=False, default="")
    parser.add_argument("--highlight_start", type=int, help="Start base of the highlight region.", required=False,
                        default=-1)
    parser.add_argument("--highlight_end", type=int, help="End base of the highlight region.", required=False,
                        default=-1)
    parser.add_argument("--unphased", action='store_true', help="Creates a unphased plot with only one haplotype.")

    args = parser.parse_args()
    return args


def parse_methylation(filename, args):
    try:
        methylation = pd.read_csv(filename, sep="\t")
        # filter to plot region
        methylation = methylation[(methylation["#chr"] == args.chr)
                                  & (methylation["start"] >= args.start)
                                  & (methylation["end"] <= args.end)]
        if args.unphased:
            #get only *_all columns:
            remaining_columns = [col for col in methylation.columns if col.endswith("_all")]
        else:
            # get only *_hp1/2 columns:
            remaining_columns = [col for col in methylation.columns if col.endswith("_hp1") or col.endswith("_hp2")]
        # slice columns
        methylation = methylation[["end"] + remaining_columns].set_index("end")
        # convert 0-100 -> 0-1
        methylation = methylation / 100.0

        return methylation
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def get_highlight_indices(positions, highlight_start, highlight_end):
    highlight_start_index = 0
    highlight_end_index = -1
    for i in range(0, len(positions)):
        pos = positions[i]
        if highlight_start > pos:
            highlight_start_index += 1
        if highlight_end > pos:
            highlight_end_index += 1
    return highlight_start_index, highlight_end_index


def smooth_data(data, window_size=20):
    data_smoothed = []
    # first step: smoothing by sliding window
    for i in range(0, len(data)):
        window_values = []
        # get window values:
        for j in range(i-int(window_size/2), i+int(window_size/2)):
            # skip indices outside of data
            if j < 0 or j >= len(data):
                continue
            window_values.append(data[j])
        # get mean value
        data_smoothed.append(np.mean(window_values))

    # second step: smoothing via pandas rolling
    return pd.DataFrame(data_smoothed).rolling(20, win_type="hamming").mean()


def main():
    # get params
    args = parse_args()

    # check parameters
    # check plot data file
    if not os.path.isfile(args.plot_data):
        raise FileNotFoundError("Plot data file '" + args.plot_data + "' not found!")
    # check intervals
    if args.start < 0 or args.end < 0:
        raise ValueError("Coordinates must be greater than 0!")
    if args.start > args.end or args.highlight_start > args.highlight_end:
        raise ValueError("Start coordinates must be smaller than end coordinates!")
    if (0 < args.highlight_start < args.start) or (args.highlight_end > 0 and args.highlight_end > args.end):
        raise ValueError("Highlight coordinates must be contained in plot coordinates!")

    methylation = parse_methylation(args.plot_data, args)

    # get cohort columns
    if args.unphased:
        cohort_columns_all = [col for col in methylation.columns if col.endswith("_all") and col != "index_all"]
    else:
        cohort_columns_hp1 = [col for col in methylation.columns if col.endswith("_hp1") and col != "index_hp1"]
        cohort_columns_hp2 = [col for col in methylation.columns if col.endswith("_hp2") and col != "index_hp2"]

    # get highlight region
    if args.highlight_start > 0 and args.highlight_end > 0:
        plot_highlight = True
        positions = methylation.index.to_list()
        highlight_start_index, highlight_end_index = get_highlight_indices(positions, args.highlight_start, args.highlight_end)
    else:
        plot_highlight = False

    # smooth data
    if len(methylation.index) > 0:
        methylation.reset_index(inplace=True)
        for col in methylation.columns:
            if col == "end":
                continue
            methylation[col] = smooth_data(methylation[col].to_list())
        #methylation.set_index("end", inplace=True)

    # drop CpG sites which are not covered in index
    if args.unphased:
        methylation = methylation.dropna(subset=["index_all"])
    else:
        methylation = methylation.dropna(subset=["index_hp1", "index_hp2"], how="all")

    # create plots
    colors = sns.color_palette("tab10")
    if not args.unphased:
          fig, axes = plt.subplots(2, 1)
          fig.set_size_inches(13.4, 6.4)
    else:
        fig, axes = plt.subplots(1, 1)
        fig.set_size_inches(13.4, 3.2)
        axes = [axes]


    fig.suptitle("Cohort plot of " + args.site + " (" + args.sample_name + ")")

    if args.unphased:
        # plot combined plot
        if plot_highlight:
            axes[0].add_patch(patches.Rectangle((highlight_start_index, -0.1), (highlight_end_index - highlight_start_index), 1.2, edgecolor=colors[1], facecolor=colors[1], alpha=0.1))
        axes[0].set_ylim(-0.1, 1.1)
        for col in cohort_columns_all:
            sns.lineplot(methylation[col], ax=axes[0], color="grey", alpha=0.2)
        if "index_all" in methylation.columns:
            sns.lineplot(methylation["index_all"], ax=axes[0], color=colors[0], linewidth=3)
        axes[0].set(xlabel=None, ylabel=None)

    else:
        # plot 2 haplotypes
        # hp1
        axes[0].set_title("allele 1")
        if plot_highlight:
            axes[0].add_patch(patches.Rectangle((highlight_start_index, -0.1), (highlight_end_index - highlight_start_index), 1.2, edgecolor=colors[1], facecolor=colors[1], alpha=0.1))
        axes[0].set_ylim(-0.1, 1.1)
        for col in cohort_columns_hp1:
            sns.lineplot(methylation[col], ax=axes[0], color="grey", alpha=0.2)
        sns.lineplot(methylation["index_hp1"], ax=axes[0], color=colors[0], linewidth=3)
        axes[0].set(xlabel=None, ylabel=None)
        # hp2
        axes[1].set_title("allele 2")
        if plot_highlight:
            axes[1].add_patch(patches.Rectangle((highlight_start_index, -0.1), (highlight_end_index - highlight_start_index), 1.2, edgecolor=colors[0], facecolor=colors[0], alpha=0.1))
        axes[1].set_ylim(-0.1, 1.1)
        for col in cohort_columns_hp2:
            sns.lineplot(methylation[col], ax=axes[1], color="grey", alpha=0.2)
        sns.lineplot(methylation["index_hp2"], ax=axes[1], color=colors[1], linewidth=3)
        axes[1].set(xlabel=None, ylabel=None)

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
    plt.savefig(args.out, format="png")
    # plt.show()


if __name__ == '__main__':
    sys.exit(main())

