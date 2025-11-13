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
    parser.add_argument("--hp1_file", type=str, help="modkit file containing the sample's methylation of the first allele.", required=True)
    parser.add_argument("--out", type=str, help="output file path of plot.", required=True)
    # optional
    parser.add_argument("--hp2_file", type=str, help="modkit file containing the sample's methylation of the second allele.", required=False, default="")
    parser.add_argument("--hp1_cohort_files", nargs='+', type=str, help="modkit files containing the cohort's methylation of the first allele (continously phased).",
                        required=False, default=[])
    parser.add_argument("--hp2_cohort_files", nargs='+', type=str, help="modkit files containing the cohort's methylation of the second allele (continously phased).",
                        required=False, default=[])
    parser.add_argument("--hp1_cohort_files_fallback", nargs='+', type=str, help="modkit files containing the cohort's methylation of the first allele (not continously phased).",
                        required=False, default=[])
    parser.add_argument("--hp2_cohort_files_fallback", nargs='+', type=str, help="modkit files containing the cohort's methylation of the second allele (not continously phased).",
                        required=False, default=[])
    parser.add_argument("--site", type=str, help="Name of the methylation site.", required=False, default="")
    parser.add_argument("--sample_name", type=str, help="Name of the sample.", required=False, default="")
    parser.add_argument("--highlight_start", type=int, help="Start base of the highlight region.", required=False, default=-1)
    parser.add_argument("--highlight_end", type=int, help="End base of the highlight region.", required=False, default=-1)
    parser.add_argument("--unphased", action='store_true', help="Creates a unphased plot with only one haplotype.")

    args = parser.parse_args()
    return args


def parse_methylation(filename, sample_name="", column_suffix=""):
    if sample_name == "":
        # extract sample name from file name
        sample_name = "_".join(str(os.path.basename(filename)).split("_")[0:2])

    try:
        methylation = pd.read_csv(filename, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    # slice to keep only methyl fraq and pos
    methylation["position"] = methylation.iloc[:, 2].astype(int)
    methylation = methylation.set_index("position")
    # methylation["pos"] = methylation.iloc[:, 2].astype(int)
    methylation = methylation.iloc[:, 10].astype(float).rename(sample_name + column_suffix) / 100.0
    return methylation


def get_highlight_indices(positions, highlight_start, highlight_end):
    highlight_start_index = 0
    highlight_end_index = -1
    for i in range(0, len(positions)):
        pos = positions[i]
        if highlight_start > pos:
            highlight_start_index = i + 1
        if highlight_end > pos:
            highlight_end_index = i + 1
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
    if not args.unphased and args.hp2_file == "":
        raise ValueError("Either 'hp2_file' or 'unphased' must be provided!")

    # check files
    if not os.path.isfile(args.hp1_file):
        raise FileNotFoundError("hp1_file '" + args.hp1_file + "' not found!")
    for file in args.hp1_cohort_files:
        if not os.path.isfile(file):
            raise FileNotFoundError("cohort file '" + file + "' not found!")
    for file in args.hp1_cohort_files_fallback:
        if not os.path.isfile(file):
            raise FileNotFoundError("cohort file '" + file + "' not found!")

    if not args.unphased:
        if not os.path.isfile(args.hp2_file):
            raise FileNotFoundError("hp2_file '" + args.hp2_file + "' not found!")
        for file in args.hp2_cohort_files:
            if not os.path.isfile(file):
                raise FileNotFoundError("cohort file '" + file + "' not found!")
        for file in args.hp2_cohort_files_fallback:
            if not os.path.isfile(file):
                raise FileNotFoundError("cohort file '" + file + "' not found!")

    # get sample methylation
    hp1_methylation = pd.DataFrame(parse_methylation(args.hp1_file, "case_hp1"))
    if args.unphased:
        methylation = hp1_methylation
    else:
        hp2_methylation = pd.DataFrame(parse_methylation(args.hp2_file, "case_hp2"))
        # create combined DataFrame to keep plots in sync
        methylation = pd.concat([hp1_methylation, hp2_methylation], axis=1)

    # get cohort methylation
    for file in args.hp1_cohort_files:
        methylation = methylation.join(parse_methylation(file, column_suffix="_hp1"))
    for file in args.hp1_cohort_files_fallback:
        methylation = methylation.join(parse_methylation(file, column_suffix="_hp1_fallback"))
    if not args.unphased:
        for file in args.hp2_cohort_files:
            methylation = methylation.join(parse_methylation(file, column_suffix="_hp2"))
        for file in args.hp2_cohort_files_fallback:
            methylation = methylation.join(parse_methylation(file, column_suffix="_hp2_fallback"))
    methylation.sort_index(inplace=True)

    # get highlight region
    if args.highlight_start > 0 and args.highlight_end > 0:
        plot_highlight = True
        positions = methylation.index.to_list()
        highlight_start_index, highlight_end_index = get_highlight_indices(positions, args.highlight_start, args.highlight_end)
    else:
        plot_highlight = False

    # update index
    methylation["index"] = range(1, methylation.shape[0] + 1)
    methylation = methylation.set_index("index")

    # smooth data
    if len(methylation.index) > 0:
        for col in methylation.columns:
            methylation[col] = smooth_data(methylation[col].to_list())

    # get hp1/hp2 columns
    columns_hp1 = []
    columns_hp2 = []
    columns_hp1_fallback = []
    columns_hp2_fallback = []
    for col in methylation.columns:
        if not col.startswith("case_"):
            if col.endswith("_hp1"):
                columns_hp1.append(col)
            if col.endswith("_hp2"):
                columns_hp2.append(col)
            if col.endswith("_hp1_fallback"):
                columns_hp1_fallback.append(col)
            if col.endswith("_hp2_fallback"):
                columns_hp2_fallback.append(col)

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

    if not args.unphased:
        axes[0].set_title("allele 1")
    if plot_highlight:
        axes[0].add_patch(patches.Rectangle((highlight_start_index, -0.1), (highlight_end_index - highlight_start_index), 1.2, edgecolor=colors[1], facecolor=colors[1], alpha=0.1))
    axes[0].set_ylim(-0.1, 1.1)
    for col in columns_hp1:
        sns.lineplot(methylation[col], ax=axes[0], color="grey", alpha=0.2)
    for col in columns_hp1_fallback:
        sns.lineplot(methylation[col], ax=axes[0], color="grey", alpha=0.1, linestyle='dashed')
    if "case_hp1" in methylation.columns:
        sns.lineplot(methylation["case_hp1"], ax=axes[0], color=colors[0], linewidth=3)
    axes[0].set(xlabel=None, ylabel=None)

    if not args.unphased:
        axes[1].set_title("allele 2")
        if plot_highlight:
            axes[1].add_patch(patches.Rectangle((highlight_start_index, -0.1), (highlight_end_index - highlight_start_index), 1.2, edgecolor=colors[0], facecolor=colors[0], alpha=0.1))
        axes[1].set_ylim(-0.1, 1.1)
        for col in columns_hp2:
            sns.lineplot(methylation[col], ax=axes[1], color="grey", alpha=0.2)
        for col in columns_hp2_fallback:
            sns.lineplot(methylation[col], ax=axes[1], color="grey", alpha=0.1, linestyle='dashed')
        if "case_hp2" in methylation.columns:
            sns.lineplot(methylation["case_hp2"], ax=axes[1], color=colors[1], linewidth=3)
        axes[1].set(xlabel=None, ylabel=None)

    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
    plt.savefig(args.out, format="png")

if __name__ == '__main__':
    sys.exit(main())

