import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from sensitivity_testing import K_RESULT_PATH, ETH_RESULT_PATH, DIVIDE_BY
from cluster_files import ClusterFile, RAW_DATASETS
from results_files import ResultFile

with open("../paths.json", "r") as f:
    paths = json.load(f)

SENSITIVITY_FIG_FOLDER = paths["figs"]["sensitivity"]

pd.set_option('display.float_format', '{:.20f}'.format)


def get_results_df(dataset_results_dir):
    df = pd.DataFrame(columns=["K", "ETH", "Inaccuracy", "Time"])
    result_paths = map(lambda file_name: os.path.join(dataset_results_dir, file_name), os.listdir(dataset_results_dir))

    for file_path in result_paths:

        raw_path = os.path.join(RAW_CLUSTERS_DIR, os.path.basename(dataset_results_dir) + ".txt")
        curr_path = os.path.dirname(dataset_results_dir)

        # find raw path
        while not os.path.exists(raw_path) and os.path.basename(curr_path) != "noam":
            raw_path = os.path.join(RAW_CLUSTERS_DIR, os.path.basename(curr_path) + ".txt")
            print(raw_path)
            curr_path = os.path.dirname(curr_path)

        # make sure the path was found
        if os.path.basename(curr_path) == "noam":
            print("no raw directory found: " + raw_path)
            exit(1)

        result_file = ResultFile(file_path)
        result_file.complete_df(raw_path)
        length = len(result_file.df.index) // DIVIDE_BY
        inaccuracy = 1 - result_file.df.loc[: length, "classification"].sum() / length
        df.loc[len(df)] = {"K": int(result_file.k), "ETH": int(result_file.eth), "Inaccuracy": inaccuracy, "Time": result_file.total_runtime}

    return df


def sensitivity_graphs(path: str, k: bool = False, eth: bool = False, fontsize=12, fontsize_amp=2):
    if k and eth:
        print("can't use both sensitivity params!")
        return
    if not (k or eth):
        print("needs a sensitivity param!")
        return
    if k:
        column = "K"
        results_dir = K_RESULT_PATH
    else:
        column = "ETH"
        results_dir = ETH_RESULT_PATH
    plt.rcParams.update({'font.size': fontsize})
    title = f"{column} Sensitivity Test"
    fig, axes = plt.subplots(2, 2, figsize=(17, 10))
    fig.suptitle(title, fontsize=fontsize + fontsize_amp * 3)
    axes = axes.flatten()
    for i, dir_name in enumerate(os.listdir(results_dir)):
        df = get_results_df(os.path.join(results_dir, dir_name))
        df = df.sort_values(by=column).reset_index(drop=True)
        print(dir_name)
        print(df)
        ax = axes[i]
        # create Inaccuracy lineplot
        sns.lineplot(x=column, y="Inaccuracy", data=df, label="Inaccuracy", ax=ax, marker="o")
        # fix legend
        ax.legend(loc="upper left")
        ax.set_ylabel("Inaccuracy", fontsize=fontsize + fontsize_amp)
        ax.set_xlabel(column, fontsize=fontsize + fontsize_amp)
        ax.ticklabel_format(style='plain', axis='y')
        ax.locator_params(axis='x', nbins=len(df.index))
        # create Time lineplot with same axis
        ax2 = ax.twinx()
        sns.lineplot(x=column, y='Time', data=df, label='Time', ax=ax2, color='orange', marker="o")
        ax2.set_xlabel(column, fontsize=fontsize + fontsize_amp)
        ax2.lines[0].set_linestyle("--")
        # fix plot
        ax2.legend(loc='upper right')
        ax2.set_ylabel("Time(Seconds)", fontsize=fontsize + fontsize_amp)
        ax2.set_yscale("log")
        # set plot title
        ax.set_title(dir_name, fontsize=fontsize + fontsize_amp * 2)

    plt.ylim(bottom=0)
    plt.tight_layout(pad=2)
    plt.savefig(path)
    plt.clf()


def complex_sensitivity(path: str, result: str, fontsize=12, fontsize_amp=2):
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(17, 10))

    fig.suptitle("ETH Vs. K Comparisons", fontsize=fontsize + fontsize_amp * 3)

    # Flatten axes for easy indexing
    axes = axes.flatten()
    dataset_paths = map(lambda file_name: os.path.join(path, file_name), os.listdir(path))
    for i, dataset_path in enumerate(dataset_paths):
        # Read CSV files in each subfolder
        k_paths = map(lambda file_name: os.path.join(dataset_path, file_name), os.listdir(dataset_path))

        ax = axes[i]
        ax.set_title(os.path.basename(dataset_path))
        combined_df = pd.DataFrame()
        for k_path in k_paths:
            df = get_results_df(k_path)
            df.sort_values(by="ETH").reset_index(drop=True)

            combined_df = pd.concat([combined_df, df])

        # color_palette_inaccuracy = sns.color_palette("Blues", n_colors=len(combined_df['K'].unique()))
        # color_palette_time = sns.color_palette("Oranges", n_colors=len(combined_df['K'].unique()))
        # create Inaccuracy lineplot
        sns.lineplot(x="ETH", y="Inaccuracy", hue="K", palette="bright", data=combined_df, ax=ax, marker="o")
        # fix legend

        for handle in ax.get_legend_handles_labels()[0]:
            handle.set_label(f'K={handle.get_label()}')
        ax.legend(loc="upper left", title="Inaccuracy")
        ax.set_ylabel("Inaccuracy", fontsize=fontsize + fontsize_amp)
        ax.set_xlabel("ETH", fontsize=fontsize + fontsize_amp)
        ax.locator_params(axis='x', nbins=combined_df["ETH"].nunique())
        print(len(combined_df.index))
        print(ax.get_legend_handles_labels())
        # create Time lineplot with same axis
        ax2 = ax.twinx()
        sns.lineplot(x="ETH", y='Time', hue="K", palette="bright", data=combined_df, ax=ax2, marker="o")
        for line in ax2.lines:
            line.set_linestyle("--")
        # fix plot
        for handle in ax2.get_legend_handles_labels()[0]:
            handle.set_label(f'K={handle.get_label()}')
        print(ax2.get_legend_handles_labels())
        ax2.legend(loc='upper right', title="Time")
        ax2.set_ylabel("Time(Seconds)", fontsize=fontsize + fontsize_amp)
        ax2.set_yscale("log")
        ax.set_title(os.path.basename(dataset_path), fontsize=fontsize + fontsize_amp * 2)

    plt.ylim(bottom=0)
    plt.tight_layout(pad=2)
    plt.savefig(result)
    plt.clf()


sensitivity_graphs(os.path.join(SENSITIVITY_FIG_FOLDER, "K.png"), k=True, fontsize=15)
sensitivity_graphs(os.path.join(SENSITIVITY_FIG_FOLDER, "ETH.png"), eth=True, fontsize=15)

complex_sensitivity("/home/noam/alphaProject/results/combined_dir_results", os.path.join(SENSITIVITY_FIG_FOLDER, "complex.png"), fontsize=14, fontsize_amp=1)
