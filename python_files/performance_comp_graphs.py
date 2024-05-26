import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from cluster_files import get_files, RAW_DATASETS
from results_files import ResultFile
from performance_comp_testing import PERFORMANCE_RESULTS_PATH

with open("../paths.json", "r") as f:
    paths = json.load(f)

PERFORMANCE_GRAPHS_PATH = paths["figs"]["performance"]


def inaccuracy_plot(df, path, size=(10, 6), fontsize=12, fontsize_amp=2):
    plt.figure(figsize=size)
    df = df[df["algorithm"] != "cpu -"]

    ax = sns.barplot(data=df, x="dataset", y="inaccuracy", hue="algorithm", palette="Set2")
    ax.set_ylim(bottom=0)
    ax.legend(title="Algorithm")
    ax.set_xlabel('Dataset', fontsize=fontsize + fontsize_amp)
    ax.set_ylabel('Inaccuracy', fontsize=fontsize + fontsize_amp)
    ax.set_title("Algorithm Comparisons: Inaccuracy Result", fontsize=fontsize + fontsize_amp * 2)
    plt.xticks(fontsize=fontsize)  # Increase x-axis label font size
    plt.yticks(fontsize=fontsize)  # Increase y-axis label font size

    plt.tight_layout()
    plt.savefig(path)
    plt.clf()


def runtime_plot(df, path, size=(10, 6), fontsize=12, fontsize_amp=2):
    plt.figure(figsize=size)

    ax = sns.barplot(data=df, x="dataset", y="total runtime", hue="algorithm", palette="Set2")
    ax.set_yscale("log")
    ax.set_xlabel('Dataset', fontsize=fontsize + fontsize_amp)
    ax.set_ylabel('Runtime (seconds)', fontsize=fontsize + fontsize_amp)
    ax.legend(title="Algorithm", fontsize=fontsize)
    ax.set_title("Algorithm Comparisons: Runtime Results", fontsize=fontsize + fontsize_amp * 2)
    plt.xticks(fontsize=fontsize)  # Increase x-axis label font size
    plt.yticks(fontsize=fontsize)  # Increase y-axis label font size

    plt.tight_layout()
    plt.savefig(path)
    plt.clf()


if __name__ == "__main__":
    files = get_files(RAW_CLUSTERS_DIR)
    df = pd.DataFrame(columns=["dataset", "algorithm", "total runtime", "inaccuracy"])
    for file_name, file in files.items():
        print(file_name)
        db_result_dir = os.path.join(PERFORMANCE_RESULTS_PATH, file_name)

        algorithms = {"gpu +": ResultFile(os.path.join(db_result_dir, "gpu_opt.csv")),
                      "gpu -": ResultFile(os.path.join(db_result_dir, "gpu_unopt.csv")),
                      "cpu -": ResultFile(os.path.join(db_result_dir, "cpu_unopt.csv"))}

        for name, result in algorithms.items():
            result.complete_df(file.raw)
            row = {"dataset": file_name, "algorithm": name, "total runtime": result.total_runtime * result.divide_by, "inaccuracy": result.inaccuracy()}
            df.loc[len(df.index)] = row
            print(name, result.divide_by)
        # debugging
        if "Full" in file_name:
            gpu_opt = algorithms["gpu +"]
            gpu = algorithms["gpu -"]
            length_opt = len(gpu_opt.df.index) // gpu_opt.divide_by
            length = len(gpu.df.index) // gpu.divide_by
            gpu_opt_index = pd.Index((gpu_opt.df.loc[: length_opt, "classification"][~gpu_opt.df.loc[: length_opt, "classification"]]).index)
            gpu_index = pd.Index((gpu.df.loc[: length, "classification"][~gpu.df.loc[: length, "classification"]]).index)
            print(gpu_index.difference(gpu_opt_index))
            print(gpu_opt_index.difference(gpu_index))

    inaccuracy_plot(df, os.path.join(PERFORMANCE_GRAPHS_PATH, "inaccuracy.png"), size=(9, 7), fontsize=10)
    runtime_plot(df, os.path.join(PERFORMANCE_GRAPHS_PATH, "runtime.png"), fontsize=11.5)
    print("done!")
