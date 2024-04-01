import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from cluster_files import get_files, RAW_CLUSTERS_DIR
from results_files import ResultFile
from performance_comp_testing import PERFORMANCE_RESULTS_PATH

PERFORMANCE_GRAPHS_PATH = "/home/noam/alphaProject/figs/performance_comp"


def inaccuracy_plot(df, path):
    plt.figure(figsize=(10, 6))
    title = f"Algorithm Comparisons: Inaccuracy Result"
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(title)
    df = df[df["algorithm"] != "cpu unoptimized"]
    datasets = df["dataset"].unique()
    axes = axes.flatten()
    for i, dataset in enumerate(datasets):
        dataset_df = df[df["dataset"] == dataset]
        ax = axes[i]
        sns.barplot(data=dataset_df, x="algorithm", y="inaccuracy", ax=ax)
        for j in ax.containers:
            ax.bar_label(j, )
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Algorithm')
        ax.set_ylabel('Inaccuracy')
        ax.set_title(dataset)
    plt.tight_layout()
    plt.savefig(path)


def runtime_plot(df, path):
    plt.figure(figsize=(10, 6))
    title = f"Algorithm Comparisons: Runtime Results"
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(title)
    datasets = df["dataset"].unique()
    axes = axes.flatten()
    for i, dataset in enumerate(datasets):
        dataset_df = df[df["dataset"] == dataset]
        ax = axes[i]
        sns.barplot(data=dataset_df, x="algorithm", y="total runtime", ax=ax)
        for j in ax.containers:
            ax.bar_label(j, )
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Algorithm')
        ax.set_ylabel('Runtime (seconds)')
        ax.set_title(dataset)
    plt.tight_layout()
    plt.savefig(path)


if __name__ == "__main__":
    files = get_files(RAW_CLUSTERS_DIR)
    df = pd.DataFrame(columns=["dataset", "algorithm", "total runtime", "inaccuracy"])
    for file_name, file in files.items():
        print(file_name)
        db_result_dir = os.path.join(PERFORMANCE_RESULTS_PATH, file_name)

        algorithms = {"gpu optimized": ResultFile(os.path.join(db_result_dir, "gpu_opt.csv")),
                      "gpu unoptimized": ResultFile(os.path.join(db_result_dir, "gpu_unopt.csv")),
                      "cpu unoptimized": ResultFile(os.path.join(db_result_dir, "cpu_unopt.csv"))}

        for name, result in algorithms.items():
            result.complete_df(file.raw)
            row = {"dataset": file_name, "algorithm": name, "total runtime": result.total_runtime * result.divide_by, "inaccuracy": result.inaccuracy()}
            df.loc[len(df.index)] = row
            print(name, result.divide_by)
        if "Full" in file_name:
            gpu_opt = algorithms["gpu optimized"]
            gpu = algorithms["gpu unoptimized"]
            length_opt = len(gpu_opt.df.index) // gpu_opt.divide_by
            length = len(gpu.df.index) // gpu.divide_by
            gpu_opt_index = pd.Index((gpu_opt.df.loc[: length_opt, "classification"][~gpu_opt.df.loc[: length_opt, "classification"]]).index)
            gpu_index = pd.Index((gpu.df.loc[: length, "classification"][~gpu.df.loc[: length, "classification"]]).index)
            print(gpu_index.difference(gpu_opt_index))
            print(gpu_opt_index.difference(gpu_index))

    inaccuracy_plot(df, os.path.join(PERFORMANCE_GRAPHS_PATH, "inaccuracy.png"))
    # runtime_plot(df, os.path.join(PERFORMANCE_GRAPHS_PATH, "runtime.png"))
    print("done!")
