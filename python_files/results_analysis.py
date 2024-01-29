import pandas as pd
from results_files import ResultFile, RAW_CLUSTERS_DIR
from cluster_files import ClusterFile
import seaborn as sns
import matplotlib.pyplot as plt
import os

RESULTS_DIR = "/home/noam/cuda_codes/results2"


def create_accuracy_line_plot(s):
    sns.lineplot(data=s.rename("Accuracy"), marker='o')
    for i in range(len(s.keys())):
        plt.plot([s.keys()[i], s.keys()[i]], [0, s.values[i]], color='gray', linestyle='--', linewidth=1)
    plt.title("Accuracy VS. ETH")
    plt.tight_layout()
    plt.savefig("/home/noam/cuda_codes/figs/accuracy_ETH.png")


def create_error_rate_line_plot(s: pd.Series):
    error_rate = 1 - s.rename("Error Rate")
    sns.lineplot(data=error_rate, marker='o')
    # for i in range(len(s.keys())):
    #     plt.plot([s.keys()[i], s.keys()[i]], [0, s.values[i]], color='gray', linestyle='--', linewidth=1)
    plt.yscale('log')
    plt.title("Error Rate VS. ETH")
    plt.tight_layout()
    plt.savefig("/home/noam/cuda_codes/figs/accuracy_ETH.png")
    error_rate.to_csv("error rate vs ETH.csv")


def create_ETH_time_line_plot():
    df = pd.read_csv("/home/noam/cuda_codes/python_files/time-ETH.csv")
    sns.lineplot(data=df, x="ETH", y="Time(seconds)", marker='o')
    # plt.yscale('log')
    plt.title("Time(seconds) VS. ETH")
    plt.tight_layout()
    plt.savefig("/home/noam/cuda_codes/figs/time_ETH.png")


def create_ETH_time_error_rate_plot(s):
    error_rate = 1 - s.rename("Error Rate")
    df = pd.read_csv("/home/noam/cuda_codes/python_files/time-ETH.csv")

    ax1 = sns.lineplot(data=error_rate, marker="o", label='Error Rate')
    ax1.set_yscale("log")

    ax2 = ax1.twinx()
    sns.lineplot(data=df, x='ETH', y='Time(seconds)', marker="o", color='orange', label='Time', ax=ax2)

    horizontal_line_value = 1 - calcTrueClassification()
    ax1.axhline(y=horizontal_line_value, color='gray', linestyle='--', label='Absolute Error')

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='center right')

    ax1.get_legend().remove()
    plt.title("Error Rate VS. ETH VS. Time")
    plt.tight_layout()
    plt.savefig("/home/noam/cuda_codes/figs/accuracy_ETH.png")


def get_accuracy(classifications):
    classifications = classifications
    return classifications.sum() / len(classifications) * 100


def get_results_files(n):
    results_files = {}
    for file_name in os.listdir(RESULTS_DIR):
        results_files[file_name] = []
        for i in range(1, n + 1):
            results_files[file_name].append(ResultFile(os.path.join(RESULTS_DIR, file_name, f"{i}.csv")))
    return results_files


# results = get_results_files(24)
# time_results = pd.DataFrame(columns=list(results.keys()))
# for k, v in results.items():
#     time_results[k] = pd.Series({i.duped_reads_buffer: i.mean_thread_runtime for i in v})
#
#
# for column in time_results:
#     sns.lineplot(x=time_results.index, y=time_results[column], label=column)
#
# plt.xlabel("Read Buffer Size")
# plt.ylabel("Mean Thread Runtime(Seconds)")
# plt.yscale("log")
# plt.savefig("/home/noam/cuda_codes/figs/Mean Thread Runtime-Read Buffer Size.png")

def saveDataSetsCSV():
    # results_path = "/home/noam/cuda_codes/Results"
    results_path = "/home/noam/cuda_codes/effectivenessTestResults/GPUSimple"
    datasets = {}
    results = {}
    columns = ["File Size(MB)", "Mean Cluster Length", "Number of Clusters", "Number of Reads", "Max Read Length", "Mean Thread Edit-Distance Calculations", "Mean thread runtime(seconds)", "Total runtime", "Accuracy"]

    for i in os.listdir(RAW_CLUSTERS_DIR):
        result_path = os.path.join(results_path, i.replace(".txt", ".csv"))
        if os.path.exists(result_path) and os.path.getsize(result_path) != 0:
            datasets[i.replace(".txt", "")] = ClusterFile(os.path.join(RAW_CLUSTERS_DIR, i))
            results[i.replace(".txt", "")] = ResultFile(result_path)
            # df = results[i.replace(".txt", "")].df
            # df = df.head(len(df.index)//20)
            # results[i.replace(".txt", "")].df = df
    for i in os.listdir(RAW_CLUSTERS_DIR):
        result_path = os.path.join(results_path, "notUcharOptimized", i.replace(".txt", ".csv"))
        if os.path.exists(result_path) and os.path.getsize(result_path) != 0:
            datasets["notUcharOptimized" + i.replace(".txt", "")] = ClusterFile(os.path.join(RAW_CLUSTERS_DIR, i))
            results["notUcharOptimized" + i.replace(".txt", "")] = ResultFile(result_path)
    indexes = list(results.keys())
    df = pd.DataFrame(columns=columns, index=indexes)
    df.index.name = "Dataset"
    for i in df.index.tolist():
        dataset = datasets[i]
        result = results[i]
        df.loc[i, "File Size(MB)"] = dataset.get_file_size() / (1024 * 1024)
        df.loc[i, "Mean Cluster Length"] = dataset.get_mean_cluster_length()
        df.loc[i, "Number of Clusters"] = len(dataset.get_clusters())
        df.loc[i, "Number of Reads"] = dataset.get_num_reads()
        df.loc[i, "Max Read Length"] = dataset.get_max_read_length()
        result.complete_df(dataset.raw)
        df.loc[i, "Mean Thread Edit-Distance Calculations"] = result.mean_thread_edit_calcs
        df.loc[i, "Mean thread runtime(seconds)"] = result.mean_thread_runtime
        total_runtime = result.total_runtime
        # df.loc[i, "Total runtime"] = f"{int(total_runtime / 3600)}:{int(total_runtime % 3600 / 60)}:{round(total_runtime % 3600 % 60,2)}"
        df.loc[i, "Total runtime"] = total_runtime
        df.loc[i, "Accuracy"] = get_accuracy(result.df["classification"])
    df.to_csv("/home/noam/cuda_codes/Datasets.csv")


def verify_correctness():
    # base case - GPU uchar optimization = false, research optimization = false:
    base_case_path = "/home/noam/cuda_codes/effectivenessTestResults/GPUSimple/notUcharOptimized"

    # 2nd case - GPU uchar optimization = false, research optimization = true:
    second_case_path = "/home/noam/cuda_codes/effectivenessTestResults/GPUSimple"

    # 3rd case - GPU uchar optimization = true research optimization = false:
    third_case_path = "/home/noam/cuda_codes/Results"

    # this comparison is for the first 0.1% of each file and is not fully truthful
    base_case_results = {}
    results = {}
    for i in os.listdir(RAW_CLUSTERS_DIR):
        base_result_path = os.path.join(base_case_path, i.replace(".txt", ".csv"))
        second_result_path = os.path.join(second_case_path, i.replace(".txt", ".csv"))
        third_result_path = os.path.join(third_case_path, i.replace(".txt", ".csv"))
        base_case_results[i.replace(".txt", "")] = ResultFile(base_result_path)
        base_case_results[i.replace(".txt", "")].complete_df(os.path.join(RAW_CLUSTERS_DIR, i))
        results["uchar: " + i.replace(".txt", "")] = ResultFile(second_result_path)
        results["uchar: " + i.replace(".txt", "")].complete_df(os.path.join(RAW_CLUSTERS_DIR, i))
        results["research: " + i.replace(".txt", "")] = ResultFile(third_result_path)
        results["research: " + i.replace(".txt", "")].complete_df(os.path.join(RAW_CLUSTERS_DIR, i))
        
    df = pd.DataFrame(columns=["dataset", "optimization", "cluster correctness", "read correctness", "edit distance correctness"])
    for i, (result_name, result) in enumerate(results.items()):
        optimization, database_name = tuple(result_name.split(": "))
        df.loc[i, "dataset"] = database_name
        df.loc[i, "optimization"] = optimization
        base_case = base_case_results[database_name]
        biggest_index = len(result.df.index) // 1000
        clusters = result.df["strand"]
        base_case_clusters = base_case.df["strand"]
        df.loc[i, "cluster correctness"] = get_accuracy((clusters == base_case_clusters).iloc[0:biggest_index])
        closest_read = result.df["closest read"]
        base_case_closest_read = base_case.df["closest read"]
        df.loc[i, "read correctness"] = get_accuracy((closest_read == base_case_closest_read).iloc[0:biggest_index])
        edit_distances = result.df["edit distance"]
        base_case_edit_distances = base_case.df["edit distance"]
        edit_correctness = (edit_distances == base_case_edit_distances).iloc[0:biggest_index]
        df.loc[i, "edit distance correctness"] = get_accuracy(edit_correctness)
        # print(i, end=":\n")
        # print(biggest_index, end=", ")
        false_correctness_indices = edit_correctness.iloc[0:biggest_index][~edit_correctness].index
        # print(false_correctness_indices)
        print(result_name)
        for j in false_correctness_indices:
            print(f"index: {j}")
            print(result.df.loc[j])
            print(base_case.df.loc[j])
            print("------------------------------\n")
    df.to_csv("/home/noam/cuda_codes/correctness.csv", index=False)
# saveDataSetsCSV()
# verify_correctness()

