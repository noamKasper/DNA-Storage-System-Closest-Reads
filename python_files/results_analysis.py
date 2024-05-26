import pandas as pd
from results_files import ResultFile, RAW_DATASETS
from cluster_files import ClusterFile
import seaborn as sns
import matplotlib.pyplot as plt
import os


def saveDataSetsCSV():
    # results_path = "/home/noam/cuda_codes/Results"
    results_path = "/home/noam/alphaProject/effectivenessTestResults/GPUSimple"
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
