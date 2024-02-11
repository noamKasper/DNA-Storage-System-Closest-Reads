import os

import pandas as pd
from results_files import NResultsFile
from n_closest_testing import N_CLOSEST_RESULTS_PATH
from cluster_files import ClusterFile

POST_PROCESSING_PATH = "/home/noam/alphaProject/results/n_closest_results/post_processing"
POST_PROCESSING_FULL_RESULTS_PATH = os.path.join(POST_PROCESSING_PATH, "full_results.csv")


if __name__ == "__main__":
    results_df = pd.DataFrame(columns=["ST", "SF", "LT", "LF"])
    for file_name in os.listdir(N_CLOSEST_RESULTS_PATH):
        print(file_name)
        if ".csv" not in file_name:
            continue
        cluster_file = ClusterFile(file_name.replace(".csv", ".txt"))
        result_path = os.path.join(N_CLOSEST_RESULTS_PATH, file_name)
        post_path = os.path.join(POST_PROCESSING_PATH, file_name)
        result_file = NResultsFile(result_path)
        if not os.path.exists(post_path):
            result_file.complete_df(cluster_file.raw)
            result_file.completed_df.to_csv(post_path, index=True, index_label="read")
        else:
            print("file exists, skipping...")
            result_file.completed_df = pd.read_csv(post_path)
        print("calculating category...")
        results_df.loc[file_name.replace(".csv", "")] = result_file.get_results()
        print("file is calculated!")
    results_df = results_df.fillna(0)
    results_df.to_csv(POST_PROCESSING_FULL_RESULTS_PATH, index=True, index_label="dataset")
    print("done!")