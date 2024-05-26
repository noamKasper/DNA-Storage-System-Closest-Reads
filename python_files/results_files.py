import os
import pandas as pd
from cluster_files import ClusterFile, RAW_DATASETS, PADDED_DATASETS, PREPARED_DATASETS
import seaborn as sns


class ResultFile:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.eth = None
        self.k = None
        self.duped_reads_buffer = None
        self.mean_thread_runtime = None
        self.mean_thread_edit_calcs = None
        self.total_runtime = None
        self.message = None
        self.divide_by = 1
        self.df = pd.read_csv(path, skiprows=1)
        self.get_settings()

    def get_settings(self):
        with open(self.path, "r") as f:
            settings = f.readline().replace("\n", "")
            if settings.startswith("settings: "):
                settings = settings.replace("settings: ", "").split(",")
                settings = tuple(map(lambda x: float(x) if x.replace(".", "").isnumeric() else x, settings))
                if len(settings) == 8:
                    self.eth, self.k, self.divide_by, self.duped_reads_buffer, self.mean_thread_runtime, self.mean_thread_edit_calcs, self.total_runtime, self.message = settings
                elif len(settings) == 7:
                    self.eth, self.k, self.duped_reads_buffer, self.mean_thread_runtime, self.mean_thread_edit_calcs, self.total_runtime, self.message = settings
                elif len(settings) == 6:
                    self.eth, self.k, self.duped_reads_buffer, self.mean_thread_runtime, self.mean_thread_edit_calcs, self.total_runtime = settings
                elif len(settings) == 3:
                    self.divide_by, self.total_runtime, self.message = settings
                else:
                    throw("Error while reading file settings")

    def complete_df(self, raw_file_path):

        if not os.path.exists(raw_file_path):
            print("Path", raw_file_path, "does not exist")
            return

        cluster_file = ClusterFile(raw_file_path)
        with open(cluster_file.padded, "r") as f:
            strand_read_map = []
            for strandIdx, cluster in enumerate(f.read().strip().split("\n\n\n")):
                for _ in cluster.split("\n*****************************\n")[1].split():
                    strand_read_map.append(strandIdx)

        strands = pd.Series(strand_read_map)
        self.df = self.df.reindex(strands.index, fill_value=-1)
        self.df["strand"] = strands
        self.df.fillna(-1)
        mask = self.df['closest read'] != -1

        closest_read_indices = self.df.loc[mask, 'closest read'].values
        closest_read_strands = strands[closest_read_indices].values

        classifications = pd.Series(False, index=self.df.index)  # Initialize with False for all rows
        classifications[mask] = strands[mask].values == closest_read_strands

        self.df["classification"] = classifications

    def inaccuracy(self):
        if "classification" not in self.df.columns:
            print("complete df first!")
            return -1

        length = len(self.df.index) // self.divide_by
        inaccuracy = 1 - self.df.loc[: length, "classification"].sum() / length
        return inaccuracy


class NResultsFile:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.n = None
        self.eth = None
        self.k = None
        self.total_runtime = None
        self.message = None
        self.cluster_file = None
        self.df: pd.DataFrame = pd.read_csv(path, skiprows=1)
        self.completed_df: pd.DataFrame | None = None
        self.get_settings()

    def get_settings(self):
        with open(self.path, "r") as f:
            settings = f.readline().replace("\n", "")
            if settings.startswith("settings: "):
                settings = settings.replace("settings: ", "").split(",")
                settings = tuple(map(lambda x: float(x) if x.replace(".", "").isnumeric() else x, settings))
                if len(settings) == 5:
                    self.n, self.eth, self.k, self.total_runtime, self.message = settings
                    self.n = int(self.n)
                elif len(settings) == 4:
                    self.n, self.eth, self.k, self.total_runtime = settings
                    self.n = int(self.n)
                else:
                    throw("Error while reading file settings")

    def complete_df(self, raw_file_path):
        if not os.path.basename(raw_file_path) in os.listdir(RAW_CLUSTERS_DIR):
            print("Path", raw_file_path, "does not exist")
            return

        def tolist(x: str):
            x = x.strip()
            x = x.removeprefix("[ ")
            x = x.removesuffix(" ]")
            x = x.split(" ")
            x = list(map(int, x))
            return x

        def sort_row(x):
            indices = x["closest read"]
            values = x["edit distance"]
            combined = list(zip(indices, values))
            sorted_combined = sorted(combined, key=lambda x: x[1])
            x["closest read"], x["edit distance"] = zip(*sorted_combined)
            return x

        self.completed_df = self.df

        # Convert string representations to lists
        self.completed_df['closest read'] = self.completed_df['closest read'].apply(tolist)
        self.completed_df['edit distance'] = self.completed_df['edit distance'].apply(tolist)

        # Sort DataFrame based on the first edit distance column in descending order
        self.completed_df = self.completed_df.apply(sort_row, axis=1)

        # Extract the desired columns
        closest_read_cols = self.completed_df['closest read'].apply(pd.Series).iloc[:, :self.n]
        edit_distance_cols = self.completed_df['edit distance'].apply(pd.Series).iloc[:, :self.n]

        # Rename columns for clarity
        closest_read_cols.columns = [f'closest read{i + 1}' for i in range(self.n)]
        edit_distance_cols.columns = [f'edit distance{i + 1}' for i in range(self.n)]

        # Concatenate the result
        self.completed_df = pd.concat([self.completed_df[[]], closest_read_cols, edit_distance_cols], axis=1)
        self.cluster_file = ClusterFile(raw_file_path)
        with open(self.cluster_file.raw, "r") as f:
            strand_read_map = []
            cluster_length_map = []
            for strandIdx, cluster in enumerate(f.read().strip().split("\n\n\n")):
                reads = cluster.split("\n*****************************\n")[1].split()
                for _ in reads:
                    strand_read_map.append(strandIdx)
                    cluster_length_map.append(len(reads))

        strands = pd.Series(strand_read_map)
        cluster_lengths = pd.Series(cluster_length_map)

        for i in range(1, self.n + 1):
            closest_read_col = f'closest read{i}'

            mask = self.completed_df[closest_read_col] != -1

            closest_read_indices = self.completed_df.loc[mask, closest_read_col].values
            closest_read_strands = strands[closest_read_indices].values

            classifications = pd.Series(False, index=self.completed_df.index)  # Initialize with False for all rows
            classifications[mask] = strands[mask].values == closest_read_strands

            self.completed_df[f'classification{i}'] = classifications

        self.completed_df["strand"] = strands
        self.completed_df["cluster length"] = cluster_lengths

    def get_results(self) -> dict:
        """
        short clusters - reads for which cluster length < n+1
        long clusters - reads for which cluster length >= n+1


        number of reads that are in short clusters and:
        ST(short true) - all the reads within the cluster are found
        SF(short false) - not all the reads within the cluster are found


        number of reads that are in long numbers and:
        LT(long true) - all the reads that were closest are in cluster
        LF(long false) - not all the reads that were closest are in cluster

        :return: ST ,SF, LT, LF
        """
        if self.completed_df is None:
            print("df completed does not exist, get help!")
            return {}

        def categorize(x: pd.Series):
            classifications = x.filter(regex='classification\d+')
            cluster_length = x["cluster length"]
            if cluster_length < self.n + 1:
                if classifications.sum() == cluster_length - 1:
                    return "ST"
                else:
                    return "SF"
            else:
                if classifications.sum() == self.n:
                    return "LT"
                else:
                    return "LF"

        self.completed_df["result type"] = self.completed_df.apply(categorize, axis=1)
        counts = self.completed_df["result type"].value_counts().to_dict()
        return counts


if __name__ == "__main__":
    result = NResultsFile("/home/noam/alphaProject/results.csv")
    result.complete_df("evyaPfitserPsuedo.txt")
    result.completed_df.to_csv("/home/noam/alphaProject/post_results.csv", index=True, index_label="read")
    print(result.get_results())
