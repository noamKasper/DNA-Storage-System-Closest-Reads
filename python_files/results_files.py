import os
import pandas as pd
from cluster_files import ClusterFile, RAW_CLUSTERS_DIR, PADDED_CLUSTERS_DIR, PREPARED_READS_DIR
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
        self.df = pd.read_csv(path, skiprows=1)
        self.get_settings()

    def get_settings(self):
        with open(self.path, "r") as f:
            settings = f.readline().replace("\n", "")
            if settings.startswith("settings: "):
                settings = settings.replace("settings: ", "").split(",")
                settings = tuple(map(lambda x: float(x) if x.replace(".", "").isnumeric() else x, settings))
                if len(settings) == 7:
                    self.eth, self.k, self.duped_reads_buffer, self.mean_thread_runtime, self.mean_thread_edit_calcs, self.total_runtime, self.message = settings
                elif len(settings) == 6:
                    self.eth, self.k, self.duped_reads_buffer, self.mean_thread_runtime, self.mean_thread_edit_calcs, self.total_runtime = settings
                else:
                    throw("Error while reading file settings")

    def complete_df(self, raw_file_path):

        if not os.path.exists(raw_file_path):
            print("Path", raw_file_path, "does not exist")
            return

        with open(raw_file_path, "r") as f:
            strand_read_map = []
            for strandIdx, cluster in enumerate(f.read().strip().split("\n\n\n")):
                for _ in cluster.split("\n*****************************\n")[1].split():
                    strand_read_map.append(strandIdx)

        self.df["strand"] = pd.Series(strand_read_map)

        mask = self.df['closest read'] != -1

        closest_read_indices = self.df.loc[mask, 'closest read'].values
        closest_read_strands = self.df.loc[closest_read_indices, 'strand'].values

        classifications = pd.Series(False, index=self.df.index)  # Initialize with False for all rows
        classifications[mask] = self.df.loc[mask, 'strand'].values == closest_read_strands

        self.df["classification"] = classifications
