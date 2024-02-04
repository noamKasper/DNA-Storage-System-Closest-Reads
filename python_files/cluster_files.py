import os

RAW_CLUSTERS_DIR = "/home/noam/alphaProject/reads/raw_files"
PADDED_CLUSTERS_DIR = "/home/noam/alphaProject/reads/padded_clusters"
PREPARED_READS_DIR = "/home/noam/alphaProject/reads/prepared_reads"


def get_files(path) -> dict:
    files_ = {}
    for file_name in os.listdir(path):
        files_[file_name.replace(".txt", "")] = ClusterFile(os.path.join(RAW_CLUSTERS_DIR, file_name))
    return files_


class ClusterFile:
    def __init__(self, path):
        for i in [RAW_CLUSTERS_DIR, PADDED_CLUSTERS_DIR, PREPARED_READS_DIR]:
            if not os.path.exists(i):
                os.makedirs(i)

        if os.path.basename(path) in os.listdir(RAW_CLUSTERS_DIR):
            self.name: str = os.path.basename(path)
            self.raw = os.path.join(RAW_CLUSTERS_DIR, self.name)
            print(self.name, self.raw)
        else:
            print(f"{os.path.basename(path)} wasn't found in {RAW_CLUSTERS_DIR}!")
        self.padded = os.path.join(PADDED_CLUSTERS_DIR, self.name)
        self.prepared = os.path.join(PREPARED_READS_DIR, self.name)

        self.pad_clusters(silent=True)
        self.prepare_clusters(silent=True)

    def get_clusters(self):
        with open(self.raw, "r") as f:
            clusters = {i.strip().split("\n*****************************\n")[0]: [j for j in i.strip().split(
                "\n*****************************\n")[1].split("\n")] for i in f.read().strip().split("\n\n\n")}
            return clusters

    def pad_clusters(self, force: bool = False, silent: bool = False):
        if self.name in os.listdir(PADDED_CLUSTERS_DIR) and not force:
            if not silent:
                print(
                    f"{self.name} is already in {PADDED_CLUSTERS_DIR} if you want to force set force attribute to true.")
            return

        with open(self.raw, 'r') as f:
            lines = f.read().strip().split('\n')
            max_length = len(max(lines, key=len))

            for i in range(len(lines)):
                if lines[i] not in ['', '*****************************']:
                    lines[i] = lines[i] + 'P' * (max_length - len(lines[i]))

        with open(self.padded, 'w') as f:
            f.write('\n'.join(lines).strip())

    def prepare_clusters(self, force: bool = False, silent: bool = False):
        if self.name in os.listdir(PREPARED_READS_DIR) and not force:
            if not silent:
                print(
                    f"{self.name} is already in {PREPARED_CLUSTERS_DIR} if you want to force set force attribute to true.")
            return

        with open(self.padded, "r") as f:
            mutations = []
            for i in f.read().split("\n\n\n"):
                for j in i.split("\n*****************************\n")[1].split("\n"):
                    mutations.append(j)

        with open(self.prepared, "w") as f:
            f.write(''.join(mutations))

    def get_mean_cluster_length(self):
        with open(self.raw, "r") as f:
            clusters = [len(i.split("\n*****************************\n")[1].split("\n")) for i in f.read().strip().split("\n\n\n")]
            return sum(clusters) / len(clusters)

    def get_max_read_length(self):
        with open(self.raw, "r") as f:
            lines = f.read().split('\n')
            return len(max(lines, key=len))

    def get_file_size(self):
        if not os.path.exists(self.prepared):
            self.prepare_clusters()
        return os.path.getsize(self.prepared)

    def get_num_reads(self):
        with open(self.raw, "r") as f:
            reads = [j for i in f.read().strip().split("\n\n\n") for j in i.split("\n*****************************\n")[1].split("\n")]
            return len(reads)


if __name__ == "__main__":
    for i in os.listdir(RAW_CLUSTERS_DIR):
        cluster = ClusterFile(os.path.join(RAW_CLUSTERS_DIR, i))
        print(i)
        print("number of reads:", cluster.get_num_reads())
        print("mean cluster length:", cluster.get_mean_cluster_length())
        print("max read length:", cluster.get_max_read_length())
        print("file size: ", cluster.get_file_size() / (1024 * 1024), "MB")
