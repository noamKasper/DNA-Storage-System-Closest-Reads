from code_compiler import Compiler, run_program_and_save_output, N_READS_GPU_OPT_PATH, EXECUTABLE_PATH
from cluster_files import ClusterFile, RAW_DATASETS, get_files
import os

with open("../paths.json", "r") as f:
    paths = json.load(f)

N_CLOSEST_RESULTS_PATH = paths["results"]["n closest"]


def n_sensitivity_testing(path_, n_range: iter = range(1, 10)):
    for n in n_range:
        file_results_path = os.path.join(path_, f"{n}.csv")
        n_testing(file_results_path, n)


def n_testing(path_, n: int = 5):
    compiler = Compiler(read_length=file.get_max_read_length(), params={"K_CLOSEST": n})
    compiler.compile(code_path=N_READS_GPU_OPT_PATH, output_path=EXECUTABLE_PATH)
    run_program_and_save_output(file.prepared, path_, exec_path=EXECUTABLE_PATH)


if __name__ == "__main__":
    files = get_files(RAW_CLUSTERS_DIR)
    for file_name, file in files.items():
        print(file_name)
        path = os.path.join(N_CLOSEST_RESULTS_PATH, file_name)
        if not os.path.exists(path):
            os.mkdir(path)
        n_sensitivity_testing(path, [1, 3, 5, 10, 20])
