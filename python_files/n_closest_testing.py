from automate_reruns import Compiler, run_program_and_save_output, GPU_OPT_PATH, EXECUTABLE_PATH
from cluster_files import ClusterFile, RAW_CLUSTERS_DIR, get_files
import os

N_CLOSEST_RESULTS_PATH = "/home/noam/alphaProject/results/n_closest_results"


if __name__ == "__main__":
    files = get_files(RAW_CLUSTERS_DIR)
    for file_name, file in files.items():
        print(file_name)
        file_results_path = os.path.join(N_CLOSEST_RESULTS_PATH, file_name + ".csv")
        compiler = Compiler(read_length=file.get_max_read_length())
        compiler.compile(code_path=GPU_OPT_PATH,  output_path=EXECUTABLE_PATH)
        run_program_and_save_output(file.prepared, file_results_path, exec_path=EXECUTABLE_PATH)
