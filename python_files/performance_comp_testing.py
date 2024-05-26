from code_compiler import Compiler, run_program_and_save_output, GPU_OPT_PATH, GPU_UNOPT_PATH, CPU_UNOPT_PATH, EXECUTABLE_PATH
from cluster_files import ClusterFile, RAW_DATASETS, get_files
import os

with open("../paths.json", "r") as f:
    paths = json.load(f)

K = 12
ETH = 10
DIVIDE_BY = 1000
PERFORMANCE_RESULTS_PATH = paths["results"]["performance"]
RUNNING_TIME = 900  # in seconds


def do_performance_comparison(file: ClusterFile, results_dir_path: str, divide_by=DIVIDE_BY, gpu_opt=True, gpu=True, cpu=True):
    params = {"DIVIDE_DATA_BY":divide_by}
    cuda_opt_compiler = Compiler(read_length=file.get_max_read_length(), params=params)
    cuda_compiler = Compiler(read_length=file.get_max_read_length(), params=params)
    cpp_compiler = Compiler(read_length=file.get_max_read_length(), params=params, compiler_type="g++", o3_opt=True)

    if gpu_opt:
        print("gpu+")
        cuda_opt_compiler.compile(code_path=GPU_OPT_PATH, output_path=EXECUTABLE_PATH)
        run_program_and_save_output(file.prepared, os.path.join(results_dir_path, "gpu_opt.csv"))

    if gpu:
        print("gpu")
        cuda_compiler.compile(code_path=GPU_UNOPT_PATH, output_path=EXECUTABLE_PATH)
        run_program_and_save_output(file.prepared, os.path.join(results_dir_path, "gpu_unopt.csv"))

    if cpu:
        print("cpu")
        cpp_compiler.compile(code_path=CPU_UNOPT_PATH, output_path=EXECUTABLE_PATH)
        run_program_and_save_output(file.prepared, os.path.join(results_dir_path, "cpu_unopt.csv"))


if __name__ == "__main__":
    files = get_files(RAW_CLUSTERS_DIR)
    for file_name, file in files.items():
        print(file_name)
        file_results_path = os.path.join(PERFORMANCE_RESULTS_PATH, file_name)

        if not os.path.exists(file_results_path):
            os.makedirs(file_results_path)
        read_length = file.get_max_read_length()
        num_reads = file.get_num_reads()
        divide_by = int(((read_length**2) * (num_reads**2)) / ((10**9) * (RUNNING_TIME/3)))
        print(divide_by)
        do_performance_comparison(file, file_results_path, divide_by=100, gpu_opt=True, gpu=True, cpu=False)
    print("done!")