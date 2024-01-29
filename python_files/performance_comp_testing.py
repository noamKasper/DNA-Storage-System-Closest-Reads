from automate_reruns import Compiler, run_program_and_save_output, GPU_OPT_BENCH_PATH, GPU_UNOPT_PATH, CPU_UNOPT_PATH, EXECUTABLE_PATH
from cluster_files import ClusterFile, RAW_CLUSTERS_DIR, get_files
import os

K = 12
ETH = 10
DIVIDE_BY = 1000
PERFORMANCE_RESULTS_PATH = '/home/noam/alphaProject/results/performance_comparison_results'


def do_performance_comparison(file: ClusterFile, results_dir_path: str):
    cuda_compiler = Compiler(read_length=file.get_max_read_length(), divide_data_by=DIVIDE_BY)
    cpp_compiler = Compiler(read_length=file.get_max_read_length(), divide_data_by=DIVIDE_BY, compiler_type="g++")

    print("gpu+")
    cuda_compiler.compile(code_path=GPU_OPT_BENCH_PATH, output_path=EXECUTABLE_PATH)
    run_program_and_save_output(file.prepared, os.path.join(results_dir_path, "gpu_opt.csv"))

    print("gpu")
    cuda_compiler.compile(code_path=GPU_UNOPT_PATH, output_path=EXECUTABLE_PATH)
    run_program_and_save_output(file.prepared, os.path.join(results_dir_path, "gpu_unopt.csv"))

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

        do_performance_comparison(file, file_results_path)
