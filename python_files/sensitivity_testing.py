from automate_reruns import Compiler, run_program_and_save_output, RAW_CLUSTERS_DIR, CODE_PATH, EXECUTABLE_PATH
from cluster_files import ClusterFile
import os

K_RESULT_PATH = "/home/noam/cuda_codes/kresults"
K_RANGE = range(9, 16)
K_DEFAULT = 12

ETH_RESULT_PATH = "/home/noam/cuda_codes/ethresults"
ETH_RANGE = range(3, 20)
ETH_DEFAULT = 10

COMBINED_RESULTS = "/home/noam/cuda_codes/combined_dir_results"
COMBINED_K_RANGE = range(10, 15, 2)
COMBINED_ETH_RANGE = range(6, 13, 2)

# in order to check that it works
DIVIDE_BY = 1


def get_files(path) -> dict:
    files_ = {}
    for file_name in os.listdir(path):
        if "evyaPfitserPsuedo.txt" == file_name:
            files_[file_name.replace(".txt", "")] = ClusterFile(os.path.join(RAW_CLUSTERS_DIR, file_name))
    return files_


def sensitivity(file: ClusterFile, result_path: str, sens_range: iter, k=False, eth=False, custom_k=None, custom_eth=None):
    if not (k or eth):
        print("sensitivity function has no objective!")
        return
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    for run in sens_range:
        output_path = os.path.join(result_path, f"{run}.csv")
        print(run)
        if k:
            if custom_eth is None:
                compiler = Compiler(read_length=file.get_max_read_length(), uchar_optimization=True, divide_data_by=DIVIDE_BY, k=run)
            else:
                compiler = Compiler(read_length=file.get_max_read_length(), uchar_optimization=True, divide_data_by=DIVIDE_BY, k=run, eth=custom_eth)
        else:
            if custom_k is None:
                compiler = Compiler(read_length=file.get_max_read_length(), uchar_optimization=True, divide_data_by=DIVIDE_BY, eth=run)
            else:
                compiler = Compiler(read_length=file.get_max_read_length(), uchar_optimization=True, divide_data_by=DIVIDE_BY, eth=run, k=custom_k)

        compiler.compile(code_path=CODE_PATH, output_path=EXECUTABLE_PATH)
        run_program_and_save_output(file.prepared, output_path)


def simple_sensitivity_testing(files_: dict, k_range, eth_range):
    for file_name, file in files_.items():
        print(file_name)
        print("k sensitivity:")
        sensitivity(file, os.path.join(K_RESULT_PATH, file_name), k_range, k=True)

        print("eth sensitivity:")
        sensitivity(file, os.path.join(ETH_RESULT_PATH, file_name), eth_range, eth=True)


def complex_sensitivity_testing(files_: dict, k_range, eth_range):
    for file_name, file in files_.items():
        print(file_name)
        database_results_dir = os.path.join(COMBINED_RESULTS, file_name)

        if not os.path.exists(database_results_dir):
            os.makedirs(database_results_dir)

        for run in k_range:
            k_result_dir = os.path.join(database_results_dir, str(run))
            if not os.path.exists(k_result_dir):
                os.makedirs(k_result_dir)

            sensitivity(file, k_result_dir, eth_range, eth=True, custom_k=run)


if __name__ == "__main__":
    files = get_files(RAW_CLUSTERS_DIR)
    print(files)
    simple_sensitivity_testing(files, K_RANGE, ETH_RANGE)
    # complex_sensitivity_testing(files, COMBINED_K_RANGE, COMBINED_ETH_RANGE)
    print("I am done!")
