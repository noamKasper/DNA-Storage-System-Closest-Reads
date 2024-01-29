import os
import subprocess
from cluster_files import ClusterFile
import time

RAW_CLUSTERS_DIR = "/home/noam/cuda_codes/reads/raw_files"
CODE_PATH = '/home/noam/cuda_codes/copy_edit_dis.cu'
EXECUTABLE_PATH = '/home/noam/cuda_codes/program'
PREPARED_DIR = '/home/noam/cuda_codes/reads/prepared_reads'
SIZES_JSON = '/home/noam/cuda_codes/reads/size.json'
RESULTS_DIR = '/home/noam/cuda_codes/Results2'
RESULTS_DIR1 = '/home/noam/cuda_codes/kresults'
RESULTS_DIR2 = '/home/noam/cuda_codes/ethresults'


class Compiler:
    def __init__(self, read_length=115, compiler_type="nvcc", k=None, eth=None, used_reads_buffer_size=None, uchar_optimization: bool = True, force_uchar=False, divide_data_by: int = 1):
        self._compiler_type = compiler_type
        self._parameters = {"READ_LENGTH": f"READ_LENGTH={read_length}",
                            "KMER": f"KMER={k}" if (k is not None) else "",
                            "ETH": f"ETH={eth}" if (eth is not None) else "",
                            "USED_READS_SIZE": f"USED_READS_SIZE={used_reads_buffer_size}" if (used_reads_buffer_size is not None) else "",
                            "UCHAR4_OPTIMIZATION": f"UCHAR4_OPTIMIZATION={str(uchar_optimization).lower()}",
                            "FORCE_UCHAR4_OPTIMIZATION": f"FORCE_UCHAR4_OPTIMIZATION={str(force_uchar).lower()}",
                            "DIVIDE_DATA_BY": f"DIVIDE_DATA_BY={divide_data_by}"}

    def compile(self, output_path, code_path):
        command = f"{self._compiler_type} {code_path} -o {output_path}"
        for i in self._parameters.values():
            if i:
                command += " -D" + i
        print(command)
        result = subprocess.run(command, shell=True)

        if result.returncode != 0:
            print(f"Error: Compilation failed with return code {result.returncode}")
            if result.stderr:
                print("Error message:")
                print(result.stderr.decode('utf-8'))
            exit(result.returncode)


# def compile_program(code_path, output_filename, N, compiler="nvcc", K=None, ETH=None, USED_READS_SIZE=None, not_uchar_optimized=None):
#     compile_command = f'{compiler} -o {output_filename} {code_path} -DREAD_LENGTH={N} {f"-DKMER={K}" if (K is not None) else ""} {f"-DETH={ETH}" if (ETH is not None) else ""} {f"-DUSED_READS_SIZE={USED_READS_SIZE}" if (USED_READS_SIZE is not None) else ""} {f"-DNON_UCHAR4_OPTIMIZED=true" if (not_uchar_optimized is not None) else ""}'
#     result = subprocess.run(compile_command, shell=True)
#
#     if result.returncode != 0:
#         print(f"Error: Compilation failed with return code {result.returncode}")
#         if result.stderr:
#             print("Error message:")
#             print(result.stderr.decode('utf-8'))
#             exit(result.returncode)


def run_program_and_save_output(padded_path, output_path, message=""):
    print(output_path)
    start_time = time.time()
    with open(output_path, 'w') as log_file, open(padded_path, 'r') as input_file:
        subprocess.run(f'./{EXECUTABLE_PATH}', shell=True, stdout=log_file, stderr=subprocess.STDOUT, stdin=input_file)
    message = "," + message
    with open(output_path, "r+") as f:
        lines = f.readlines()
        lines[0] = lines[0][:-1] + f",{time.time() - start_time}{message}\n"
        f.seek(0)
        f.writelines(lines)


if __name__ == "__main__":

    files = {}
    for file_name in os.listdir(PREPARED_DIR):
        files[file_name.replace(".txt", "")] = ClusterFile(os.path.join(RAW_CLUSTERS_DIR, file_name))

    for i, (dir_name, file) in enumerate(files.items()):
        output_dir = "/home/noam/cuda_codes/effectivenessTestResults/GPUSimple/notUcharOptimized"
        output_path = os.path.join(output_dir, dir_name + ".csv")
        print(f"{dir_name}({f'{i + 1}/{len(files)}'}):")
        compiler = Compiler(read_length=file.get_max_read_length(), uchar_optimization=False, divide_data_by=1000)
        compiler.compile(code_path="/home/noam/cuda_codes/trueResults.cu", output_path=EXECUTABLE_PATH)
        # compile_program("/home/noam/cuda_codes/trueResults.cu", EXECUTABLE_PATH, file.get_max_read_length(), not_uchar_optimized=True)
        run_program_and_save_output(file.prepared, output_path)

    # for i, (dir_name, file) in enumerate(files.items()):
    #     output_dir = "/home/noam/cuda_codes/effectivenessTestResults/GPUSimple"
    #     output_path = os.path.join(output_dir, dir_name + ".csv")
    #     print(f"{dir_name}({f'{i + 1}/{len(files)}'}):")
    #     compiler = Compiler(read_length=file.get_max_read_length(), force_uchar=True, divide_data_by=1000)
    #     compiler.compile(code_path="/home/noam/cuda_codes/trueResults.cu", output_path=EXECUTABLE_PATH)
    #     # compile_program("/home/noam/cuda_codes/trueResults.cu", EXECUTABLE_PATH, file.get_max_read_length(), not_uchar_optimized=True)
    #     run_program_and_save_output(file.prepared, output_path)
    print("I am over!")
    # for i, (dir_name, file) in enumerate(sorted(files.items(), key=lambda x: x[0], reverse=False)):
    #     output_dir = "/home/noam/cuda_codes/effectivenessTestResults/CPUSimple"
    #     output_path = os.path.join(output_dir, dir_name.replace(".txt", ".csv"))
    #     print(f"{dir_name}({f'{i + 1}/{len(files)}'}):")
    #     compile_program("/home/noam/cuda_codes/get_true_results.cpp", EXECUTABLE_PATH, file.get_max_read_length(), compiler="g++")
    #     run_program_and_save_output(file.prepared, output_path)

    # for i, (dir_name, file) in enumerate(files.items()):
    #     output_dir = os.path.join(RESULTS_DIR, dir_name)
    #     # output_file = output_dir+".csv"
    #     print(f"{dir_name}({f'{i + 1}/{len(os.listdir(PREPARED_DIR))}'}):")
    #     if not os.path.exists(output_dir):
    #         os.makedirs(output_dir)
    #     sensitivity_test(12, 12, file, output_dir)
