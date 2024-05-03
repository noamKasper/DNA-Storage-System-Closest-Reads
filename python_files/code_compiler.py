import os
import subprocess
from cluster_files import ClusterFile, RAW_CLUSTERS_DIR
import time

N_READS_GPU_OPT_PATH = '/home/noam/alphaProject/edit_distance/gpu_opt.cu'
GPU_OPT_PATH = '/home/noam/alphaProject/edit_distance/gpu_opt_bench.cu'
GPU_UNOPT_PATH = '/home/noam/alphaProject/edit_distance/gpu_unopt.cu'
CPU_UNOPT_PATH = '/home/noam/alphaProject/edit_distance/cpu_unopt.cpp'
EXECUTABLE_PATH = '/home/noam/alphaProject/edit_distance/program'


class Compiler:
    def __init__(self, read_length, compiler_type="nvcc", params=None, o3_opt=False):
        if params is None:
            params = {}
        self.__compiler_type = compiler_type
        self.params = params
        self.params["READ_LENGTH"] = read_length
        self.o3 = o3_opt

    def compile(self, output_path, code_path):
        command = f"{self.__compiler_type} {code_path} -o {output_path} "
        param_list = []
        for param, value in self.params.items():
            param_list.append(f"-D{param}={value}")
        command += " ".join(param_list)
        if self.o3:
            command += " -O3"
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


def run_program_and_save_output(padded_path, output_path, message="", exec_path=EXECUTABLE_PATH):
    start_time = time.time()
    with open(output_path, 'w') as log_file, open(padded_path, 'r') as input_file:
        subprocess.run(f'./{exec_path}', shell=True, stdout=log_file, stderr=subprocess.STDOUT, stdin=input_file)
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
        output_dir = "/home/noam/alphaProject/effectivenessTestResults/GPUSimple/notUcharOptimized"
        output_path = os.path.join(output_dir, dir_name + ".csv")
        print(f"{dir_name}({f'{i + 1}/{len(files)}'}):")
        compiler = Compiler(read_length=file.get_max_read_length(), uchar_optimization=False, divide_data_by=1000)
        compiler.compile(code_path=GPU_UNOPT_PATH, output_path=EXECUTABLE_PATH)
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
