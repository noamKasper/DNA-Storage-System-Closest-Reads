import os
import argparse
import subprocess
from cluster_files import ClusterFile, RAW_DATASETS
import time

with open("../paths.json", "r") as f:
    paths = json.load(f)

ALGORITHMS = paths["algorithms"]

N_READS_GPU_OPT_PATH = ALGORITHMS["n reads gpu opt"]
GPU_OPT_PATH = ALGORITHMS["gpu opt"]
GPU_UNOPT_PATH = ALGORITHMS["gpu unopt"]
CPU_UNOPT_PATH = ALGORITHMS["cpu unopt"]
EXECUTABLE_PATH = ALGORITHMS["executable"]


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
