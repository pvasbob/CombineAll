from mpi4py import MPI
import os
import shutil
import subprocess

def create_working_dir():
    working_dir = "Working_dir"
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    return working_dir

def create_subprocess_folder(rank):
    folder_name = f"folder_{rank}"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    return folder_name

def copy_input_files(src_folder, dest_folder):
    input_files_folder = "./InputFileTemplate"
    for file_name in os.listdir(input_files_folder):
        src_path = os.path.join(input_files_folder, file_name)
        dest_path = os.path.join(dest_folder, file_name)
        shutil.copy2(src_path, dest_path)
    print(f"Input files copied to {dest_folder}")

def run_fortran_executable(folder_name, rank):
    executable_path = "./Exes/hfbtho_main"  # Replace with the actual path to your Fortran executable
    subprocess.run([executable_path], cwd=folder_name)
    print(f"Process {rank} finished running Fortran executable in {folder_name}")

def copy_output_files(src_folder, dest_folder):
    for file_name in os.listdir(src_folder):
        src_path = os.path.join(src_folder, file_name)
        dest_path = os.path.join(dest_folder, f"{os.path.basename(src_folder)}_{file_name}")
        shutil.copy2(src_path, dest_path)
    print(f"Output files from {src_folder} copied to {dest_folder}")

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    working_dir = create_working_dir()
    subprocess_folder = create_subprocess_folder(rank)

    copy_input_files(subprocess_folder, subprocess_folder)
    run_fortran_executable(subprocess_folder, rank)

    comm.barrier()  # Synchronize processes

    if rank == 0:
        output_folder = "Output"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        for i in range(comm.Get_size()):
            src_folder = f"folder_{i}"
            copy_output_files(src_folder, output_folder)

if __name__ == "__main__":
    main()
