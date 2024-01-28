from mpi4py import MPI
import os
import sys
import shutil
import subprocess


from Input.hfbtho_NAMELIST import hfbtho_NAMELIST
from Input.qrpa_inp import qrpa_input

lpfam_dir = {
   "working" : "Working",
   "output" : "Output",
   "Exes" : "Exes",
   "hfbtho_allpart" : "hfbtho_allpart",
   "hfbtho_diffpart" : "hfbtho_diffpart",
   "InputTemplate" : "InputTemplate",
}

lpInput_files = {
    "hfbtho_namelist" : "hfbtho_NAMELIST.dat",
    "qrpa.inp": "qrpa.inp",
}

rbm_dir = {
    "trainingData_dir": "TrainingData",
    "H11ToPPnfam:" : "H11ToPPnfam",
}


def mkdirWorking(top, label):
    working_dir = os.path.join(top, label)
    if os.path.exists(working_dir):
        print(f"file {label} exists, check and delete it manually. ")
        sys.exit()
    else:
        os.mkdir(working_dir)
        return working_dir
        

def mkdirJobs(top, label):
    job_dir = os.path.join(top, label.zfill(4))
    os.mkdir(job_dir)
    return job_dir

def copyExeToJob(src_label, src_top, dest_dir):
    shutil.copy2(os.path.join(src_top, src_label), os.path.join(dest_dir))

def copyInputToJob(src_label, src_top, dest_dir):
    if not os.path.exists(os.path.join(src_top, src_label)):
        print(f"error, no input file {src_label}")
        sys.exit()
    shutil.copy2(os.path.join(src_top, src_label), os.path.join(dest_dir))


def run_lpfam(job_dir, exe, rank):
    os.chdir(job_dir)
    # print(os.listdir())
    if not os.path.exists(os.path.join(job_dir, exe)):
        sys.exit()
    # 
    #  discard output in case too big out.
    if rank != 0:
        with open(os.devnull, 'w') as nullfile:
            subprocess.run("./" + exe, stdout = nullfile, stderr= nullfile)
    else: 
        subprocess.run("./" + exe)




# Function to convert nested dictionary to string
def dict_to_str(d, indent=0):
    result = ''
    for key, value in d.items():
        if isinstance(value, dict):
            result += ' ' * indent + key + '\n'
            result += dict_to_str(value, indent + 4)
            result += ' ' * indent + '/\n'
        elif isinstance(value, list):
            result += ' ' * (indent + 4) + f"{key} = {', '.join(map(str, value))}\n"
        else:
            result += ' ' * (indent + 4) + f"{key} = {value}\n"
    return result

# # Write the qrpa_inp dict to a file
def write_qrpa_inp_ToJobs(job_dir, qrpa_input, exe):
    # Writing the dictionary to the output file with more spaces and aligned comments
    with open(os.path.join(job_dir, exe), 'w') as output_file:
        output_file.write("============== qrpa.inp =====================================\n")

        max_length = max(len(comment) for comment in qrpa_input.keys())
        for comment, values in qrpa_input.items():
            values_str = ', '.join(map(str, values))
            spaces = ' ' * (max_length - len(comment))
            output_file.write(f"{values_str}  !{spaces}{comment}\n")

        output_file.write("=============================================================\n")

    # print(f"Output written to {output_file_path}")



def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()


    # record number of omega points, since we are gonna to change the value in dictionary. 
    # -1 because python start with 0.   
    maxNumJobs = qrpa_input['qrpa_points'][0] -1
    qrpaOmegaStart_re = qrpa_input["line parameters"][0]
    qrpaOmegaStart_im = qrpa_input["line parameters"][1]
    qrpaOmegaStep_re = qrpa_input["line parameters"][2]
    qrpaOmegaStep_im = qrpa_input["line parameters"][3]

    operatorType = str(qrpa_input["O,T,L,K"][1]) + str(qrpa_input["O,T,L,K"][2]) + str(qrpa_input["O,T,L,K"][3])

    currentDir = os.getcwd()	
    # print(currentDir)

    if(rank == 0):
        working_dir = mkdirWorking(currentDir, lpfam_dir["working"])
        output_dir = mkdirWorking(currentDir, lpfam_dir["output"])
        trainingDataOpeator_dir = mkdirWorking(output_dir, rbm_dir["trainingData_dir"]+operatorType)

    # get working_dir once one time for every subprocess. 
    # working_dir = os.path.join(currentDir, lpfam_dir["working"])
        for dest_rank in range(1, maxNumJobs + 1):
            comm.send(working_dir, dest_rank)
            comm.send(output_dir, dest_rank)
            comm.send(trainingDataOpeator_dir, dest_rank)
    else:
        if (rank <=maxNumJobs):
            working_dir = comm.recv(source = 0)
            output_dir = comm.recv(source = 0)
            trainingDataOpeator_dir = comm.recv(source = 0)


    comm.barrier()  # insure working dir exists before every subprocess create folder inside working dir.


    # exclude core which beyond the total amount of jobs.
    if (rank > maxNumJobs): sys.exit()

    job_dir = mkdirJobs(working_dir, str(rank))
    #
    # copy executables hfbtho_main to job_dir
    if rank  == 0:
        copyExeToJob(lpfam_dir["hfbtho_allpart"], os.path.join(currentDir, lpfam_dir["Exes"]), job_dir)
    else:
        copyExeToJob(lpfam_dir["hfbtho_diffpart"], os.path.join(currentDir, lpfam_dir["Exes"]), job_dir)
        
    # construct input file hfbtho_NAMELIST.dat, qrpa.in for hfbtho_main. Right now just simply copy them to job_dir.
    #Later need to design hfbtho_NAMELIST.dat according to input data.
    # copyInputToJob(lpInput_files["hfbtho_namelist"], os.path.join(currentDir, lpfam_dir["InputTemplate"]), job_dir)
    with open(os.path.join(job_dir, lpInput_files["hfbtho_namelist"]), 'w') as file:
        file.write(dict_to_str(hfbtho_NAMELIST))

    # Later we need to design different qrpa.in for each job.
    # copyInputToJob(lpInput_files["qrpa.inp"], os.path.join(currentDir, lpfam_dir["InputTemplate"]), job_dir)
    # modify qrpa_int dictionary for each subroutine so that it has the correct one omega qrpa_inp dict.
    qrpa_input["line parameters"][0] = qrpaOmegaStart_re + rank*qrpaOmegaStep_re
    qrpa_input["line parameters"][1] = qrpaOmegaStart_im + rank*qrpaOmegaStep_im
    qrpa_input["qrpa_points"][0] = 1
    write_qrpa_inp_ToJobs(job_dir, qrpa_input, lpInput_files["qrpa.inp"])


    if rank == 0:
        run_lpfam(job_dir, lpfam_dir["hfbtho_allpart"], rank)
    else:
        run_lpfam(job_dir, lpfam_dir["hfbtho_diffpart"], rank)
    # comm.barrier()  # insure emulatortraining.dat exists before every subprocess moved output to the destination.
        
    print("list jobdir: " ,os.listdir())
    if not os.path.exists("emulatortraining.dat"):
        print(f"emulatortraining.dat NOT exists in job {rank}, exits! ")
        sys.exit()
    else:
        print(rank, "done!")
        trainingDataName = os.path.join(trainingDataOpeator_dir, "emulatortraining.dat" + str(rank).zfill(3))
        shutil.copy2(os.path.join(job_dir, "emulatortraining.dat"), trainingDataName)

    

if __name__ == "__main__":
    main()
