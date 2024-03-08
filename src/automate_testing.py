import subprocess
import traceback
import shutil
import os
import re

REF_FILE = '/storage1/fs1/jin.zhang/Active/HLA-EM/Ensemble_Genome/reference/hla_gen.fasta'
RESULTS_DIR = '/storage1/fs1/jin.zhang/Active/HLA-EM/Ensemble_Genome/results/automated_runs'

# Read in completed samples
# with open("completed_sample_paths.txt") as c:
#     completed = []
#     for line in c:
#         completed.append(line.strip())

samples_dir = "/Users/zacheliason/Documents/Work/zhang/2024/new/reference/samples"
output_dir = "/Users/zacheliason/Documents/Work/zhang/2024/new/reference/output"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sample_paths = os.listdir(samples_dir)
# sample_paths = [ for x in sample_runs]

output_runs = os.listdir(output_dir)
output_paths = [os.path.split(x)[-1] for x in output_runs]

samples_todo = [os.path.join(samples_dir, x) for x in sample_paths if os.path.split(x)[-1] not in output_paths]
print(f"length of todo: {len(samples_todo)}")

bash_lines = ""
cwd = "/Users/zacheliason/Documents/Work/zhang/2024/new"
hla_em_path = os.path.join(cwd, 'HLA_EM.py')
star_index = os.path.join(cwd, 'EnsembleGenome_STAR_without_scaffolds')
ref_path = os.path.join(cwd, 'hla_gen.fasta')
star_viral = os.path.join(cwd, 'HLA_Allele_Sorter_Output_Tier2.fa_STAR')

for sample_path in samples_todo:
    basename = os.path.basename(sample_path)
    out_dir = os.path.join(output_dir, basename)
    fq_files = [os.path.join(sample_path, x) for x in os.listdir(sample_path)]

    hla_em_cmd = ['python3', hla_em_path, '-t', '4', '-o', out_dir, '-s', star_index, '-r', ref_path] + fq_files

    result = subprocess.run(hla_em_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")

#     if len(sample_path.strip()) == 0:
#         print(sample_path)
#         continue
#
#     if "reads_02" in sample_path or "SRR701474_2.fastq" in sample_path:
#         paired_reads = True
#         # continue
#     elif "reads_01" in sample_path or "SRR701474_1.fastq" in sample_path:
#         continue
#         paired_reads = False
#
#     try:
#         if "/samples/" not in sample_path:
#             print('ERR: ' + sample_path)
#         results_dir = "-".join(sample_path.split("/samples/")[-1].split("/")).replace(".", "_").replace(".fq", "_results")
#         bash_path = "bash/" + "-".join(sample_path.split("/samples/")[-1].split("/")).replace(".", "_").replace("fq", "_hla-em.sh")
#     except:
#         print(traceback.format_exc())
#         print(sample_path)
#         raise Exception
#
#     if paired_reads:
#         results_dir = f"paired_reads_{results_dir}"
#         bash_path = bash_path.replace("bash/", "bash/paired_reads_")
#
#     # Set results directory
#     full_results_dir = os.path.join(RESULTS_DIR, results_dir)
#     if not os.path.isdir(full_results_dir):
#         os.makedirs(full_results_dir)
#         os.chmod(full_results_dir, 777)
#         print(f"MADE: {full_results_dir}")
#     else:
#         shutil.rmtree(full_results_dir)
#         os.makedirs(full_results_dir)
#         os.chmod(full_results_dir, 777)
#         print(f"REPLACED: {full_results_dir}")
#
#     # Add current .sh file to cumulative runbash.sh file
#     bash_lines += "./" + bash_path + "\n"
#
#     if not os.path.isdir("bash"):
#         os.mkdir("bash")
#
#     # Write bash file, including information specific to this run
#     b = bash_text.replace("[RESULTS_NUM]", f"/{results_dir}")
#
#     if paired_reads:
#         b = b.replace("[FQ_FILE]", f"{sample_path.replace('2.f', '1.f')} {sample_path}").replace("[OUT_DIR]", RESULTS_DIR)
#
#     m = re.search(r"Trial(\d+)", sample_path)
#     if m:
#         name = f"Trial{m.group(1)}_hla-em"
#     else:
#         name = f"hla-em"
#
#     b = b.replace("[FQ_FILE]", sample_path)
#     b = b.replace("[REF_FILE]", REF_FILE)
#     b = b.replace("[JOB_NAME]", name.upper())
#
#     with open(bash_path, "w") as bash:
#         bash.write(b)
#
#     with open("completed_sample_paths.txt", "a") as csp:
#         csp.write(sample_path + "\n")
#
#     with open("results_dirs.txt", "a") as rd:
#         rd.write(results_dir + "\n")
#
# with open(f"runbash.sh", "w") as bash:
#     bash.write("chmod 700 bash/*\n\n")
#     bash.write(bash_lines)
#
# os.chmod("runbash.sh", 700)
