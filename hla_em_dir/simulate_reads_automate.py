from random import sample
import pandas as pd
import subprocess
import os
import re


NUM_TO_GENERATE = 5
NUM_READS = 1500
ERROR_RATE = .01
best_case = True

# Path to the local directory open_docker.sh is in and that run_all_bash.sh writes to
reference_dir = os.path.join(os.getcwd(), 'reference')
hla_gen_path = os.path.join(os.getcwd(), 'hla_gen.fasta')
patient_filepath = "TCGA_HLA_alleles.tsv"
allele_record_filepath = os.path.join(reference_dir, 'allele_record.csv')
run_all_script = os.path.join(reference_dir, "run_all_bash.sh")

if not os.path.exists(reference_dir):
    os.mkdir(reference_dir)

bash_text = """
#!/bin/bash
REFDIR=[REF_DIR]/fastas
OUTDIR=[REF_DIR]/samples/[TRIAL_NUM]

echo 'simulating data' && wgsim -e [ERROR_RATE] -r 0 -h  -N [NUM_READS] -1 151 -2 151 -d 500 $REFDIR/reference_[TRIAL_NUM].fa $OUTDIR/sim.HLA.reads_01.fq $OUTDIR/sim.HLA.reads_02.fq || exit 1;
"""

if os.path.exists(allele_record_filepath):
    full_allele_record_df = pd.read_csv(allele_record_filepath)
else:
    full_allele_record_df = None

with open(run_all_script, "w") as run_all_bash_script:
    bash_dir = os.path.join(reference_dir, "bash")

    # run_all_bash_script.write(f"chmod 700 {os.path.join(reference_dir, '*.sh')}\n\n")
    run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', '*')}\n")
    run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', 'fastas', '*')}\n")
    run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', 'bash', '*')}\n")

    update_df = {}
    update_df["reference_fasta"] = []
    update_df["sample_01"] = []
    update_df["sample_02"] = []

    # Read in master fasta file
    with open(hla_gen_path) as hla_gen:
        contents = hla_gen.read()

    for i in range(NUM_TO_GENERATE):
        fasta_dir = os.path.join(reference_dir, "fastas")

        if not os.path.isdir(fasta_dir):
            os.makedirs(fasta_dir)

        num_fastas = len(os.listdir(fasta_dir))
        dir_name = f"trial_{num_fastas}"

        fasta_path = os.path.join(fasta_dir, f"reference_{dir_name}.fa")
        samples_dir = os.path.join(reference_dir, "samples", dir_name)
        read01_path = os.path.join(samples_dir, "sim.HLA.reads_01.fq")
        read02_path = os.path.join(samples_dir, "sim.HLA.reads_02.fq")

        if not os.path.isdir(samples_dir):
            os.makedirs(samples_dir)
        if not os.path.isdir(bash_dir):
            os.makedirs(bash_dir)

        update_df["reference_fasta"].append(fasta_path)
        update_df["sample_01"].append(read01_path)
        update_df["sample_02"].append(read02_path)

        print(update_df['sample_01'][i])
        print(update_df['sample_02'][i])

        # Read in patient HLA data
        df = pd.read_csv(patient_filepath, sep="\t")
        df.drop(columns=["Reads", "Objective", "Aliquot", "Cancer"], inplace=True)

        # Choose random patient from the list
        random_patient = df.sample().iloc[0]

        hla_alleles = []
        for hla_gene, allele in random_patient.iteritems():
            if hla_gene not in update_df:
                update_df[hla_gene] = []
                update_df[f"{hla_gene}_allele"] = []

            update_df[hla_gene].append(allele)
            pattern = r'(>HLA:\w*\s' + allele.replace("*", "\*") + r'.*?\s.*\n[ACGT\n]*)'
            all_matching_hla_alleles = re.findall(pattern, contents)

            if best_case:
                # To sample a random allele (worst case scenario) based on the TCGA patient use the following line
                hla_allele = all_matching_hla_alleles[0]
            else:
                # To sample the most common allele (best case scenario) based on the TCGA patient use the following line
                hla_allele = sample(all_matching_hla_alleles, 1)[0]

            hla_alleles.append(hla_allele)
            update_df[f"{hla_gene}_allele"].append(re.split(r"\s\d+\sbp\n", hla_allele)[0])

        # Write each reference gene to fasta file
        with open(fasta_path, 'w') as fasta:
            for hla in hla_alleles:
                fasta.write(hla)

        # Write bash script file to simulate reads
        bash_path = os.path.join(bash_dir, f"sim_reads_{dir_name}.sh")
        with open(bash_path, 'w') as bash:
            bash.write(bash_text.replace("[TRIAL_NUM]", dir_name).replace("[NUM_READS]", str(NUM_READS)).replace("[ERROR_RATE]", str(ERROR_RATE)).replace("[REF_DIR]", 'scripts'))

        # Keep track of paths to simulated reads to a text file
        with open(os.path.join(reference_dir, "sample_paths.txt"), "a") as rd:
            rd.write(read01_path + "\n")
            rd.write(read02_path + "\n")

        # Because the docker mounts the reference directory to a new directory inside the container (called "scripts"),
        # the paths to the bash script are altered here
        bash_path = os.path.join('scripts', 'bash', f"sim_reads_{dir_name}.sh")
        run_all_bash_script.write(f"bash {bash_path}" + "\n")

    # Save the allele record for each randomly selected patient to a csv file
    allele_record_update = pd.DataFrame(update_df)

    if full_allele_record_df is not None:
        full_allele_record_df = pd.concat([full_allele_record_df, allele_record_update])
    else:
        full_allele_record_df = allele_record_update

    full_allele_record_df.to_csv(allele_record_filepath, index=False)

subprocess.run(["chmod", "700", run_all_script])

# Pull wgsim docker image
subprocess.run(["docker", "pull", "pegi3s/wgsim"])

# Create a container for the docker image mounted to the reference directory
# Inside this container, run the bash script that will simulate reads
wgsim_command = ["docker", "run", "--rm", "-v", f"{reference_dir}:/scripts", "pegi3s/wgsim:latest", "/bin/bash", "-c", "bash ./scripts/run_all_bash.sh"]
result = subprocess.run(wgsim_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

if result.returncode == 0:
    print("finished simulating reads")
else:
    print("failed to simulate reads")
    print(result.stdout)
    print(result.stderr)
    exit(1)
