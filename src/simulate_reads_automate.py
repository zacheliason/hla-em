from random import sample
import pandas as pd
import subprocess
import random
import shutil
import os
import re


# NUM_TO_GENERATE = 5
# NUM_READS = 1500
# ERROR_RATE = .01
# best_case = True


# Path to the local directory open_docker.sh is in and that run_all_bash.sh writes to
PARENT_DIR = os.path.join(os.path.dirname(os.getcwd()))
REFERENCE_DIR = os.path.join(PARENT_DIR, 'reference')


def apply_loss_of_heterozygousity(alleles, num_lost_alleles):
    # Group the Series by the first character of the index
    allele_groups = alleles.groupby(alleles.index.str[0])

    random_allele_groups = random.sample(sorted(allele_groups.groups.keys()), k=min(num_lost_alleles, len(allele_groups)))

    updated_alleles = alleles.copy()
    for random_group in random_allele_groups:
        # Get the indices of the selected groups
        indices = [index for index in allele_groups.get_group(random_group).index]

        # Select a random index from the selected groups
        index_to_drop = random.choice(indices)

        # Drop the randomly selected index
        # updated_alleles = updated_alleles.drop(index_to_drop)
        updated_alleles[index_to_drop] = None

    return updated_alleles


def simulate_hla_reads(num_to_generate, num_reads, error_rate=.01, best_case=True, loh_number=0):
    hla_gen_path = os.path.join(PARENT_DIR, 'hla_gen.fasta')
    patient_filepath = os.path.join(REFERENCE_DIR, "TCGA_HLA_alleles.tsv")
    allele_record_filepath = os.path.join(REFERENCE_DIR, 'allele_record.csv')
    run_all_script = os.path.join(REFERENCE_DIR, "run_all_bash.sh")

    if not os.path.exists(REFERENCE_DIR):
        os.mkdir(REFERENCE_DIR)

    hla_bash_text = """
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
        bash_dir = os.path.join(REFERENCE_DIR, "bash")

        run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', '*')}\n")
        run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', 'fastas', '*')}\n")
        run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', 'bash', '*')}\n")

        update_df = {}
        update_df["best_case_used"] = [best_case] * num_to_generate
        update_df['num_hla_reads'] = [num_reads] * num_to_generate
        update_df['hla_reads_percentage'] = [num_reads / 1000000] * num_to_generate
        update_df["reference_fasta"] = []
        update_df["sample_01"] = []
        update_df["sample_02"] = []

        # Read in master fasta file
        with open(hla_gen_path) as hla_gen:
            contents = hla_gen.read()

        for i in range(num_to_generate):
            fasta_dir = os.path.join(REFERENCE_DIR, "fastas")

            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)

            num_fastas = len(os.listdir(fasta_dir))
            dir_name = f"trial_{num_fastas}"

            fasta_path = os.path.join(fasta_dir, f"reference_{dir_name}.fa")
            samples_dir = os.path.join(REFERENCE_DIR, "samples", dir_name)
            read01_path = os.path.join(samples_dir, "sim.HLA.reads_01.fq")
            read02_path = os.path.join(samples_dir, "sim.HLA.reads_02.fq")

            if not os.path.isdir(samples_dir):
                os.makedirs(samples_dir)
            if not os.path.isdir(bash_dir):
                os.makedirs(bash_dir)

            update_df["reference_fasta"].append(fasta_path)
            update_df["sample_01"].append(read01_path)
            update_df["sample_02"].append(read02_path)

            # Read in patient HLA data
            df = pd.read_csv(patient_filepath, sep="\t")
            df.drop(columns=["Reads", "Objective", "Aliquot", "Cancer"], inplace=True)

            # Choose random patient from the list
            random_patient = df.sample().iloc[0]

            # Apply loss of heterozygosity
            random_patient = apply_loss_of_heterozygousity(random_patient, loh_number)

            hla_alleles = []
            for hla_gene, allele in random_patient.items():

                if hla_gene not in update_df:
                    update_df[hla_gene] = []
                    update_df[f"{hla_gene}_allele"] = []

                if allele is None:
                    update_df[f"{hla_gene}"].append(allele)
                    update_df[f"{hla_gene}_allele"].append(allele)
                    continue

                update_df[hla_gene].append(allele)
                pattern = r'(>HLA:\w*\s' + allele.replace("*", r"\*") + r'.*?\s.*\n[ACGT\n]*)'
                all_matching_hla_alleles = re.findall(pattern, contents)

                if best_case:
                    # sample a random allele (worst case scenario) based on the TCGA patient
                    hla_allele = all_matching_hla_alleles[0]
                else:
                    # sample the most common allele (best case scenario) based on the TCGA patient
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
                bash.write(hla_bash_text.replace("[TRIAL_NUM]", dir_name).replace("[NUM_READS]", str(num_reads)).replace("[ERROR_RATE]", str(error_rate)).replace("[REF_DIR]", 'scripts'))

            # Keep track of paths to simulated reads to a text file
            with open(os.path.join(REFERENCE_DIR, "sample_paths.txt"), "a") as rd:
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

    # Create a container for the docker image mounted to the reference directory
    # Inside this container, run the bash script that will simulate reads
    wgsim_command = ["docker", "run", "--rm", "-v", f"{REFERENCE_DIR}:/scripts", "pegi3s/wgsim:latest", "/bin/bash", "-c", "bash ./scripts/run_all_bash.sh"]
    result = subprocess.run(wgsim_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode == 0:
        # clean up
        shutil.rmtree(bash_dir)
        os.remove(run_all_script)
    else:
        print(result.stdout)
        print(result.stderr)
        exit(1)


def simulate_masked_reads(num_to_generate, num_reads, error_rate=.01):
    masked_human_genome_path = os.path.join(os.getcwd(), 'Genome.fa')
    run_all_script = os.path.join(REFERENCE_DIR, "run_all_bash.sh")

    masked_human_genome_bash_text = """
    #!/bin/bash
    REFDIR=[REF_DIR]
    OUTDIR=[REF_DIR]/reference/samples/[TRIAL_NUM]

    echo 'simulating data' && wgsim -e [ERROR_RATE] -r 0 -h  -N [NUM_READS] -1 151 -2 151 -d 500 $REFDIR/MaskedHumanGenome.fa $OUTDIR/sim.masked.reads_01.fq $OUTDIR/sim.masked.reads_02.fq || exit 1;
    cat $OUTDIR/sim.masked.reads_01.fq >> $OUTDIR/sim.HLA.reads_01.fq
    cat $OUTDIR/sim.masked.reads_02.fq >> $OUTDIR/sim.HLA.reads_02.fq
    # rm $OUTDIR/sim.masked.reads_01.fq
    # rm $OUTDIR/sim.masked.reads_02.fq
    """

    if not os.path.exists(REFERENCE_DIR):
        os.mkdir(REFERENCE_DIR)

    with open(run_all_script, "w") as run_all_bash_script:
        bash_dir = os.path.join(REFERENCE_DIR, "bash")

        run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', 'reference', '*')}\n")
        run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', 'reference', 'fastas', '*')}\n")
        run_all_bash_script.write(f"chmod 700 {os.path.join('scripts', 'reference', 'bash', '*')}\n")

        sample_dirs = []
        for i in range(num_to_generate):
            curr_sample = len(os.listdir(os.path.join(REFERENCE_DIR, 'fastas'))) - num_to_generate + i
            dir_name = f"trial_{curr_sample}"

            samples_dir = os.path.join(REFERENCE_DIR, "samples", dir_name)
            sample_dirs.append(samples_dir)

            if not os.path.isdir(samples_dir):
                raise Exception(f"Directory {samples_dir} does not exist")
            if not os.path.isdir(bash_dir):
                os.makedirs(bash_dir)

            # Write bash script file to simulate reads
            bash_path = os.path.join(bash_dir, f"sim_reads_{dir_name}.sh")
            with open(bash_path, 'w') as bash:
                bash.write(masked_human_genome_bash_text.replace("[TRIAL_NUM]", dir_name).replace("[NUM_READS]", str(num_reads)).replace("[ERROR_RATE]", str(error_rate)).replace("[REF_DIR]", 'scripts'))

            # Because the docker mounts the reference directory to a new directory inside the container (called "scripts"),
            # the paths to the bash script are altered here
            bash_path = os.path.join('scripts', 'reference', 'bash', f"sim_reads_{dir_name}.sh")
            run_all_bash_script.write(f"bash {bash_path}" + "\n")

    subprocess.run(["chmod", "700", run_all_script])

    # Create a container for the docker image mounted to the reference directory
    # Inside this container, run the bash script that will simulate reads
    wgsim_command = ["docker", "run", "--rm", "-v", f"{os.getcwd()}:/scripts", "pegi3s/wgsim:latest", "/bin/bash", "-c", "bash ./scripts/reference/run_all_bash.sh"]
    result = subprocess.run(wgsim_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode == 0:
        # clean up
        for samples_dir in sample_dirs:
            if os.path.exists(os.path.join(samples_dir, 'sim.masked.reads_01.fq')):
                os.remove(os.path.join(samples_dir, 'sim.masked.reads_01.fq'))
            else:
                print(os.path.join(samples_dir, 'sim.masked.reads_01.fq'), " doesn't exist")
            if os.path.exists(os.path.join(samples_dir, 'sim.masked.reads_02.fq')):
                os.remove(os.path.join(samples_dir, 'sim.masked.reads_02.fq'))
            else:
                print(os.path.join(samples_dir, 'sim.masked.reads_01.fq'), " doesn't exist")
        shutil.rmtree(bash_dir)
        os.remove(run_all_script)
    else:
        print(result.stdout)
        print(result.stderr)
        exit(1)


def create_large_tests():
    num_masked_reads = 1000000
    test_cases = [100000, 50000, 10000, 5000, 1500]
    num_per_test_case = 10

    ctr = 0
    for num_hla_reads in test_cases:
        num_human_reads = num_masked_reads - num_hla_reads
        # if ctr == 0:
        #     ctr += 1
        #     continue

        # simulate hla reads using best case HLA alleles from TCGA data
        simulate_hla_reads(num_per_test_case // 2, num_hla_reads, best_case=True)
        # simulate hla reads using random HLA alleles from TCGA data
        simulate_hla_reads(num_per_test_case // 2, num_hla_reads, best_case=False)
        # simulate reads from the rest of the genome and concatenate them with the hla reads
        simulate_masked_reads(num_per_test_case, num_human_reads)

def create_dummy_tests():
    num_per_test_case = 3
    num_hla_reads = 1500

    current_path = os.environ.get("PATH", "")
    home_dir = os.path.expanduser("~")
    new_path = f"{current_path}:{home_dir}/.docker/bin"
    os.environ['PATH'] = new_path

    loh_number = 1
    num_reads = (num_hla_reads // 6) * (6 - loh_number)
    # simulate hla reads using best case HLA alleles from TCGA data
    simulate_hla_reads(num_per_test_case, num_reads, best_case=True, loh_number=loh_number)
    # simulate hla reads using random HLA alleles from TCGA data
    simulate_hla_reads(num_per_test_case, num_reads, best_case=False, loh_number=loh_number)

    loh_number = 2
    num_reads = (num_hla_reads // 6) * (6 - loh_number)
    simulate_hla_reads(num_per_test_case, num_reads, best_case=True, loh_number=loh_number)
    simulate_hla_reads(num_per_test_case, num_reads, best_case=False, loh_number=loh_number)

create_dummy_tests()