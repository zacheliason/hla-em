from collections import Counter
import pandas as pd
import traceback
import copy
import re
import os


def split_hla(df):
	series = df.transpose().iloc[:, 0]
	split_values = series.str.split()
	return split_values.str[0], split_values.str[1]

def calculate_match_percentage(list1, list2):
	if len(list1) != len(list2):
		raise ValueError("Input lists must have the same length")

	total_items = len(list1)
	match_counts = [0] * 4

	for allele1, allele2 in zip(list1, list2):
		tiers1 = allele1.split(':')
		tiers2 = allele2.split(':')

		for i in range(min(len(tiers1), len(tiers2))):
			if tiers1[i] == tiers2[i]:
				match_counts[i] += 1
			else:
				break

		if i < 3 and tiers1[i] == tiers2[i]:
			for j in range(i + 1, len(tiers1) + 2):
				match_counts[j] += 1

	match_percentages = [count / total_items * 100 for count in match_counts]
	return match_percentages



alleles_path = "/Users/zeliason/Desktop/hla-em/src/reference/allele_record.csv"
output_dir = "/Users/zeliason/Desktop/hla-em/src/reference/output"

# Use the following directories if testing locally
# ALLELES_PATH = "/Users/zacheliason/HLA/reference/allele_record.csv"
# RESULTS_DIR = "/Users/zacheliason/HLA"

# df contains real references simulated data is based on for each run
alleles_df = pd.read_csv(alleles_path)
alleles_df['Trial'] = alleles_df.index
hla_df = alleles_df[["A1_allele", "A2_allele", "B1_allele", "B2_allele", "C1_allele", "C2_allele"]]
hla_df['Trial'] = hla_df.index

# clean df
# for col in hla_df.columns[1:]:
#     if col == 'Trial':
#         continue
#     hla_df[col] = hla_df[col].str.replace(">", "")

for trial_dir in os.listdir(output_dir):
	trial_dir = os.path.join(output_dir, trial_dir)
	if not os.path.isdir(trial_dir):
		continue

	# Get the trial number
	trial_num = int(re.search(r"trial_(\d+)", trial_dir).group(1))
	trial_results = os.listdir(trial_dir)

	# find results.tsv in trial_results concisely
	results_tsv = [x for x in trial_results if x.endswith("results.tsv")]
	if len(results_tsv) != 1:
		print(f"Error: {trial_dir} does not contain exactly one results.tsv file")
		continue
	else:
		results_tsv = results_tsv[0]

	results = pd.read_csv(os.path.join(trial_dir, results_tsv), sep="\t")

	hla_obs = hla_df[hla_df['Trial'] == trial_num].drop(columns={"Trial"})
	hla_exp = alleles_df[alleles_df['Trial'] == trial_num][hla_obs.columns]

	hla_type_obs, hla_code_obs = split_hla(hla_obs)
	hla_type_exp, hla_code_exp = split_hla(hla_exp)

	look = calculate_match_percentage(hla_code_obs, hla_code_exp)

	print()

print()



#
# DELIM = "\t"
# with open(f"{RESULTS_DIR}/test_results.tsv", "w") as tsv:
#     # Begin writing test_results.tsv
#
#     tsv.write("Trial" + DELIM + "Converged" + DELIM + "Accuracy" + DELIM + "First Tier Count" + DELIM + "Second Tier Count" + DELIM + "Third Tier Count" + DELIM + "Perfect Match Count" +  DELIM + "A1 obs" + DELIM + "A2 obs" + DELIM + "B1 obs" + DELIM + "B2 obs" + DELIM + "C1 obs" + DELIM + "C2 obs" + DELIM + "A1 real" + DELIM + "A2 real" + DELIM + "B1 real" + DELIM + "B2 real" + DELIM + "C1 real" + DELIM + "C2 real" + DELIM + "Notes" + DELIM + "Number of Reads" + DELIM + "Number of Ref Alleles" + DELIM + "Output Path" + DELIM + "CPU Time" + DELIM + "Run Time" + DELIM + "Average Memory" + DELIM + "Max Memory" + DELIM + "Total Memory" + DELIM + "Runs" + DELIM + "Path\n")
#     for results_dir in sorted(os.listdir(RESULTS_DIR)):
#         try:
#             orig_results = copy.copy(results_dir)
#             results_dir = os.path.join(RESULTS_DIR, results_dir)
#
#             if not os.path.isdir(results_dir):
#                 continue
#
#             # Find sample path based on current result
#             sample_path = "/storage1/fs1/jin.zhang/Active/HLA-EM/Ensemble_Genome/samples/" + orig_results.replace("_fq", ".fq").replace("-", "/").replace("paired_reads", "").replace("sim_HLA_reads", "sim.HLA.reads").replace("_simulated_data", "simulated_data")
#
#             # Find current trial number
#             try:
#                 trial_num = int(re.findall(r"trial_(\d+).*", results_dir)[0])
#             except:
#                 trial_num = "R"
#
#             # # TODO zach
#             # trial_num = 16
#
#             # Initialize values
#             runs = []
#             real_hlas_sorted = ["","","","","",""]
#             obs_hlas_sorted = ["","","","","",""]
#             notes = ""
#             accuracy = 0
#             perfect_matches = 0
#             accuracy_dict = {
#                 1: 0,
#                 2: 0,
#                 3: 0,
#                 'perfect': 0
#             }
#             second_tier_matches = 0
#             imperfect_match_list = []
#             cpu_time = ""
#             run_time = ""
#             average_memory = ""
#             max_memory = ""
#             total_memory = ""
#             top_genes = ""
#             reads = ""
#             genes_num = ""
#             converged = "TRUE"
#             lettersdict = {
#                 "A1": "",
#                 "A2": "",
#                 "B1": "",
#                 "B2": "",
#                 "C1": "",
#                 "C2": "",
#             }
#
#             results_path = os.path.join(results_dir, "hla_em_fullrun", "v1.results.tsv")
#             out_path = os.path.join(results_dir, f"TRIAL{trial_num}_HLA-EM.out")
#             err_path = os.path.join(results_dir, f"TRIAL{trial_num}_HLA-EM.err")
#
#             # Use following outpath if testing locally
#             # out_path = "/Users/zacheliason/HLA/output.txt"
#
#             try:
#                 with open(out_path) as out:
#                     pass
#             except:
#                 out_path = os.path.join(results_dir, f"HLA_EM_fullrun.out")
#                 err_path = os.path.join(results_dir, f"HLA_EM_fullrun.err")
#
#             # Begin parsing outpath for results
#             with open(out_path) as out:
#                 contents = out.read()
#
#                 m = re.search(r'CPU time\s+:\s+(\d.*) sec', contents)
#                 if m:
#                     cpu_time = m.group(1)
#                 m = None
#
#                 m = re.search(r'Max Memory\s+:\s+(\d.*) MB', contents)
#                 if m:
#                     max_memory = m.group(1)
#                 m = None
#
#                 m = re.search(r'Average Memory\s*:\s+(\d.*) MB', contents)
#                 if m:
#                     average_memory = m.group(1)
#                 m = None
#
#                 m = re.search(r'Total Requested Memory\s+:\s+(\d.*) MB', contents)
#                 if m:
#                     total_memory = m.group(1)
#                 m = None
#
#                 m = re.search(r'Run time\s*:\s*(\d.*)\s*s', contents)
#                 if m:
#                     run_time = m.group(1)
#                 m = None
#
#                 m = re.search(r'Total reads: (\d+)', contents)
#                 if m:
#                     reads = str(m.group(1))
#
#                 pattern = re.compile(r"Iteration ([abc]{1,3}-\d\**)")
#                 for match in pattern.finditer(contents):
#                     runs.append(match.group(1))
#
#                 pattern = re.compile(r"Results from iteration: ([ABC]{1,3}-\d\**)")
#                 for match in pattern.finditer(contents):
#                     passed_run = match.group(1).lower()
#                     runs = [x.upper() if x == passed_run else x for x in runs]
#
#                 pattern = re.compile(r'num genes after filtering by breadth: (\d+)')
#                 gene_lines = []
#                 for match in pattern.finditer(contents):
#                     gene_lines.append(match.group(1))
#
#                 if len(gene_lines) > 0:
#                     genes_num = str(gene_lines[-1])
#
#                 if "Successfully completed" in contents:
#                     converged = "TRUE"
#                 elif "Exited with exit code 1." in contents:
#                     converged = "FALSE"
#                 else:
#                     converged = "PENDING"
#
#                 try:
#                     real_hlas = hla_df[hla_df['Trial'] == trial_num].values.tolist()[0][1:]
#                 except:
#                     print(traceback.format_exc())
#                     print(trial_num)
#                     print()
#                     raise
#
#                 hlas_sorted = []
#                 gene_groups = []
#                 classifications = []
#                 for allele in real_hlas:
#                     m = re.search(r"HLA:HLA\d+\s(([ABC])\*[\d:\w]*)", allele)
#                     if m:
#                         gg = m.group(2)
#                         gene_groups.append(gg)
#                         classifications.append(m.group(1))
#                         hlas_sorted.append(allele)
#
#                 duplicates = []
#                 dup_genes = []
#                 dup_class = []
#                 counter = Counter(gene_groups)
#                 for k, v in counter.items():
#                     if v < 2:
#                         i = gene_groups.index(k)
#                         duplicates.append(hlas_sorted[i])
#                         dup_class.append(classifications[i])
#                         dup_genes.append(k)
#
#                 hlas_sorted.extend(duplicates)
#                 gene_groups.extend(dup_genes)
#                 classifications.extend(dup_class)
#
#                 gene_groups, classifications, hlas_sorted = zip(*sorted(zip(gene_groups, classifications, hlas_sorted)))
#
#                 gene_groups = list(gene_groups)
#                 real_hlas_sorted = list(hlas_sorted)
#
#             if converged == "TRUE":
#                 m = re.search(r"((FINAL RESULTS|Completing iterations.*)\n(.*\n)*?((HLA.*\n)*)\n*-------------------------)", contents)
#                 if m:
#                     final_results = m.group(1)
#
#                     pattern = re.compile(r'(HLA:HLA\d*\s((?:\(([ABC]).*\))|(?:([ABC]).*==))).*')
#
#                     allele_ids = []
#                     letters = []
#
#                     for match in pattern.finditer(final_results):
#                         allele_id = match.group(1).replace("(", "").replace(")", "").replace("=", "")
#                         n = re.search(r"\(*([ABC])\*.*", allele_id)
#                         if n:
#                             letter = n.group(1)
#
#                         letters.append(letter)
#                         allele_ids.append(allele_id)
#
#                     letters, allele_ids = zip(*sorted(zip(letters, allele_ids)))
#
#                     seen_letters = []
#                     for i in range(len(letters)):
#                         letter = letters[i]
#                         allele_id = allele_ids[i]
#
#                         if seen_letters.count(letter) == 0:
#                             lettersdict[f"{letter}1"] = allele_id
#                         elif seen_letters.count(letter) == 1:
#                             if allele_id != lettersdict[f"{letter}1"]:
#                                 lettersdict[f"{letter}2"] = allele_id
#                             else:
#                                 continue
#
#                         seen_letters.append(letter)
#
#                     # Copy values over if needed
#                     seconds = ['A2', 'B2', 'C2']
#                     for s in seconds:
#                         if lettersdict[s] == "":
#                             lettersdict[s] = lettersdict[s.replace("2","1")]
#
#                     obs_hlas = list(lettersdict.values())
#
#                     hlas_sorted = []
#                     gene_groups = []
#                     classifications = []
#                     for allele in obs_hlas:
#                         m = re.search(r"HLA:HLA\d+\s(([ABC])\*[\d:\w]*)", allele)
#                         if m:
#                             gg = m.group(2)
#                             gene_groups.append(gg)
#                             classifications.append(m.group(1))
#                             hlas_sorted.append(allele)
#
#                     gene_groups, classifications, hlas_sorted = zip(*sorted(zip(gene_groups, classifications, hlas_sorted)))
#                     obs_hlas_sorted = list(hlas_sorted)
#
#                     if len(obs_hlas_sorted) < 6:
#                         continue
#
#                     print(f"NOW WORKING ON TRIAL {trial_num}")
#                     for i in range(0, 6, 2):
#                         o = obs_hlas_sorted[i:i+2]
#                         r = real_hlas_sorted[i:i+2]
#
#                         o0 = [x for x in re.split(r"([:*])", re.split(r"HLA.*\s", o[0])[1]) if x != "" and x != "*" and x != ":"][1:]
#                         r0 = [x for x in re.split(r"([:*])", re.split(r"HLA.*\s", r[0])[1]) if x != "" and x != "*" and x != ":"][1:]
#                         o1 = [x for x in re.split(r"([:*])", re.split(r"HLA.*\s", o[1])[1]) if x != "" and x != "*" and x != ":"][1:]
#                         r1 = [x for x in re.split(r"([:*])", re.split(r"HLA.*\s", r[1])[1]) if x != "" and x != "*" and x != ":"][1:]
#
#                         for j in range(min(len(o0), len(r0))):
#                             o0j = o0[j]
#                             r0j = r0[j]
#                             if o0[j] != r0[j] and o1[j] != r1[j]:
#                                 temp = obs_hlas_sorted[i]
#                                 obs_hlas_sorted[i] = obs_hlas_sorted[i + 1]
#                                 obs_hlas_sorted[i + 1] = temp
#                                 break
#                             else:
#                                 break
#
#
#
#                     for i in range(len(obs_hlas_sorted)):
#                         o_hla = obs_hlas_sorted[i]
#                         r_hla = real_hlas_sorted[i]
#                         if o_hla == r_hla:
#                             r = re.search(r"HLA.*([ABC].*)", r_hla)
#                             if r:
#                                 r_hla_tiers = [x for x in re.split(r"[:\*]", r.group(1)) if x not in ["A", "B", "C"]]
#                                 for z in range(len(r_hla_tiers)):
#                                     if z + 1 in accuracy_dict:
#                                         accuracy_dict[z + 1] += 1
#
#                             perfect_matches += 1
#                             accuracy_dict['perfect'] += 1
#                             second_tier_matches += 1
#                             imperfect_match_list.append(1)
#                         else:
#                             current_progress = 0
#                             o = re.search(r"HLA.*([ABC].*)", o_hla)
#                             if o:
#                                 o_hla_tiers = [x for x in re.split(r"[:\*]", o.group(1)) if x not in ["A", "B", "C"]]
#                             r = re.search(r"HLA.*([ABC].*)", r_hla)
#                             if r:
#                                 r_hla_tiers = [x for x in re.split(r"[:\*]", r.group(1)) if x not in ["A", "B", "C"]]
#
#                             for j in range(min(len(r_hla_tiers), len(o_hla_tiers))):
#                                 o = o_hla_tiers[j]
#                                 r = r_hla_tiers[j]
#
#                                 if r == o:
#                                     if len(r_hla_tiers) != j + 1:
#                                         accuracy_dict[j + 1] += 1
#                                     if j == 1:
#                                         # second tier matches!
#                                         second_tier_matches += 1
#                                     current_progress += 1
#                                 else:
#                                     break
#
#                             imperfect_match_list.append(current_progress / max(len(o_hla_tiers), len(r_hla_tiers)))
#
#                     if len(imperfect_match_list) > 0:
#                         accuracy = sum(imperfect_match_list) / len(imperfect_match_list)
#
#             with open(err_path) as err:
#                 if converged == "FALSE":
#                     err_message = err.readlines()[-1]
#                     notes += err_message.replace("\n", " ")
#
#             if type(trial_num) != str:
#                 trial_num = str(trial_num)
#
#             tsv.write(f"Trial {trial_num}" +
#                 DELIM + converged +
#                 DELIM + str(accuracy * 100) + "%" +
#                 DELIM + str(accuracy_dict[1]) +
#                 DELIM + str(accuracy_dict[2]) +
#                 DELIM + str(accuracy_dict[3]) +
#                 DELIM + str(accuracy_dict['perfect']) +
#                 DELIM + DELIM.join(obs_hlas_sorted) +
#                 DELIM + DELIM.join(real_hlas_sorted) +
#                 DELIM + notes +
#                 DELIM + reads +
#                 DELIM + genes_num +
#                 DELIM + out_path +
#                 DELIM + cpu_time + " sec" +
#                 DELIM + run_time + " sec" +
#                 DELIM + average_memory + " MB" +
#                 DELIM + max_memory + " MB" +
#                 DELIM + total_memory + " MB" +
#                 DELIM + ", ".join(runs) +
#                 DELIM + results_path + "\n")
#
#
#         except:
#             print(traceback.format_exc())
#             print(f"{trial_num} FAILED")
#             tsv.write(f"Trial {trial_num}" + DELIM + "PENDING" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + DELIM + "" + "\n")
#
# df = pd.read_csv(f"{RESULTS_DIR}/test_results.tsv", sep="\t")
#
# new_c = list(df.columns.values.tolist()[1:])
# new_c.append("")
#
# nc = dict(zip(df.columns.values.tolist(), new_c))
# df = df.rename(columns=nc)
# df = df.dropna(subset=['Accuracy'])
#
# try:
#     accuracy_avg = str(sum([float(x) for x in df["Accuracy"].astype('str').str.replace("%", "").values.tolist()]) / len(df['Accuracy'].values.tolist()))
# except:
#     print("ERR in accuracy")
#     print(df["Accuracy"].values.tolist())
#
# first_avg = str(df["First Tier Count"].astype('int').mean())
# second_avg = str(df["Second Tier Count"].astype('int').mean())
# third_avg = str(df["Third Tier Count"].astype('int').mean())
# perfect_avg = str(df["Perfect Match Count"].astype('int').mean())
#
# cpu_avg = sum([float(x.strip()) for x in df["CPU Time"].str.replace("sec", "").values.tolist()]) / len(df['CPU Time'].values.tolist())
# runtime_avg = sum([float(x.strip()) for x in df["Run Time"].str.replace("sec", "").values.tolist()]) / len(df['Run Time'].values.tolist())
# avg_memory_avg = sum([float(x.strip()) for x in df["Average Memory"].str.replace("MB", "").values.tolist()]) / len(df['Average Memory'].values.tolist())
# max_memory_avg = sum([float(x.strip()) for x in df["Max Memory"].str.replace("MB", "").values.tolist()]) / len(df['Max Memory'].values.tolist())
# total_memory_avg = sum([float(x.strip()) for x in df["Total Memory"].str.replace("MB", "").values.tolist()]) / len(df['Total Memory'].values.tolist())
#
# with open(f"{RESULTS_DIR}/test_results.tsv", "a") as tsv:
#     tsv.write(f"AVERAGES" +
#           DELIM + "" +
#           DELIM + accuracy_avg + "%" +
#           DELIM + first_avg +
#           DELIM + second_avg +
#           DELIM + third_avg +
#           DELIM + perfect_avg +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM +
#           DELIM + str(round(cpu_avg,2)) + " sec" +
#           DELIM + str(round(runtime_avg,2)) + " sec" +
#           DELIM + str(round(avg_memory_avg,2)) + " MB" +
#           DELIM + str(round(max_memory_avg,2)) + " MB" +
#           DELIM + str(round(total_memory_avg,2)) + " MB" +
#           DELIM +
#           DELIM + "\n")
