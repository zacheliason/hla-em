import pandas as pd
import numpy as np
import re
import os


def count_matching_characters(code1, code2):
    arr1 = np.array(list(code1))
    arr2 = np.array(list(code2))

    # Find the minimum length of the two input strings
    min_length = min(len(code1), len(code2))

    # Compare the characters and create a boolean array
    matches = arr1[:min_length] == arr2[:min_length]

    # Count the True values (matches) in the boolean array
    count = np.sum(matches)

    return count


def clean_output(path):
    # Step 1: Read the TSV file into a DataFrame
    df = pd.read_csv(path, sep='\t')

    # Step 2: Extract the HLA type and MLE Probability columns
    df['HLAletter'] = df['HLAtype'].str.extract(r'\((\w+\d*)\*.*\)')  # Extract HLA type (A, B, or C)
    df[['HLA_type', 'HLA_code']] = df['HLAtype'].str.split(' ', expand=True)

    df = df[df['HLAletter'].isin(['A', 'B', 'C'])]


    # Step 3: Filter out alleles based on a threshold probability
    threshold = 0.001  # Adjust the threshold as needed
    df_filtered = df[df['MLE_Probability'] >= threshold]

    # Step 4: Group alleles by their HLA type (A, B, or C)
    groups = df_filtered.groupby('HLAletter')

    # Step 5: Sort alleles within each group based on MLE Probability
    sorted_groups = {key: group.sort_values(by='MLE_Probability', ascending=False) for key, group in groups}

    # Step 6: Select the top two alleles for each HLA type
    hla_predictions = []
    for group in sorted_groups.values():
        top_values = group['MLE_Probability'].nlargest(2).values.tolist()
        filtered_group = group[group['MLE_Probability'].isin(top_values)]
        for mle_prob in filtered_group['MLE_Probability'].unique():
            ng = filtered_group[filtered_group['MLE_Probability'] == mle_prob]
            if len(ng) < 2:
                hla_predictions.append(ng)
            else:
                # Define a custom sorting key function
                def custom_sort_key(code):
                    # Sort by the number of matching characters with other codes
                    matching_counts = [count_matching_characters(code, other_code) for other_code in ng['HLA_code']]
                    max_matching_count = sorted(matching_counts, reverse=True)[1]

                    # Second, sort alphabetically
                    return (-max_matching_count, code)

                # Sort the DataFrame based on the custom sorting key
                df_sorted = ng.iloc[ng['HLA_code'].map(custom_sort_key).argsort()]
                hla_predictions.append(df_sorted.head(1))

    hla_prediction_df = pd.concat(hla_predictions)
    hla_prediction_df.index = hla_prediction_df['HLAtype']
    hla_prediction_df = hla_prediction_df.drop(columns=['HLAletter', "HLA_code", "HLA_type", "HLAtype"])

    output_dir, _ = os.path.split(path)

    hla_prediction_df.to_csv(os.path.join(output_dir, "final_predictions.csv"), index=True)


def calculate_match_percentage(obs, exp):
	if len(obs) != len(exp):
		raise ValueError("Input lists must have the same length")

	total_items = len(obs)
	match_counts = [0] * 4

	for allele1, allele2 in zip(obs, exp):
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


def count_common_sections(allele1, allele2, include_length=False):
	"""
	Count the number of comma-separated sections that are common between two alleles.
	"""

	common_sections = 0
	sections1 = allele1.split(':')
	sections2 = allele2.split(':')

	for i in range(min(len(sections1), len(sections2))):
		if sections1[i] != sections2[i]:
			break
		common_sections += 1

	if include_length:
		return common_sections, len(sections2)
	else:
		return common_sections


def score_allele_lists(obs, exp):
	"""
	Score how well two lists of HLA alleles correspond to each other.
	"""

	def make_dict(lst):
		new_dict = {}
		for x in lst:
			if x[0] not in new_dict:
				new_dict[x[0]] = []
			new_dict[x[0]].append(x)
		return new_dict

	obs_dict = make_dict(obs)
	exp_dict = make_dict(exp)
	ordered_alleles = {}
	allele_group_match = {}
	protein_match = {}
	coding_match = {}
	noncoding_match = {}

	for key in obs_dict:
		obs_letter = obs_dict[key]
		exp_letter = exp_dict[key]

		if len(obs_letter) == 1 and len(exp_letter) == 1:
			ordered_alleles[f'{key}1_obs'], ordered_alleles[f'{key}1_exp'] = obs_letter[0], exp_letter[0]

		elif len(obs_letter) == 1 and len(exp_letter) > 1:
			coll = (
				count_common_sections(obs_letter[0], exp_letter[0]),
				count_common_sections(obs_letter[0], exp_letter[1]),
			)
			ordered_alleles[f'{key}1_obs'] = obs_letter[0]
			if coll[0] > coll[1]:
				ordered_alleles[f'{key}1_exp'] = exp_letter[0]
				ordered_alleles[f'{key}2_exp'] = exp_letter[1]
			else:
				ordered_alleles[f'{key}1_exp'] = exp_letter[1]
				ordered_alleles[f'{key}2_exp'] = exp_letter[0]

		elif len(obs_letter) > 1 and len(exp_letter) == 1:
			coll = (
				count_common_sections(obs_letter[0], exp_letter[0]),
				count_common_sections(obs_letter[1], exp_letter[0]),
			)
			ordered_alleles[f'{key}1_exp'] = exp_letter[0]
			if coll[0] > coll[1]:
				ordered_alleles[f'{key}1_obs'] = obs_letter[0]
				ordered_alleles[f'{key}2_obs'] = obs_letter[1]
			else:
				ordered_alleles[f'{key}1_obs'] = obs_letter[1]
				ordered_alleles[f'{key}2_obs'] = obs_letter[0]

		elif len(obs_letter) > 1 and len(exp_letter) > 1:
			coll = (
				sum([count_common_sections(obs_letter[0], exp_letter[0]),
					 count_common_sections(obs_letter[1], exp_letter[1])]),
				sum([count_common_sections(obs_letter[0], exp_letter[1]),
					 count_common_sections(obs_letter[1], exp_letter[0])])
			)

			if coll[0] > coll[1]:
				ordered_alleles[f'{key}1_obs'], ordered_alleles[f'{key}1_exp'] = obs_letter[0], exp_letter[0]
				ordered_alleles[f'{key}2_obs'], ordered_alleles[f'{key}2_exp'] = obs_letter[1], exp_letter[1]
			else:
				ordered_alleles[f'{key}1_obs'], ordered_alleles[f'{key}1_exp'] = obs_letter[0], exp_letter[1]
				ordered_alleles[f'{key}2_obs'], ordered_alleles[f'{key}2_exp'] = obs_letter[1], exp_letter[0]

		# TODO take out binary
		for num in [1, 2]:
			if f'{key}{num}_obs' not in ordered_alleles and f'{key}{num}_exp' not in ordered_alleles:
				# allele_group_match[f"{key}{num}"] = 1
				# protein_match[f"{key}{num}"] = 1
				# coding_match[f"{key}{num}"] = 1
				# noncoding_match[f"{key}{num}"] = 1
				continue
			elif f'{key}{num}_obs' not in ordered_alleles or f'{key}{num}_exp' not in ordered_alleles:
				allele_group_match[f"{key}{num}"] = 0
				protein_match[f"{key}{num}"] = 0
				coding_match[f"{key}{num}"] = 0
				noncoding_match[f"{key}{num}"] = 0
				continue

			obs_allele = ordered_alleles[f'{key}{num}_obs']
			exp_allele = ordered_alleles[f'{key}{num}_exp']

			num_matching, length_exp = count_common_sections(obs_allele, exp_allele, include_length=True)
			if length_exp == 1:
				allele_group_match[f"{key}{num}"] = 1 if num_matching > 0 else 0
				protein_match[f"{key}{num}"] = 1 if num_matching > 0 else 0
				coding_match[f"{key}{num}"] = 1 if num_matching > 0 else 0
				noncoding_match[f"{key}{num}"] = 1 if num_matching > 0 else 0
			if length_exp == 2:
				allele_group_match[f"{key}{num}"] = 1 if num_matching > 0 else 0
				protein_match[f"{key}{num}"] = 1 if num_matching > 1 else 0
				coding_match[f"{key}{num}"] = 1 if num_matching > 1 else 0
				noncoding_match[f"{key}{num}"] = 1 if num_matching > 1 else 0
			elif length_exp == 3:
				allele_group_match[f"{key}{num}"] = 1 if num_matching > 0 else 0
				protein_match[f"{key}{num}"] = 1 if num_matching > 1 else 0
				coding_match[f"{key}{num}"] = 1 if num_matching > 2 else 0
				noncoding_match[f"{key}{num}"] = 1 if num_matching > 2 else 0
			elif length_exp == 4:
				allele_group_match[f"{key}{num}"] = 1 if num_matching > 0 else 0
				protein_match[f"{key}{num}"] = 1 if num_matching > 1 else 0
				coding_match[f"{key}{num}"] = 1 if num_matching > 2 else 0
				noncoding_match[f"{key}{num}"] = 1 if num_matching > 3 else 0
			else:
				print(f"Matching sections for {key}{num}: {num_matching}/{length_exp}")

	allele_group_score = sum(allele_group_match.values()) / len(allele_group_match)
	protein_score = sum(protein_match.values()) / len(protein_match)
	coding_score = sum(coding_match.values()) / len(coding_match)
	noncoding_score = sum(noncoding_match.values()) / len(noncoding_match)

	ordered_alleles['allele_group_score'] = allele_group_score
	ordered_alleles['protein_score'] = protein_score
	ordered_alleles['coding_score'] = coding_score
	ordered_alleles['noncoding_score'] = noncoding_score

	return ordered_alleles


def filter_fasta(input_file, output_file):
	def extract_abc_alleles(header):
		"""
		Extract A, B, and C alleles from the header line.
		Returns a list of alleles or an empty list if no A, B, or C alleles are found.
		"""
		pattern = r"\s[ABC]\*\d+(?::\d+){0,3}"
		alleles = re.findall(pattern, header)
		return alleles

	total_alleles = 0
	alleles_kept = 0

	with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
		sequence = ''
		for line in infile:
			line = line.strip()
			if line.startswith('>'):
				if sequence:
					total_alleles += 1
					alleles = extract_abc_alleles(header)
					if alleles:
						alleles_kept += 1
						outfile.write('>%s\n%s\n' % (header, sequence))
				sequence = ''
				header = line[1:]
			else:
				sequence += line
		if sequence:
			alleles = extract_abc_alleles(header)
			if alleles:
				outfile.write('>%s\n%s\n' % (' '.join(alleles), sequence))

	print(f"filtered {alleles_kept}/{total_alleles} alleles")


def score_output(output_dir, alleles_path, tool="HLA_EM"):
	hla_exp_df = pd.read_csv(alleles_path)
	# TODO move Trial to the sim script
	hla_exp_df['Trial'] = hla_exp_df.index
	hla_exp_df = hla_exp_df[["A1_allele", "A2_allele", "B1_allele", "B2_allele", "C1_allele", "C2_allele", "Trial"]]

	output_dirs = sorted(os.listdir(output_dir))
	scores = []
	for trial_dir in output_dirs:
		trial_dir = os.path.join(output_dir, trial_dir)
		if not os.path.isdir(trial_dir):
			continue

		# Get the trial number
		trial_num = int(re.search(r"trial_(\d+)", trial_dir).group(1))
		trial_results = os.listdir(trial_dir)


		if tool == "HLA_EM":
			final_predictions_path = [os.path.join(trial_dir, x) for x in trial_results if x.endswith("final_predictions.csv")]
			if len(final_predictions_path) != 1:
				print(f"Error: {trial_dir} does not contain exactly one final_predictions.csv file")
				continue
			else:
				final_predictions_path = final_predictions_path[0]

			hla_obs = pd.read_csv(final_predictions_path)['HLAtype'].values.tolist()
			hla_code_obs = [re.sub(r'[()]', '', x.split()[1]) for x in hla_obs]
		elif tool == "Optitype":
			# TODO fix
			pass
		else:
			raise ValueError("Invalid tool name")

		hla_code_exp = [x[1] for x in hla_exp_df[hla_exp_df['Trial'] == trial_num].drop(columns={'Trial'}).transpose().iloc[:, 0].str.split().values.tolist()]

		scores.append(score_allele_lists(hla_code_obs, hla_code_exp))

	return pd.DataFrame(scores)


# scores = score_output("/Users/zeliason/Desktop/hla-em/output", "/Users/zeliason/Desktop/hla-em/reference/allele_record.csv")
#
# print()