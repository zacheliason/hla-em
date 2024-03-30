from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import pandas as pd
import numpy as np
import json
import re
import os


def split_into_tiers(input_value):
	if type(input_value) == str:
		result = re.match(r'^(.*?)([a-zA-Z]*)$', input_value)

		# Extracting captured groups
		hla_code, expression = [result.group(1), result.group(2)] if result else [input_value, None]

		# Splitting the string based on a regex pattern and expanding into separate parts
		tiers = re.split(r'[\*:]', hla_code, maxsplit=4)

		tiers += [''] * (5 - len(tiers)) + [expression]

		tiers = [x if x != "" else None for x in tiers]

		return tiers

	else:
		df = input_value

		df[['HLA_code', 'expression']] = df['HLA_code'].str.extract(r'^(.*?)([a-zA-Z])*$')
		expanded_columns = df['HLA_code'].str.split(r'\*|:', n=5, expand=True)
		expanded_columns.columns = ['HLAletter', 'two_digit', 'four_digit', 'six_digit', 'eight_digit']

		# Concatenate the expanded columns with the original DataFrame
		df = pd.concat([df, expanded_columns], axis=1)
		front_columns = ['HLAletter', 'two_digit', 'four_digit', 'six_digit', 'eight_digit', 'expression']
		df = df[front_columns + [x for x in df.columns if x not in front_columns]]
		df = df.where(df.notnull(), None)

		return df


def group_by_protein(df):
	df = df.sort_values(by=['MLE_Reads', 'HLA_code'], ascending=[False, True])

	grouped = df.groupby(['HLAletter', 'two_digit', 'four_digit']).agg({
		'HLAtype': 'first',
		'HLA_type': 'first',
		'HLA_code': 'first',
		'MappedReads': 'sum',
		'MappedProportion': 'sum',
		'MLE_Reads': 'sum',
		'MLE_Probability': 'sum'
	}).reset_index().sort_values(by='MLE_Reads', ascending=False)

	return grouped[['HLAtype', 'HLA_type', 'HLA_code', 'HLAletter', 'MappedReads', 'MappedProportion', 'MLE_Reads', 'MLE_Probability']]


def replace_smaller_allele(df):
	doubled_ratio_threshold = 7
	absent_ratio_threshold = 22
	absent_threshold = .18

	# Group the alleles by their "HLAletter" values
	allele_groups = df.groupby("HLAletter")

	for group_name, group in allele_groups:
		# If the group has more than one allele
		if len(group) > 1:
			group = group.sort_values(by='MLE_Probability', ascending=False)
			allele1 = group.iloc[0].copy()
			allele2 = group.iloc[1].copy()

			# Calculate the ratio between the higher and lower MLE_Probability values
			ratio = allele1['MLE_Probability'] / allele2['MLE_Probability']

			# If the ratio is greater than absent_ratio_threshold and neither allele's MLE_Probability is greater than
			# the absent_threshold, delete the smaller allele
			if ratio > absent_ratio_threshold:
				# delete smaller allele
				if allele1['MLE_Probability'] < absent_threshold and allele1['MLE_Probability'] < absent_threshold:
					df = df.drop(allele2.name)
					continue

			# If the ratio is greater than doubled_ratio_threshold, replace the smaller row with the larger one
			if ratio > doubled_ratio_threshold:
				# Adjust MLE_Probability accordingly
				allele1['MLE_Probability'] = allele1['MLE_Probability'] / 2
				df.loc[allele1.name, 'MLE_Probability'] = allele1['MLE_Probability']
				df.loc[allele2.name] = allele1

	return df


def predict_genotype_from_MLE(path, group_by_p=True):
	# Step 1: Read the TSV file into a DataFrame
	df = pd.read_csv(path, sep='\t')

	# Step 2: Extract the HLA type and MLE Probability columns
	df[['HLA_type', 'HLA_code']] = df['HLAtype'].str.split(' ', expand=True)
	df['HLA_code'] = df['HLA_code'].apply(lambda x: re.sub(r'[\(\)]', '', x))

	df = split_into_tiers(df)
	if group_by_p:
		df = group_by_protein(df)


	df = df[df['HLAletter'].isin(['A', 'B', 'C'])]

	# Step 3: Filter out alleles based on a threshold MLE_Probability
	threshold = 0.001  # Adjust the threshold as needed
	df_filtered = df[df['MLE_Probability'] >= threshold]

	# Step 4: Group alleles by their HLA type (A, B, or C)
	groups = df_filtered.groupby('HLAletter')

	# Step 5: Sort alleles within each group based on MLE Probability
	sorted_groups = {key: group.sort_values(by='MLE_Probability', ascending=False) for key, group in groups}

	# Step 6: Select the top two alleles for each HLA type
	hla_predictions = []
	for group in sorted_groups.values():
		# Select top 2 MLE_Probability values
		top_values = group['MLE_Probability'].nlargest(2).values.tolist()

		# Filter the group to only include rows equal to either of the top 2 MLE_Probability values
		filtered_group = group[group['MLE_Probability'].isin(top_values)]

		for mle_prob in filtered_group['MLE_Probability'].unique():
			mle_prob_group = filtered_group[filtered_group['MLE_Probability'] == mle_prob]
			# If only one allele is present in the group, add it to the final predictions
			if len(mle_prob_group) < 2:
				hla_predictions.append(mle_prob_group)
			# Else, sort alleles and choose top allele based on custom sorting key
			else:
				# sorting key prioritizes alleles with more matching sections relative to other alleles in mle_prop_group
				# As a secondary sorting key, the alleles are sorted alphabetically
				def custom_sort_key(code):
					matching_counts = sorted([count_common_sections(code, other_code) for other_code in mle_prob_group['HLA_code']], reverse=True)
					sum_matching = sum(matching_counts[1:])

					return (-sum_matching, code)

				# Sort the DataFrame based on the custom sorting key
				df_sorted = mle_prob_group.iloc[mle_prob_group['HLA_code'].map(custom_sort_key).argsort()]

				# Select the top allele
				hla_predictions.append(df_sorted.head(1))

	hla_prediction_df = pd.concat(hla_predictions)
	hla_prediction_df.index = hla_prediction_df['HLAtype']

	hla_prediction_df = replace_smaller_allele(hla_prediction_df)

	hla_prediction_df = hla_prediction_df.drop(columns=['HLAletter', "HLA_code", "HLA_type", "HLAtype"])

	output_dir, _ = os.path.split(path)
	hla_prediction_df.to_csv(os.path.join(output_dir, "final_predictions.csv"), index=True)

	return hla_prediction_df


def count_common_sections(pred_allele, exp_allele, include_length=False):
	"""
	Count the number of tiers are common between two alleles (excluding hla letter and expression).
	"""

	common_sections = 0

	# Split the alleles into their respective sections, excluding letter and expression
	pred_sections = split_into_tiers(pred_allele)[1:-1]
	exp_sections = split_into_tiers(exp_allele)[1:-1]

	pred_sections = [x for x in pred_sections if x is not None]
	exp_sections = [x for x in exp_sections if x is not None]

	for i in range(min(len(pred_sections), len(exp_sections))):
		if pred_sections[i] != exp_sections[i]:
			break
		common_sections += 1

	if include_length:
		return common_sections, len(pred_sections), len(exp_sections)
	else:
		return common_sections


def score_allele_lists(pred, exp, tool_typing_level=8):
	"""
	Score how well two lists of HLA alleles correspond to each other.
	"""

	# Organizes the HLA allele list in a dictionary by their first letter
	def organize_alleles_by_letter(allele_list):
		organized_alleles = {}
		for x in allele_list:
			if not pd.isna(x):
				# Check if letter (first character of allele) is in dict
				if x[0] not in organized_alleles:
					organized_alleles[x[0]] = []
				organized_alleles[x[0]].append(x)

		# Add empty lists for missing letters
		missing_letters = [item for item in ['A', 'B', 'C'] if item not in organized_alleles]
		for letter in missing_letters:
			organized_alleles[letter] = []

		return organized_alleles

	pred_dict = organize_alleles_by_letter(pred)
	exp_dict = organize_alleles_by_letter(exp)
	ordered_alleles = {}
	two_digit = {}
	four_digit = {}
	six_digit = {}
	eight_digit = {}

	# For letter in ['A', 'B', 'C']
	for hla_letter in pred_dict:
		pred_alleles_by_letter = pred_dict[hla_letter]
		exp_alleles_by_letter = exp_dict[hla_letter]

		# Assign proper hla order for comparison
		############################################

		# If no HLA alleles for a letter were predicted when some were expected
		if len(pred_alleles_by_letter) == 0 and len(exp_alleles_by_letter) > 1:
			for i, exp_allele in enumerate(exp_alleles_by_letter):
				ordered_alleles[f'{hla_letter}{i + 1}_exp'] = exp_allele

		# If one HLA allele within a letter was predicted and one was expected
		elif len(pred_alleles_by_letter) == 1 and len(exp_alleles_by_letter) == 1:
			ordered_alleles[f'{hla_letter}1_pred'], ordered_alleles[f'{hla_letter}1_exp'] = pred_alleles_by_letter[0], exp_alleles_by_letter[0]

		# If one HLA alleles within a letter was predicted and two were expected
		elif len(pred_alleles_by_letter) == 1 and len(exp_alleles_by_letter) > 1:
			# Score the two possible pairings (pred allele 1 with exp allele 1, pred allele 1 with exp allele 2, etc)
			pairing_orders = (
				count_common_sections(pred_alleles_by_letter[0], exp_alleles_by_letter[0]),
				count_common_sections(pred_alleles_by_letter[0], exp_alleles_by_letter[1]),
			)
			ordered_alleles[f'{hla_letter}1_pred'] = pred_alleles_by_letter[0]
			if pairing_orders[0] > pairing_orders[1]:
				ordered_alleles[f'{hla_letter}1_exp'] = exp_alleles_by_letter[0]
				ordered_alleles[f'{hla_letter}2_exp'] = exp_alleles_by_letter[1]
			else:
				ordered_alleles[f'{hla_letter}1_exp'] = exp_alleles_by_letter[1]
				ordered_alleles[f'{hla_letter}2_exp'] = exp_alleles_by_letter[0]

		# If two HLA alleles within a letter were predicted and one was expected
		elif len(pred_alleles_by_letter) > 1 and len(exp_alleles_by_letter) == 1:
			pairing_orders = (
				count_common_sections(pred_alleles_by_letter[0], exp_alleles_by_letter[0]),
				count_common_sections(pred_alleles_by_letter[1], exp_alleles_by_letter[0]),
			)
			ordered_alleles[f'{hla_letter}1_exp'] = exp_alleles_by_letter[0]
			if pairing_orders[0] > pairing_orders[1]:
				ordered_alleles[f'{hla_letter}1_pred'] = pred_alleles_by_letter[0]
				ordered_alleles[f'{hla_letter}2_pred'] = pred_alleles_by_letter[1]
			else:
				ordered_alleles[f'{hla_letter}1_pred'] = pred_alleles_by_letter[1]
				ordered_alleles[f'{hla_letter}2_pred'] = pred_alleles_by_letter[0]

		# If two HLA alleles within a letter were predicted and two were expected
		elif len(pred_alleles_by_letter) > 1 and len(exp_alleles_by_letter) > 1:
			pairing_orders = (
				sum([count_common_sections(pred_alleles_by_letter[0], exp_alleles_by_letter[0]),
					 count_common_sections(pred_alleles_by_letter[1], exp_alleles_by_letter[1])]),
				sum([count_common_sections(pred_alleles_by_letter[0], exp_alleles_by_letter[1]),
					 count_common_sections(pred_alleles_by_letter[1], exp_alleles_by_letter[0])])
			)

			if pairing_orders[0] > pairing_orders[1]:
				ordered_alleles[f'{hla_letter}1_pred'], ordered_alleles[f'{hla_letter}1_exp'] = pred_alleles_by_letter[0], exp_alleles_by_letter[0]
				ordered_alleles[f'{hla_letter}2_pred'], ordered_alleles[f'{hla_letter}2_exp'] = pred_alleles_by_letter[1], exp_alleles_by_letter[1]
			else:
				ordered_alleles[f'{hla_letter}1_pred'], ordered_alleles[f'{hla_letter}1_exp'] = pred_alleles_by_letter[0], exp_alleles_by_letter[1]
				ordered_alleles[f'{hla_letter}2_pred'], ordered_alleles[f'{hla_letter}2_exp'] = pred_alleles_by_letter[1], exp_alleles_by_letter[0]

		# Now that the alleles are ordered, score them
		for num in [1, 2]:
			# If no allele is predicted or expected, continue without penalizing
			if f'{hla_letter}{num}_pred' not in ordered_alleles and f'{hla_letter}{num}_exp' not in ordered_alleles:
				continue

			# If an expected allele is missing or one is predicted when not expected, penalize
			elif f'{hla_letter}{num}_pred' not in ordered_alleles or f'{hla_letter}{num}_exp' not in ordered_alleles:
				two_digit[f"{hla_letter}{num}"] = 0
				four_digit[f"{hla_letter}{num}"] = 0
				six_digit[f"{hla_letter}{num}"] = 0
				eight_digit[f"{hla_letter}{num}"] = 0
				continue

			pred_allele = ordered_alleles[f'{hla_letter}{num}_pred']
			exp_allele = ordered_alleles[f'{hla_letter}{num}_exp']

			# Count the number of matching sections between the predicted and expected alleles
			num_matching, length_pred, length_exp = count_common_sections(pred_allele, exp_allele, include_length=True)

			# Handle an expected allele with 1 section
			if length_exp == 1:
				two_digit[f"{hla_letter}{num}"] = 1 if num_matching > 0 else 0

				# Handle a tool that types to 4 or 6 or 8 digits
				if length_pred > 1 and tool_typing_level > 2:
					four_digit[f"{hla_letter}{num}"] = 0
				if length_pred > 2 and tool_typing_level > 4:
					six_digit[f"{hla_letter}{num}"] = 0
				if length_pred > 3 and tool_typing_level > 6:
					eight_digit[f"{hla_letter}{num}"] = 0

			# Handle an expected allele with 2 sections
			if length_exp == 2:
				two_digit[f"{hla_letter}{num}"] = 1 if num_matching > 0 else 0
				four_digit[f"{hla_letter}{num}"] = 1 if num_matching > 1 else 0

				# Handle a tool that types to 6 or 8 digits
				if length_pred > 2 and tool_typing_level > 4:
					six_digit[f"{hla_letter}{num}"] = 0
				if length_pred > 3 and tool_typing_level > 6:
					eight_digit[f"{hla_letter}{num}"] = 0

			# Handle an expected allele with 3 sections
			elif length_exp == 3:
				two_digit[f"{hla_letter}{num}"] = 1 if num_matching > 0 else 0
				four_digit[f"{hla_letter}{num}"] = 1 if num_matching > 1 else 0
				six_digit[f"{hla_letter}{num}"] = 1 if num_matching > 2 else 0

				# Handle a tool that types to 8 digits
				if length_pred > 3 and tool_typing_level > 6:
					eight_digit[f"{hla_letter}{num}"] = 0

			# Handle an expected allele with 4 sections
			elif length_exp == 4:
				two_digit[f"{hla_letter}{num}"] = 1 if num_matching > 0 else 0
				four_digit[f"{hla_letter}{num}"] = 1 if num_matching > 1 else 0
				six_digit[f"{hla_letter}{num}"] = 1 if num_matching > 2 else 0
				eight_digit[f"{hla_letter}{num}"] = 1 if num_matching > 3 else 0

			# Catch any unexpected lengths
			else:
				print(f"Matching sections for {hla_letter}{num}: {num_matching}/{length_exp}")
				raise ValueError("Invalid number of sections")

	two_digit_score = sum(two_digit.values()) / len(two_digit)
	four_digit_score = sum(four_digit.values()) / len(four_digit)
	six_digit_score = sum(six_digit.values()) / len(six_digit)
	eight_digit_score = sum(eight_digit.values()) / len(eight_digit)

	ordered_alleles['two_digit_score'] = two_digit_score
	ordered_alleles['four_digit_score'] = four_digit_score
	ordered_alleles['six_digit_score'] = six_digit_score
	ordered_alleles['eight_digit_score'] = eight_digit_score

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


def plot_coverage_maps(json_args_path, predicted_hla_types, covMapYmax=None):
	with open(json_args_path) as json_args_file:
		json_args = json.load(json_args_file)

	hlaRefIdCovDict = json_args['hlaRefIdCovDict']
	hlaRefID_to_seq = json_args['hlaRefID_to_seq']
	hlaRefID_to_type = json_args['hlaRefID_to_type']
	hlaRefIdGeneDict = json_args['hlaRefIdGeneDict']
	outputName = json_args['outputName']

	# Clean args
	os.remove(json_args_path)

	annotColorDict = {'E1': 'g', 'E2': 'gray', 'E3': 'y', 'E4': 'r', 'E5': 'orange', 'E6': 'b', 'E7': 'm', 'E8': 'c', 'L1': 'indigo', 'L2': 'brown'}
	annotColors = ['maroon', 'navy', 'pink', 'g', 'gray', 'k', 'y', 'r', 'orange', 'b', 'm', 'c', 'indigo']

	annotScale = 1.3

	for refId in predicted_hla_types:
		if " " in refId:
			refId = refId.split(" ")[0]

		cov = hlaRefIdCovDict[refId]
		seq_len = len(hlaRefID_to_seq[refId])
		hlaName = refId.replace(' ', '')
		hla_type = hlaRefID_to_type[refId]

		fig, ax = plt.subplots(figsize=(9, 4))
		ax.plot(np.arange(seq_len), cov, 'k', lw=0.8)
		ax.set_ylabel('Read coverage', fontsize=14, color='black')
		ax.set_title(f"{hlaName} ({hla_type})")

		if covMapYmax:
			ax.set_ylim(top=covMapYmax)

		ypos1 = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) / 12 * annotScale
		ypos2 = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) / 7.9 * annotScale
		ypos3 = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) / 5.8 * annotScale
		yposlab1 = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) / 8.5 * annotScale
		yposlab2 = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) / 6.2 * annotScale
		yposlab3 = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) / 4.8 * annotScale

		glines = []
		glabels = []
		y1end = 0
		y2end = 0
		ic = 0
		gNameLast = ''

		if refId in hlaRefIdGeneDict:
			for gene in hlaRefIdGeneDict[refId]:
				gName, gStart, gEnd = gene[0], int(gene[1]), int(gene[2])

				tname1 = gName[:2].upper()
				tname2 = gName[-2:].upper()
				if (tname1 in annotColorDict and (len(gName) < 3 or gName[2] not in '^*')):
					gc = annotColorDict[tname1]
				elif (tname2 in annotColorDict and (len(gName) < 3 or gName[-3] not in '^*')):
					gc = annotColorDict[tname2]
				else:
					if gName != gNameLast:
						ic += 1
					gc = annotColors[ic % len(annotColors)]

				if gStart >= y1end:
					ypos = ypos1
					yposlab = yposlab1
				elif gStart >= y2end:
					ypos = ypos2
					yposlab = yposlab2
				else:
					ypos = ypos3
					yposlab = yposlab3

				gline = Line2D([gStart, gEnd], [ypos, ypos], color=gc, linewidth=2)
				glines.append(ax.add_line(gline))
				glabel = Text(gStart, yposlab, gName)
				glabels.append(ax.add_artist(glabel))

				if ypos == ypos1:
					y1end = max(gEnd, glabel.get_window_extent().x1)
				elif ypos == ypos2:
					y2end = max(gEnd, glabel.get_window_extent().x1)

				gNameLast = gName

		fig.tight_layout()
		fig.savefig(f"{outputName}.{hlaName}.{hla_type}.cov.pdf", bbox_inches='tight', metadata=None)
		plt.close(fig)


def score_optitype_output(output_dir, alleles_path):
	# Load expected HLA alleles in from spreadsheet
	hla_exp_df = pd.read_csv(alleles_path)
	hla_exp_df['Trial'] = hla_exp_df.index
	hla_exp_df = hla_exp_df[["A1_allele", "A2_allele", "B1_allele", "B2_allele", "C1_allele", "C2_allele", "Trial"]]

	output_dirs = sorted(os.listdir(output_dir))
	scores = []
	# For trial in output directory
	for trial_dir_name in output_dirs:
		trial_dir = os.path.join(output_dir, trial_dir_name)
		if not os.path.isdir(trial_dir):
			continue

		# Get the trial number
		trial_num = int(re.search(r"results_run_(\d+)", trial_dir).group(1))
		optitype_dir = [os.path.join(trial_dir, x) for x in os.listdir(trial_dir) if os.path.isdir(os.path.join(trial_dir, x))][0]

		# Retrieve expected HLA alleles
		filtered_df = hla_exp_df[hla_exp_df['Trial'] == trial_num]
		dropped_df = filtered_df.drop(columns={'Trial'})
		transposed_df = dropped_df.transpose()
		transposed_df[['type', 'code']] = transposed_df[transposed_df.columns[0]].str.split(expand=True)
		transposed_df = transposed_df.dropna()
		hla_code_exp = transposed_df['code'].values.tolist()

		# Retrieve predicted HLA alleles
		try:
			trial_results_path = [os.path.join(optitype_dir, x) for x in os.listdir(optitype_dir) if x.endswith("result.tsv")][0]
		except:
			print(f"Could not find result.tsv for {trial_num}")
			continue

		hla_pred = pd.read_csv(trial_results_path, sep='\t')
		hla_pred = hla_pred[['A1', 'A2', 'B1', 'B2', 'C1', 'C2']].values.tolist()[0]

		trial_scores = score_allele_lists(hla_pred, hla_code_exp, tool_typing_level=4)
		trial_scores['Trial'] = trial_dir_name
		trial_scores['Trial_Num'] = trial_num

		scores.append(trial_scores)

	# Create, clean, and return a DataFrame from the scores
	results = pd.DataFrame(scores)
	results = results.sort_values(by='Trial_Num')
	results = results.drop(columns='Trial_Num')
	results = results.reset_index(drop=True)

	front_cols = ['Trial', 'two_digit_score', 'four_digit_score', 'six_digit_score', 'eight_digit_score']
	results = results[front_cols + [col for col in results.columns if col not in front_cols]]

	return results


def score_output(output_dir, alleles_path):
	# Load expected HLA alleles in from spreadsheet
	hla_exp_df = pd.read_csv(alleles_path)
	hla_exp_df['Trial'] = hla_exp_df.index
	hla_exp_df = hla_exp_df[["A1_allele", "A2_allele", "B1_allele", "B2_allele", "C1_allele", "C2_allele", "Trial"]]

	output_dirs = sorted(os.listdir(output_dir))
	scores = []
	# For each trial within output directory
	for trial_dir_name in output_dirs:
		trial_dir = os.path.join(output_dir, trial_dir_name)
		if not os.path.isdir(trial_dir):
			continue

		# Get the trial number
		trial_num = int(re.search(r"trial_(\d+)", trial_dir).group(1))
		trial_results = os.listdir(trial_dir)

		# Retrieve expected HLA alleles
		filtered_df = hla_exp_df[hla_exp_df['Trial'] == trial_num]
		dropped_df = filtered_df.drop(columns={'Trial'})
		transposed_df = dropped_df.transpose()
		transposed_df[['type', 'code']] = transposed_df[transposed_df.columns[0]].str.split(expand=True)
		transposed_df = transposed_df.dropna()
		hla_code_exp = transposed_df['code'].values.tolist()

		# Retrieve predicted HLA alleles
		try:
			final_predictions_path = [os.path.join(trial_dir, x) for x in trial_results if x.endswith("final_predictions.csv")][0]
		except:
			print(f"Could not find final_predictions.csv for {trial_num}")
			continue

		hla_pred = pd.read_csv(final_predictions_path)['HLAtype'].values.tolist()
		hla_code_pred = [re.sub(r'[()]', '', x.split()[1]) for x in hla_pred]

		# Score the predicted and expected HLA alleles
		trial_scores = score_allele_lists(hla_code_pred, hla_code_exp, tool_typing_level=8)
		trial_scores['Trial'] = trial_dir_name
		trial_scores['Trial_Num'] = trial_num

		scores.append(trial_scores)

	# Create, clean, and return a DataFrame from the scores
	results = pd.DataFrame(scores)
	results = results.sort_values(by='Trial_Num')
	results = results.drop(columns='Trial_Num')
	results = results.reset_index(drop=True)

	front_cols = ['Trial', 'two_digit_score', 'four_digit_score', 'six_digit_score', 'eight_digit_score']
	results = results[front_cols + [col for col in results.columns if col not in front_cols]]

	return results


# hla_em_paired_scores = score_output("/Users/zacheliason/Downloads/wetransfer_hla-thirsday_2024-03-28_2324/hla-em/output_paired", "/Users/zacheliason/Downloads/wetransfer_hla-thirsday_2024-03-28_2324/hla-em/reference/allele_record.csv")
# optitype_paired_scores = score_optitype_output("/Users/zacheliason/Downloads/wetransfer_hla-thirsday_2024-03-28_2324/hla-em/optitype_paired_output", "/Users/zacheliason/Downloads/wetransfer_hla-thirsday_2024-03-28_2324/hla-em/reference/allele_record.csv")
