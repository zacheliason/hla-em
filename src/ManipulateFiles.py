from matplotlib import pyplot as plt
from matplotlib import font_manager
import pandas as pd
import numpy as np
import json
import re
import os


def set_font(font_name="Work Sans"):
	installDir = os.path.dirname(os.path.abspath(__file__))
	parent_dir = os.path.dirname(installDir)
	font_dir = os.path.join(parent_dir, "_".join(font_name.split()))

	font_info = {
		'dir': font_dir,
		'name': font_name
	}

	fonts = font_manager.findSystemFonts(fontpaths=font_info['dir'])
	for font in fonts:
		font_manager.fontManager.addfont(font)

	fontname = font_info['name']
	return fontname


def split_string_into_tiers(input_value):
	result = re.match(r'^(.*?)([a-zA-Z]*)$', input_value)

	# Extracting captured groups
	hla_code, expression = [result.group(1), result.group(2)] if result else [input_value, None]

	# Splitting the string based on a regex pattern and expanding into separate parts
	tiers = re.split(r'[\*:]', hla_code, maxsplit=4)

	tiers += [''] * (5 - len(tiers)) + [expression]

	tiers = [x if x != "" else None for x in tiers]

	return tiers

def split_into_tiers(df):
	df[['hla_type', 'hla_code']] = df['HLAtype'].str.split(' ', expand=True)
	df['hla_code'] = df['hla_code'].apply(lambda x: re.sub(r'[\(\)]', '', x))

	df[['hla_code', 'expression']] = df['hla_code'].str.extract(r'^(.*?)([a-zA-Z])*$')
	expanded_columns = df['hla_code'].str.split(r'\*|:', n=5, expand=True)
	expanded_columns.columns = ['hla_letter', 'two_digit', 'four_digit', 'six_digit', 'eight_digit']

	# Concatenate the expanded columns with the original DataFrame
	df = pd.concat([df, expanded_columns], axis=1)
	front_columns = ['hla_letter', 'two_digit', 'four_digit', 'six_digit', 'eight_digit', 'expression']
	df = df[front_columns + [x for x in df.columns if x not in front_columns]]
	df = df.where(df.notnull(), None)

	return df


# def plot_pie_charts(typesAll, lOrdAll, readProps, emProps, outputName):
def plot_pie_charts(predicted_hla_path, em_results_path, outname):
	def assign_color_dict(df):
		green = (143 / 255, 194 / 255, 33 / 255)
		green2 = (0.7764705882352941, 0.9686274509803922, 0.35294117647058826)
		blue = (70 / 255, 183 / 255, 232 / 255)
		blue2 = (0.4745098039215686, 0.9568627450980393, 1.0)
		yellow = (252 / 255, 222 / 255, 50 / 255)
		yellow2 = (1.0, 1.0, 0.4392156862745098)

		alleles_dict = {}
		for i, r in df[df['HLA_category'] != "Other Alleles"].iterrows():
			if r['hla_letter'] not in alleles_dict:
				alleles_dict[r['hla_letter']] = []

			alleles_dict[r['hla_letter']].append(r['HLAtype'])

		color_dict = {}
		for letter, alleles in alleles_dict.items():
			if letter == 'A':
				for i in range(len(alleles)):
					color_dict[alleles[i]] = green if i == 0 else green2
			elif letter == 'B':
				for i in range(len(alleles)):
					color_dict[alleles[i]] = yellow if i == 0 else yellow2
			elif letter == 'C':
				for i in range(len(alleles)):
					color_dict[alleles[i]] = blue if i == 0 else blue2

		return color_dict

	predicted_hla_df = pd.read_csv(predicted_hla_path)
	em_results_df = pd.read_csv(em_results_path, sep='\t')

	predicted_hla_df = split_into_tiers(predicted_hla_df)
	em_results_df = split_into_tiers(em_results_df)

	em_results_df = em_results_df.fillna(0)
	em_results_df = em_results_df.sort_values("MLE_Probability", ascending=False)
	em_results_df = group_by_protein(em_results_df)

	em_results_df['HLA_category'] = em_results_df['HLAtype']
	em_results_df.loc[~em_results_df['HLA_category'].isin(predicted_hla_df['HLAtype']), 'HLA_category'] = 'Other Alleles'

	em_results_df = em_results_df.reset_index(drop=True)
	em_results_df['pie_label_filtered_reads'] = em_results_df.apply(lambda row: '{:.1f}%'.format(row['MappedProportion'] * 100) if row['HLA_category'] != "Other Alleles" else "", axis=1)
	em_results_df['pie_label_mle'] = em_results_df.apply(lambda row: '{:.1f}%'.format(row['MLE_Probability'] * 100) if row['HLA_category'] != "Other Alleles" else "", axis=1)

	plt.rcParams['font.family'] = set_font()
	fig, axes = plt.subplots(1, 2, figsize=(10, 5), subplot_kw=dict(aspect="equal"))

	color_dict = assign_color_dict(em_results_df)

	greys = em_results_df[em_results_df['HLA_category'] == "Other Alleles"]
	greys = greys.reset_index()
	num_grey = len(greys)
	grey_linspace = np.linspace(0.35, 0.65, num_grey)
	sum_grey_mapped_proportion = '{:.1f}%'.format(greys['MappedProportion'].sum() * 100)
	sum_grey_mle_probability = '{:.1f}%'.format(greys['MLE_Probability'].sum() * 100)

	# Assign a color to each Other Allele
	for i, r in greys.iterrows():
		if i % 2 == 0:
			color_dict[r['HLAtype']] = plt.cm.Greys(.4)
		else:
			color_dict[r['HLAtype']] = plt.cm.Greys(.475)

	for ax in axes:
		ax.set_prop_cycle('color', [color_dict[x] for x in em_results_df['HLAtype']])

	plot_pie(axes[0], em_results_df['MappedReads'], em_results_df['pie_label_filtered_reads'], sum_grey_mapped_proportion)
	plot_pie(axes[1], em_results_df['MLE_Probability'], em_results_df['pie_label_mle'], sum_grey_mle_probability)

	axes[1].legend(labels=em_results_df['HLA_category'])
	legend = axes[1].get_legend()
	handles = legend.legendHandles
	labels = legend.texts

	# Remove duplicate legend entries
	seen = []
	new_labels = []
	new_handles = []
	for i, label in enumerate(labels):
		if label.get_text() not in seen:
			seen.append(label.get_text())
			new_labels.append(label.get_text())
			new_handles.append(handles[i])

	# Reset legend
	axes[1].legend(new_handles, new_labels, title="HLA Types", loc="center left", frameon=False, bbox_to_anchor=(1, 0, 0.5, 1))

	# Set titles for subplots
	axes[0].set_title("Mapped Read Proportions")
	axes[1].set_title("Maximum Likelihood Estimate")

	# Adjust layout and save the figure
	fig.subplots_adjust(wspace=0, right=0.8)
	# fig.tight_layout(rect=[0, 0, 0.9, 0.9])
	fig.tight_layout()
	fig.savefig(outname + '.props.pdf')
	plt.close(fig)


def plot_pie(ax, props, labels, other_val=None):
	wedges, texts = ax.pie(props)
	other_alleles = []
	for i, p in enumerate(wedges):
		if labels[i] != "":
			ang = (p.theta2 - p.theta1) / 2. + p.theta1
			y = 0.8 * np.sin(np.deg2rad(ang))
			x = 0.8 * np.cos(np.deg2rad(ang))
			ax.annotate(labels[i], xy=(x, y), ha='center', va='center')
		else:
			other_alleles.append(p)

	if other_val is not None:
		ang = (other_alleles[-1].theta2 - other_alleles[0].theta1) / 2. + other_alleles[0].theta1
		y = 0.7 * np.sin(np.deg2rad(ang))
		x = 0.7 * np.cos(np.deg2rad(ang))
		ax.annotate(other_val, xy=(x, y), ha='center', va='center')


def group_by_protein(df):
	df = df.sort_values(by=['MLE_Reads', 'hla_code'], ascending=[False, True])

	grouped = df.groupby(['hla_letter', 'two_digit', 'four_digit']).agg({
		'HLAtype': 'first',
		'hla_type': 'first',
		'hla_code': 'first',
		'MappedReads': 'sum',
		'MappedProportion': 'sum',
		'MLE_Reads': 'sum',
		'MLE_Probability': 'sum'
	}).reset_index().sort_values(by='MLE_Reads', ascending=False)

	return grouped[['HLAtype', 'hla_type', 'hla_code', 'hla_letter', 'MappedReads', 'MappedProportion', 'MLE_Reads', 'MLE_Probability']]


def replace_smaller_allele(df, trial_name=None, training_path=""):
	# Group the alleles by their "hla_letter" values
	allele_groups = df.groupby("hla_letter")

	for group_name, group in allele_groups:
		# If the group has more than one allele
		if len(group) > 1:
			group = group.sort_values(by='MLE_Probability', ascending=False)
			allele1 = group.iloc[0].copy()
			allele2 = group.iloc[1].copy()

			if training_path != "":
				with open(training_path, 'a') as f:
					f.write("\t".join([trial_name, group_name, str(allele1['MLE_Probability']), str(allele2['MLE_Probability'])]) + "\n")

			ml1 = allele1['MLE_Probability']
			ml2 = allele2['MLE_Probability']
			ratio = ml1 / ml2

			# Choose whether to keep, drop, or copy the smaller allele from the larger based on the pair's relative MLE probabilities
			action = predict_allele_action(ml1, ml2, ratio)

			if action == "DROP":
				df = df.drop(allele2.name)
			elif action == "COPY":
				allele1[['MLE_Probability', 'MappedProportion']] /= 2
				allele1[['MLE_Reads', 'MappedReads']] //= 2

				columns_to_update = ['MLE_Probability', 'MappedProportion', 'MLE_Reads', 'MappedReads']
				# print()
				# print(allele1.index)
				# print(allele1.name)
				# print(df.columns)
				# print(df.index)
				# print()
				df.loc[allele1.name, columns_to_update] = allele1[columns_to_update]

				df.loc[allele2.name] = allele1

	return df


def predict_genotype_from_MLE(em_results_path, outname, trial_name, training_csv, group_by_p=True):
	if training_csv != "":
		if not os.path.exists(training_csv):
			with open(training_csv, 'w+') as f:
				f.write("trial_name\thla_letter\tmle1\tmle2\n")
		else:
			sep = "," if training_csv.endswith(".csv") else "\t"
			training_df = pd.read_csv(training_csv, sep=sep).drop_duplicates()
			training_df = training_df[training_df['trial_name'] != trial_name]
			training_df.to_csv(training_csv, index=False, sep=sep)

	# Step 1: Read the TSV file into a DataFrame
	df = pd.read_csv(em_results_path, sep='\t')

	df = split_into_tiers(df)
	if group_by_p:
		old_df = df.copy()
		df = group_by_protein(df)

	df = df[df['hla_letter'].isin(['A', 'B', 'C'])]

	# Step 3: Filter out alleles based on a threshold MLE_Probability
	threshold = 0.001  # Adjust the threshold as needed
	df_filtered = df[df['MLE_Probability'] >= threshold]

	# Step 4: Group alleles by their HLA type (A, B, or C)
	groups = df_filtered.groupby('hla_letter')

	# Step 5: Sort alleles within each group based on MLE Probability
	sorted_groups = {key: group.sort_values(by='MLE_Probability', ascending=False) for key, group in groups}

	# Step 6: Select the top two alleles for each HLA type
	hla_predictions = []
	for group in sorted_groups.values():
		# Select top 2 MLE_Probability values
		top_values = group['MLE_Probability'].drop_duplicates().nlargest(2).values.tolist()

		# Filter the group to only include rows equal to either of the top 2 MLE_Probability values
		filtered_group = group[group['MLE_Probability'].isin(top_values)]

		for i, mle_prob in enumerate(filtered_group['MLE_Probability'].unique()):
			mle_prob_group = filtered_group[filtered_group['MLE_Probability'] == mle_prob]
			# If only one allele is present in the group, add it to the final predictions
			if len(mle_prob_group) < 2:
				hla_predictions.append(mle_prob_group)
			# Else, sort alleles and choose top allele based on custom sorting key
			else:
				# sorting key prioritizes alleles with more matching sections relative to other alleles in mle_prop_group
				# As a secondary sorting key, the alleles are sorted alphabetically
				def custom_sort_key(code):
					matching_counts = sorted([count_common_sections(code, other_code) for other_code in mle_prob_group['hla_code']], reverse=True)
					sum_matching = sum(matching_counts[1:])

					return (-sum_matching, code)

				# Sort the DataFrame based on the custom sorting key
				df_sorted = mle_prob_group.iloc[mle_prob_group['hla_code'].map(custom_sort_key).argsort()]

				# Select the top allele
				# If i == 0, there are multiple alleles with the highest MLE_Probability value so grab the top 2 and break
				if i == 0:
					hla_predictions.append(df_sorted.head(2))
					break
				# If i == 1, there is only one allele left to add
				else:
					hla_predictions.append(df_sorted.head(1))

	hla_prediction_df = pd.concat(hla_predictions)
	hla_prediction_df.index = hla_prediction_df['HLAtype']

	hla_prediction_df = replace_smaller_allele(hla_prediction_df, trial_name=trial_name, training_path=training_csv)

	hla_prediction_df = hla_prediction_df.drop(columns=['hla_letter', "hla_code", "hla_type", "HLAtype"])

	hla_prediction_df.to_csv(outname + ".final_predictions.csv", index=True)

	return hla_prediction_df


def count_common_sections(pred_allele, exp_allele, include_length=False):
	"""
	Count the number of tiers are common between two alleles (excluding hla letter and expression).
	"""

	common_sections = 0

	# Split the alleles into their respective sections, excluding letter and expression
	pred_sections = split_string_into_tiers(pred_allele)[1:-1]
	exp_sections = split_string_into_tiers(exp_allele)[1:-1]

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
	hlaRefID_to_len = json_args['hlaRefID_to_length']
	hlaRefID_to_type = json_args['hlaRefID_to_type']
	outputName = json_args['outputName']

	for refId in predicted_hla_types:
		if " " in refId:
			refId = refId.split(" ")[0]

		cov = hlaRefIdCovDict[refId]
		seq_len = hlaRefID_to_len[refId]
		hlaName = refId.replace(' ', '')
		hla_type = hlaRefID_to_type[refId]

		fig, ax = plt.subplots(figsize=(9, 4))
		plt.rcParams['font.family'] = set_font()

		ax.bar(np.arange(seq_len), cov, width=1.0, align='edge')
		ax.set_ylabel('Read Coverage', fontsize=14, color='black')
		ax.set_title(f"{hlaName} ({hla_type})")

		if covMapYmax:
			ax.set_ylim(top=covMapYmax)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		fig.tight_layout()
		fig.savefig(f"{outputName}.{hlaName}.{hla_type.replace(':', '.')}.cov.pdf", bbox_inches='tight', metadata=None)
		plt.close(fig)

	# Clean args
	os.remove(json_args_path)


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


#
# def predict_allele_action(mle1, mle2, ratio):
# 	if mle2 <= 0.0345:
# 		if mle1 <= 0.188:
# 			return "DROP"
# 		else:  # if mle1 > 0.188
# 			return "COPY"
# 	else:  # if mle2 > 0.0345
# 		if mle2 <= 0.0519:
# 			if ratio <= 3:
# 				return "KEEP"
# 			else:  # if mle2 > 0.0492
# 				return "COPY"
# 		else:  # if mle2 > 0.0519
# 			return "KEEP"


def predict_allele_action(mle1, mle2, ratio):
	if mle2 <= 0.0345:
		if mle1 <= 0.188:
			return "DROP"
		else:  # if mle1 > 0.188
			return "COPY"
	else:  # if mle2 > 0.0345
		if mle2 <= 0.0519:
			if ratio <= 3:
				return "KEEP"
			else:  # if mle2 > 0.0492
				return "COPY"
		else:  # if mle2 > 0.0519
			return "KEEP"


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


# optitype_paired_scores = score_optitype_output("/Users/zacheliason/Downloads/wetransfer_hla-thirsday_2024-03-28_2324/hla-em/optitype_paired_output", "/Users/zacheliason/Downloads/wetransfer_hla-thirsday_2024-03-28_2324/hla-em/reference/allele_record.csv")


# hla_em_paired_scores = score_output("/Users/zeliason/Desktop/hla-em/output_paired", "/Users/zeliason/Desktop/hla-em/reference/allele_record.csv")
# print(hla_em_paired_scores['two_digit_score'].mean())
# print(hla_em_paired_scores['four_digit_score'].mean())
# print(hla_em_paired_scores['six_digit_score'].mean())
# print(hla_em_paired_scores['eight_digit_score'].mean())
#
# print()