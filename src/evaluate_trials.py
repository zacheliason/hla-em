from ManipulateFiles import count_common_sections
import argparse as argp
import pandas as pd
import sys
import os
import re

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
			final_predictions = [os.path.join(trial_dir, x) for x in trial_results if x.endswith("final_predictions.csv")]
			# sort by length
			final_predictions = sorted(final_predictions, key=len)
			# choose the longest
			final_predictions_path = final_predictions[-1]
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



def evaluate(name, output_dir, reference_allele_record_path, outname):
	if 'optitype' in output_dir.lower() or 'optitype' in name.lower():
		scores = score_optitype_output(output_dir, reference_allele_record_path)
	else:
		scores = score_output(output_dir, reference_allele_record_path)

	two_digit_score = round(scores['two_digit_score'].mean(), 3)
	four_digit_score = round(scores['four_digit_score'].mean(), 3)
	six_digit_score = round(scores['six_digit_score'].mean(), 3)
	eight_digit_score = round(scores['eight_digit_score'].mean(), 3)

	print(f"{name} scores")
	print("-" * len(f"{name} scores"))
	print(f"two-digit score: {two_digit_score:.2%}")
	print(f"four-digit score: {four_digit_score:.2%}")
	print(f"six-digit score: {six_digit_score:.2%}")
	print(f"eight-digit score: {eight_digit_score:.2%}")
	print()

	print(f"saved scores to {outname}")
	print()

	scores.to_csv(outname, index=False)


if __name__ == "__main__":
	# add argp
	parser = argp.ArgumentParser(description='Evaluate HLA typing results')
	parser.add_argument('-n', '--name', help='Name of the tool')
	parser.add_argument('-o', '--output_dir', help='Path to the output directory')
	parser.add_argument('-r', '--reference_allele_record_path', help='Path to the reference allele record')
	parser.add_argument('--outname', help='path of the scores csv to save')

	args = parser.parse_args()

	# if missing any quit
	if not args.name or not args.output_dir or not args.reference_allele_record_path or not args.outname:
		parser.print_help()
		sys.exit(1)

	evaluate(args.name, args.output_dir, args.reference_allele_record_path, args.outname)
