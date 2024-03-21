import pandas as pd
import numpy as np
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