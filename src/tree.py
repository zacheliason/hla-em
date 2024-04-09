from sklearn.tree import DecisionTreeClassifier, _tree
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import pandas as pd
import numpy as np


def add_labels(filepath, training_path):
    df = pd.read_csv(filepath)
    # extract trial num from reference fasta column
    df['trial'] = df['reference_fasta'].str.extract(r'(trial_\d+)')

    labels_dict = []
    for i, r in df.iterrows():
        trial_name = r['trial']
        for letter in ['A', 'B', 'C']:
            # extract allele from reference fasta column
            if r[f'{letter}{1}_allele'] == r[f'{letter}{2}_allele']:
                label = "COPY"
            elif pd.isna(r[f'{letter}{1}_allele']) or pd.isna(r[f'{letter}{2}_allele']):
                label = "DROP"
            else:
                label = "KEEP"

            labels_dict.append({
                'trial_name': trial_name,
                'hla_letter': letter,
                'label': label
            })


    labels_df = pd.DataFrame(labels_dict)

    sep = "," if training_path.endswith('.csv') else "\t"
    training_df = pd.read_csv(training_path, sep=sep)

    # merge on trial_name and letter
    training_df = training_df.merge(labels_df, on=['trial_name', 'hla_letter'])
    training_df = training_df.sort_values(['trial_name', 'hla_letter'])
    training_df.to_csv(training_path, sep=sep, index=False)

    print()



def tree_to_code(tree, feature_names):
    index_to_class = {i: f'"{class_name}"' for i, class_name in enumerate(tree.classes_)}
    tree_ = tree.tree_
    feature_name = [
        feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
        for i in tree_.feature
    ]

    code_string = "\n# COPY THIS CODE INTO YOUR SCRIPT\n"
    code_string += "# Choose whether to keep, drop, or copy the smaller allele from the larger based on the pair's relative MLE probabilities\n"
    code_string += f"def predict_allele_action({', '.join(feature_names)}):\n"

    def recurse(node, depth, code_string):
        indent = "\t" * depth
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            name = feature_name[node]
            threshold = tree_.threshold[node]
            code_string += f"{indent}if {name} <= {np.round(threshold,4)}:\n"
            code_string = recurse(tree_.children_left[node], depth + 1, code_string)
            code_string += f"{indent}else:  # if {name} > {np.round(threshold,4)}\n"
            code_string = recurse(tree_.children_right[node], depth + 1, code_string)
        else:
            node_value = tree_.value[node]
            value_index = index_to_class[np.argmax(node_value)]

            code_string += f"{indent}return {value_index}\n"

        return code_string

    code_string = recurse(0, 1, code_string)
    return code_string


def train_tree_model(training_data_path):
    # Load and preprocess data
    sep = "," if training_data_path.endswith('.csv') else "\t"
    data = pd.read_csv(training_data_path, sep=sep)
    data['ratio'] = data['mle1'] / data['mle2']
    feature_names = ['mle1', 'mle2', 'ratio']

    # Split data into features and labels
    X = data[feature_names]
    y = data['label']

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    clf = DecisionTreeClassifier(max_depth=4)
    # clf = DecisionTreeClassifier()
    # clf.fit(X_train, y_train)
    clf.fit(X, y)

    y_pred = clf.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy:", accuracy)

    python_tree_code = tree_to_code(clf, feature_names=feature_names)
    return python_tree_code

# training_path = '/Users/zacheliason/Downloads/hla-em/output_training/training.tsv'
# training_path = '/Users/zacheliason/Downloads/hla-em/output_paired/training.tsv'
training_path = '/Users/zacheliason/Downloads/hla-em/output_training/training.tsv'

# add_labels('/Users/zacheliason/Downloads/hla-em/reference_paired/allele_record.csv', training_path)
python_code = train_tree_model(training_path)
print(python_code)
print()
