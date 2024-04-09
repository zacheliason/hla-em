from ManipulateFiles import score_output, score_optitype_output
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as plt
import os


def percent_formatter(x, pos):
    return f'{x * 100:.0f}%'


def compare_within_tool(df, title):
    selected_columns = ['two_digit_score', 'four_digit_score', 'six_digit_score', 'eight_digit_score']

    if "Opti" in title:
        selected_columns = selected_columns[:2]

    mean_values = df[selected_columns].mean()
    columns = mean_values.index
    mean_scores = mean_values.values

    plt.figure(figsize=(10, 6))

    plt.bar(columns, mean_scores, color='xkcd:apple green')

    plt.ylabel(f'Mean Typing Accuracy')
    plt.title(f'Mean Prediction Accuracy for {title}')

    for i, v in enumerate(mean_scores):
        plt.text(i, v - 0.06, str(int(round(v, 2) * 100)) + "%", ha='center', va='bottom')

    buffer = 1.1
    plt.xlim(-buffer, len(columns) + buffer - 1)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylim(0, 1)

    plt.gca().yaxis.set_major_formatter(FuncFormatter(percent_formatter))
    plt.tight_layout()

    figpath = os.path.join(os.path.expanduser("~/Downloads"), f"{title.replace(' ', '_')}_comp.png")
    plt.savefig(figpath, dpi=600)


def compare_tools(opti_paired_scores, opti_single_scores, hla_em_single_scores, hla_em_paired_scores):
   for metric in ['two_digit_score', 'four_digit_score']:
        type = " ".join(metric.split("_")[:-1]).title()
        mean_scores = {
            'Optitype Paired': opti_paired_scores[metric].mean(),
            'Optitype Single': opti_single_scores[metric].mean(),
            'HLA-EM Paired': hla_em_paired_scores[metric].mean(),
            'HLA-EM Single': hla_em_single_scores[metric].mean()
        }

        sorted_mean_scores = sorted(mean_scores.items(), key=lambda item: item[1])
        sorted_keys, sorted_values = zip(*sorted_mean_scores)

        colors = ['cornflowerblue' if 'Opti' in key else 'xkcd:apple green' for key in sorted_keys]
        plt.figure(figsize=(10, 8))  # Adjust the width and height as needed
        bar_width = 0.4
        bar_positions = range(len(sorted_keys))

        plt.bar(bar_positions, sorted_values, color=colors, width=bar_width)
        plt.ylabel(f'Mean {type} Accuracy')
        plt.title(f'Mean HLA {type} Prediction Accuracy Across Tools')

        for i, v in enumerate(sorted_values):
            plt.text(i, v + 0.01, str(int(round(v, 2) * 100)) + "%", ha='center', va='bottom')

        plt.xticks(bar_positions, sorted_keys)

        buffer = 1.5
        plt.xlim(-buffer, len(sorted_keys) + buffer - 1)

        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)

        plt.gca().yaxis.set_major_formatter(FuncFormatter(percent_formatter))
        plt.tight_layout()

        figpath = os.path.join(os.path.expanduser("~/Downloads"), f"tool_comp_{metric}.png")
        plt.savefig(figpath, dpi=600)


opti_paired_scores = score_optitype_output("/Users/zacheliason/Downloads/hla-em/optitype_output_paired", "/Users/zacheliason/Downloads/hla-em/reference_paired/allele_record.csv")
opti_single_scores = score_optitype_output("/Users/zacheliason/Downloads/hla-em/optitype_output_single", "/Users/zacheliason/Downloads/hla-em/reference_paired/allele_record.csv")
hla_em_single_scores = score_output("/Users/zacheliason/Downloads/hla-em/output_single", "/Users/zacheliason/Downloads/hla-em/reference_paired/allele_record.csv")
hla_em_paired_scores = score_output("/Users/zacheliason/Downloads/hla-em/output_paired2", "/Users/zacheliason/Downloads/hla-em/reference_paired/allele_record.csv")

# extract num from trial
opti_paired_scores['trial'] = opti_paired_scores['Trial'].str.extract(r'(\d+)').astype(int)
opti_single_scores['trial'] = opti_single_scores['Trial'].str.extract(r'(\d+)').astype(int)
hla_em_single_scores['trial'] = hla_em_single_scores['Trial'].str.extract(r'(\d+)').astype(int)
hla_em_paired_scores['trial'] = hla_em_paired_scores['Trial'].str.extract(r'(\d+)').astype(int)

opti_paired_scores = opti_paired_scores[opti_paired_scores['trial'] < 20]
opti_single_scores = opti_single_scores[opti_single_scores['trial'] < 20]
hla_em_single_scores = hla_em_single_scores[hla_em_single_scores['trial'] < 20]
hla_em_paired_scores = hla_em_paired_scores[hla_em_paired_scores['trial'] < 20]

compare_tools(opti_paired_scores, opti_single_scores, hla_em_single_scores, hla_em_paired_scores)

compare_within_tool(hla_em_paired_scores, 'HLA-EM Paired Reads')
compare_within_tool(hla_em_single_scores, 'HLA-EM Single Reads')

compare_within_tool(opti_paired_scores, 'Optitype Paired Reads')
compare_within_tool(opti_single_scores, 'Optitype Single Reads')
