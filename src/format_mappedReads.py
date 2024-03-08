import re
import os
import pandas as pd
import numpy as np
from collections import Counter

def unique(list1):
    x = np.array(list1)
    return list(np.unique(x))

path = "/Users/zacheliason/Documents/Work/zhang/2024/new/output/output.readsTable.tsv"

lines = []
with open(path) as table:
    i = 0
    for line in table:
        line = line.strip()
        if i == 0:
            columns = line.split("\t")
            asdf1 = columns[-5:]
            print('okay')
            columns_t = [x + "\t\t\t\t" for x in columns]
            columns_t = ["Gene\t"] + columns_t
            columns = "".join(columns_t)

            # columns_t = columns.split('\t')[:-1]
            # columns = "\t".join(columns_t)
        elif i == 1:
            ambig = line.split("\t")
            ambig_t = [x + "\t\t\t\t" for x in ambig]
            ambig_t[0] = ambig_t[0].replace("\t", "")
            ambig_t[0] += "\t"
            ambig = "".join(ambig_t)
        else:
            lines.append(line)
        i += 1

with open(path.replace(".tsv", "_formatted.tsv"), "w") as table:
    columns = "\t".join(columns.split("\t")[:-2])
    ambig = "\t".join(ambig.split("\t")[:-2])

    table.write(columns + "\n")
    table.write(ambig + "\n")

    c = columns.split("\t")
    a = ambig.split("\t")

    for line in lines:
        x = line.split("\t")
        if len(x) != len(c) or len(x) != len(a):
            print('okay')
        else:
            table.write(line + "\n")

interest = ["HLA:HLA18924 (A*02:736)", "HLA:HLA00005 (A*02:01:01:01)"]

df = pd.read_csv(path.replace(".tsv", "_formatted.tsv"), sep="\t")
dff = df[df['Gene'].str.contains("A*02:")]

lost = []
cs = []
for c in dff.columns.tolist()[1:]:
    if "Unnamed" in c:
        cs.append(c)
    else:
        s = pd.to_numeric(dff[c]).sum()
        if s != 0:
            cs.append(c)
        else:
            lost.append(c)


# unique_lost = []
# for l in lost:
#     m = re.search(r"HLA:HLA(\d+).*", l)
#     if m:
#         unique_lost.append(m.group(1))
#
# c = Counter(unique_lost)

dff = dff[cs]
# dff = dff[dff.columns.tolist()[:3]][dff[dff.columns.tolist()[:3]].columns.tolist()[-1]].values.tolist()
matched = df.columns.tolist()
matched = [x for x in matched if not x.startswith("Un")]
df = df[df["Gene"].str.contains("HLA")]
df = df[matched]
df.index = df['Gene']
df = df.drop('Gene', axis=1)
# df = df.filter(items=interest, axis=0)
index = df.index.tolist()
sums = np.sum(df.values, dtype=np.float64, axis=1)

sums, index = zip(*sorted(zip(sums, index), reverse=True))

print("asdf")
# dataFrame = df.sum(axis=1)
