#!rusr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
import math
import sys
import os

matplotlib.use('Agg')

# Difference value between successive loglikelihoods below which algorithm is deemed to have converged
conVal = 1e-4
numIter = 5


def parse_read_info(line_data):
    line_data = line_data
    read_info = line_data.reshape(-1, 4)
    read_nums = read_info[:, 0].astype(np.int32)
    lm = read_info[:, 1].astype(np.float32)
    le = read_info[:, 2].astype(np.float32)
    # genes = read_info[:, 3] #line_data[1:][::4][read_nums == 0]
    return read_nums, lm, le # genes

def read_table(readsTable):
    arrays_list = [np.array(line.split('\t')) for line in readsTable[2:]]
    list_len = arrays_list[0].shape[0]
    arrays_list = [np.array(arr) for arr in arrays_list if arr.shape[0] == list_len]
    lines = np.vstack(arrays_list)
    hla_types = lines[:, 0]
    lines = lines[:, 1:]
    parsed_lines = [parse_read_info(line) for line in lines]

    parsed_data = np.array(parsed_lines)

    read_hits_matrix = parsed_data[:, 0]
    lm_matrix = parsed_data[:, 1]
    le_matrix = parsed_data[:, 2]

    uniq_reads, *is_read_ambig = readsTable[1].split('\t')

    return hla_types, read_hits_matrix, lm_matrix, le_matrix, is_read_ambig, int(uniq_reads)


def em_component(read_hits_matrix, lm_matrix, le_matrix, max_steps=500, conv_val=1e-4, num_iter=5):
    k = read_hits_matrix.shape[0]  # number of HLA types
    m = read_hits_matrix.shape[1]  # number of reads
    iter_out = None
    converged_once = False
    l_out = -np.inf
    err_out = 0.0
    phi_out = np.zeros(k)
    steps_out = 0

    for iteration_number in range(num_iter):
        converged = False

        if iteration_number == 0:
            err = 0.005
            phi = np.full(k, 1.0 / k)
        else:
            err = 0.05 * np.random.random()
            phi = np.random.random(k)
            phi /= phi.sum()

        l = -np.inf
        w = np.zeros((m, k))
        steps = 0

        while not converged:
            steps += 1
            if steps > max_steps:
                if not converged_once:
                    max_steps += 100
                break

            # E step
            for j in range(k):
                lm = lm_matrix[j]
                le = le_matrix[j]
                read_hits = read_hits_matrix[j]
                column_values = ((1.0 - err) ** lm * err ** le) * phi[j]
                masked_column_values = read_hits * column_values
                w[:, j] = masked_column_values

            w /= w.sum(axis=1, keepdims=True)

            # M step
            ## err
            b_num = np.sum(w[:, :] * le_matrix.T)
            b_den = np.sum(w[:, :] * lm_matrix.T)
            b = b_num / b_den
            err = b / (1.0 + b)

            ## phi
            phi[:] = w.sum(axis=0) / m

            # Calculate log-likelihood
            l_new = 0.0
            for j in range(k):
                lm = lm_matrix[j]
                le = le_matrix[j]
                read_hits = read_hits_matrix[j]

                numerator = (1.0 - err) ** lm * err ** le * phi[j]
                numerator[numerator == 0] = sys.float_info.min

                denominator = w[:, j]
                denominator[denominator == 0] = sys.float_info.min

                log_numerator = np.log(numerator)
                log_denominator = np.log(denominator)

                log_likelihood_values = w[:, j] * (log_numerator - log_denominator)

                masked_log_values = np.where(read_hits > 0, log_likelihood_values, 0)

                l_new += masked_log_values.sum()

            print(f"\titeration: {iteration_number + 1}/{num_iter}, convergence value: {round(l_new - l, 6)}, solution found: {converged_once}", end='\r')

            if l_new - l < conv_val:
                converged = True
                if iteration_number == 0 or l_new > l_out or iter_out is None:
                    converged_once = True
                    l_out = l_new
                    err_out = err
                    phi_out = phi
                    steps_out = steps
                    iter_out = iteration_number

            l = l_new

    if not converged_once:
        raise RuntimeError(f'EM algorithm failed to converge after {max_steps} steps and {num_iter} iterations; aborting.')

    return err_out, phi_out, steps_out


def EmAlgo(readsTable, outname, thresholdTpm=1.5):
    hla_types, read_hits_matrix, lm_matrix, le_matrix, is_read_ambig, uniq_reads = read_table(readsTable)
    num_total_reads = read_hits_matrix.sum()
    err_out, phi_out, steps_out = em_component(read_hits_matrix, lm_matrix, le_matrix)

    mle_reads = uniq_reads * phi_out
    tpm_values = mle_reads * 1e6 / num_total_reads

    mle_reads = np.round(mle_reads).astype(int)

    # Get number of reads that pass TPM threshold:
    passes_threshold = tpm_values > thresholdTpm
    totalOutReads = np.sum(mle_reads[passes_threshold])
    output_lines = []
    for j, hla in enumerate(hla_types):
        if tpm_values[j] > thresholdTpm:
            mapped_reads = read_hits_matrix[j].sum()

            HLA_reference = {
                'HLAtype': str(hla),
                'MappedReads': mapped_reads,
                'MappedProportion': mapped_reads / num_total_reads,
                'MLE_Reads': mle_reads[j],
                'MLE_Probability': mle_reads[j] / totalOutReads
            }
            output_lines.append(HLA_reference)

    df = pd.DataFrame(output_lines)
    df = df.sort_values('MLE_Probability', ascending=False)
    df.to_csv(outname + ".results.tsv", sep='\t', index=False)

    with open(f"{outname}.em_algorithm.Log", "w") as f:
        f.write(f"Converged to < {conVal:.1e} in {steps_out} iterations\n")
        f.write(f"err\t{err_out:.5f}\n")

    return df


class mappedRead:
    def __init__(self, inputList):
        self.hlaType = inputList[0]
        self.readInfo = []
        self.readNum = 0
        self.genes = []
        count = 0

        for val in inputList[1:]:
            if (count % 4) == 0:
                val = int(val)
                self.readInfo.append([val])
                self.readNum += val
            elif (count % 4) == 3:
                if not val:
                    self.genes.append([])
                else:
                    self.genes.append(val.split(','))
            else:
                self.readInfo[-1].append(float(val))
            count += 1


def main(argv):
    readTableFile = sys.argv[1]
    readTable = []
    if len(sys.argv) > 3:
        outputName = sys.argv[3]
    else:
        outputName = 'hlaType'
    with open(readTableFile, 'r') as inFile:
        for line in inFile:
            readTable.append(line.strip('\n'))
    # EmAlgo(readTable, allReadsNum=int(sys.argv[2]), thresholdTpm=1.48, outputName=outputName, printResult=True)
    EmAlgo(readTable, outname=outputName, thresholdTpm=1.48)


if __name__=="__main__":
    main(sys.argv)
