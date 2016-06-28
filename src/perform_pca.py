#!/nfs/goldstein/software/python2.7/bin/python2.7
"""
Perform PCA, retain the top two axes of variation, and plot cases and controls
to infer if the cases and controls are sufficiently matched to proceed with
further analyses
"""

import argparse
import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
from copy import deepcopy
import matplotlib as mpl
mpl.use("pdf")
import matplotlib.pyplot as plt

AFFECTATION_COLUMN = 5

def create_genotypes_list(ped_genotypes, arg, num_variants, line_num):
    """return a list of the genotypes for the current sample or raise an
    argparse Exception if the format is invalid
    """
    if num_variants * 2 != len(ped_genotypes):
        raise argparse.ArgumentTypeError(
            "{} has an incorrect number of columns in row {}: {}".format(
                arg, line_num, len(ped_genotypes) + 6))
    genotypes_list = []
    try:
        ped_genotypes = [int(genotype) for genotype in ped_genotypes]
    except ValueError:
        raise argparse.ArgumentTypeError(
            "{} has a non-integer encoded genotype in row {}".format(
                arg, line_num))
    for variant_num in xrange(num_variants):
        a1_index = variant_num * 2
        a2_index = a1_index + 1
        if (ped_genotypes[a1_index] == 0 and
            ped_genotypes[a2_index] == 0):
            genotypes_list.append(-1)
        elif (ped_genotypes[a1_index] <= 0 or
              ped_genotypes[a1_index] > 2 or
              ped_genotypes[a2_index] <= 0 or
              ped_genotypes[a2_index] > 2):
            raise argparse.ArgumentTypeError(
                "{} has an improper genotype endoding at genotype #{}: {} "
                "in row {}".format(
                    arg, variant_num + 1, "{} {}".format(
                        *ped_genotypes[a1_index:a1_index + 2]), line_num))
        elif (ped_genotypes[a1_index] == 1 and
              ped_genotypes[a2_index] == 1):
            genotypes_list.append(0)
        elif (ped_genotypes[a1_index] == 2 and
              ped_genotypes[a2_index] == 2):
            genotypes_list.append(2)
        else:
            genotypes_list.append(1)
    return genotypes_list

def PED_file(arg):
    """return a list of affectation statuses and the two-dimensional
    numpy array of genotypes by sample
    """
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(arg + " does not exist")
    affectation_statuses = []
    genotypes_matrix = []
    with open(arg) as ped:
        first_sample = ped.next().split()
        num_fields = len(first_sample)
        if num_fields < 8 or num_fields % 2:
            raise argparse.ArgumentTypeError(
                "{} has an incorrect number of columns in row 1: {}".format(
                    arg, num_fields))
        affectation_statuses.append(int(first_sample[AFFECTATION_COLUMN]) - 1)
        num_variants = (num_fields - 6) / 2
        genotypes_matrix.append(create_genotypes_list(
            first_sample[6:], arg, num_variants, 1))
        for line_num, line in enumerate(ped):
            ped_fields = line.split()
            if len(ped_fields) < 8:
                raise argparse.ArgumentTypeError(
                    "{} has an incorrect number of columns in row {}: {}".format(
                        arg, line_num + 2, len(ped_fields)))
            affectation_statuses.append(int(ped_fields[AFFECTATION_COLUMN]) - 1)
            genotypes_matrix.append(create_genotypes_list(
                ped_fields[6:], arg, num_variants, line_num + 2))
    return (np.array(affectation_statuses),
            np.array(genotypes_matrix, dtype=np.int8), arg)

def main(affectation_statuses, genotypes_matrix, ped_path):
    # impute missing values
    g = Imputer(missing_values=-1, strategy="most_frequent", axis=0)
    g.fit(genotypes_matrix)
    genotypes = g.transform(genotypes_matrix) 
    pca = PCA(n_components=2)
    del genotypes_matrix
    pcs = pca.fit_transform(genotypes)
    fig = plt.figure()
    ctrls = plt.scatter(pcs[affectation_statuses == 0, 0],
                        pcs[affectation_statuses == 0, 1],
                        c="b", s=4)
    cases = plt.scatter(pcs[affectation_statues == 1, 0],
                        pcs[affectation_statuses == 1, 1],
                        c="y", s=4)
    plt.title("PC1 vs. PC2 post-pruning")
    plt.ylabel("PC2")
    plt.xlabel("PC1")
    plt.legend((ctrls, cases), ("Cases", "Controls"), scatterpoints=1)
    fig.savefig("{ped_path}.pdf".format(ped_path=ped_path))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("PED", type=PED_file, help="the PED file to operate on")
    args = parser.parse_args()
    main(*args.PED)
