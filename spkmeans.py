import argparse
from random import seed
import numpy as np
from enum import Enum

import spkmeansmodule as spk  # Our API

seed(0)


# Creating Enum of Goal
class Goal(Enum):
    SPK = 'spk'
    WAM = 'wam'
    DDG = 'ddg'
    LNORM = 'lnorm'
    JACOBI = 'jacobi'
    MISSING = 'missing'

    @classmethod
    def _missing_(cls, value):
        return Goal.MISSING


############################
# input: error message to print
# prints the error message and terminates the program
############################
def term(string_to_print):
    print(string_to_print)
    exit(1)


############################
# input: file_name as string
# returns list of vectors as floats
############################
def read_file(input_filename):
    file = open(input_filename)
    lines = file.readlines()
    vectors = []
    for line in lines:
        vec_strings = line.split(',')
        vec_floats = [float(s) for s in vec_strings]
        vectors.append(vec_floats)
    file.close()
    return vectors


############################
# Reading CMD argguments
############################

parser = argparse.ArgumentParser()
parser.add_argument('a', type=str, nargs='+')
args = parser.parse_args()

# Confirm we received exactly 3 arguments and first argument is an integer
if len(args.a) != 3 or not (args.a[0].isnumeric()): term("Invalid Input!")

k = int(args.a[0])
goal = Goal(args.a[1])  # Enum Goal
input_filename = args.a[2]
max_iter = 300

vectors = read_file(input_filename)
N = len(vectors)
d = len(vectors[0])

if goal == Goal.SPK:
    if k <= 1: term("Invalid Input!")
    if k >= N: term("Invalid Input!")

    ## Call my_spk from our module spkmeans
    ## my_spk: input: vectors, k, len(vectors) = N, len(vectors[0]) = d, int goal
    ##         return: data_points
    ## cluster data_points in to k clusters:
    ## k_means++ init from hw2
    ## k_means from hw1

elif goal == Goal.WAM:
    pass
elif goal == Goal.DDG:
    pass
elif goal == Goal.LNORM:
    pass
elif goal == Goal.JACOBI:
    eigvals_and_eigvecs = spk.general_capi(vectors, 0, N, d, "jacobi")
    print(eigvals_and_eigvecs[0])
    print(eigvals_and_eigvecs[1:])
else:
    term("Invalid Input!")

################################
## FROM HW2

#
# def print_cent(final_centroids, initial_centroind_ind):
#     # printing initial centroids indices
#     n = len(initial_centroind_ind)
#     for i in range(n):
#         ind = initial_centroind_ind[i]
#         if i == n - 1:
#             print(ind)
#         else:
#             print(ind, end=", ")
#
#     # printing final centroids
#     n = len(final_centroids[0])
#     for i in range(k):
#         for j in range(n):
#             x = final_centroids[i][j]
#             if j == n - 1:
#                 print('%.4f' % x)
#             else:
#                 print('%.4f' % x, end=", ")
