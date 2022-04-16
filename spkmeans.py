import argparse
from random import seed, choice
import numpy as np
from enum import Enum

import spkmeansmodule as spk  # Our API

GOAL_SPK = 0
GOAL_WAM = 1
GOAL_DDG = 2
GOAL_LNORM  = 3
GOAL_JACOBI = 4

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
# prints the error message in case of invalid input and terminates the program
############################
def err():
    print("Invalid Input!")
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

####################
# prints indecies and
# then matrix.
# taken from hw2
####################
def print_cent(final_centroids, initial_centroid_ind):
    k = len(initial_centroid_ind)
    n = len(final_centroids)

    for i in range(k):
        ind = initial_centroid_ind[i]
        if i == k - 1:
            print(str(ind), end="\n")
        else:
            print(str(ind)+',', end="")

    # printing final centroids
    for i in range(n):
        for j in range(k):
            x = final_centroids[i][j]
            if j == k - 1:
                print(str(format(x, ".4f")), end="\n")
            else:
                print(str(format(x, ".4f")) + ',', end="")


####################
# prints matrix
####################
def print_mat(mat):
    for line in mat:
        for ind, elem in enumerate(line):
            if ind == len(line) - 1:
                print(str(format(elem, ".4f")), end="\n")
            else:
                print(str(format(elem, ".4f")) + ',', end="")

############################
# Kmeans++ initialization from HW2
############################
def delta_norm_pow2(v1, v2):
    res = 0
    for i in range(len(v1)):
        delta = v1[i] - v2[i]
        res = res + delta * delta
    return res

def dist_of_closest_cent(centroids_ind, j, vectors):
    res = float('inf')
    for cent_ind in centroids_ind:
        curr_dist = delta_norm_pow2(vectors[cent_ind], vectors[j])
        if curr_dist < res:
            res = curr_dist

    return res

def kmeanspp_algo(vectors, k):
    np.random.seed(0)
    DP = [0 for x in vectors]
    centroids_ind = []
    centroids_ind.append(np.random.choice(range(len(vectors))))

    for i in range(1,k):

        for j in range(len(vectors)):
            DP[j] = (dist_of_closest_cent(centroids_ind, j, vectors))
        sumD = sum(DP)

        for j in range(len(vectors)):
            DP[j] /= sumD

        centroids_ind.append(np.random.choice(range(len(vectors)), p=DP))

    return centroids_ind


############################
# Reading CMD argguments
############################

parser = argparse.ArgumentParser()
parser.add_argument('a', type=str, nargs='+')
args = parser.parse_args()

# Confirm we received exactly 3 arguments and first argument is an integer
if len(args.a) != 3 or not (args.a[0].isnumeric()): err()

k = int(args.a[0])
goal = Goal(args.a[1])  # Enum Goal
input_filename = args.a[2]
max_iter = 300

vectors = read_file(input_filename)

N = len(vectors)
d = len(vectors[0])

if goal == Goal.SPK:
    if k < 0 or k == 1: err()
    if k >= N: err()

    t_mat = spk.general_capi(vectors, k, N, d, GOAL_SPK)

    k = len(t_mat[0])

    initial_centroid_ind = kmeanspp_algo(t_mat, k)

    initial_centroids = [t_mat[cent_ind] for cent_ind in initial_centroid_ind]

    final_centroids = spk.fit(k, 300, initial_centroids, t_mat, k, N, 0.0)

    print_cent(final_centroids, initial_centroid_ind)

elif goal == Goal.WAM:
    w_mat = spk.general_capi(vectors, 0, N, d, GOAL_WAM)
    print_mat(w_mat)
elif goal == Goal.DDG:
    d_mat = spk.general_capi(vectors, 0, N, d, GOAL_DDG)
    print_mat(d_mat)
elif goal == Goal.LNORM:
    l_mat = spk.general_capi(vectors, 0, N, d, GOAL_LNORM)
    print_mat(l_mat)
elif goal == Goal.JACOBI:
    eigvals_and_eigvecs = spk.general_capi(vectors, 0, N, d, GOAL_JACOBI)
    for ind, elem in enumerate(eigvals_and_eigvecs[0]):
        if ind != len(eigvals_and_eigvecs[0])-1:
            s = str(format(elem, ".4f"))
            if s == "-0.0000":
                s = "0.0000"
            print(s+',', end="")
        else:
            print(str(format(elem, ".4f")), end="\n")

    print_mat(eigvals_and_eigvecs[1:])
else:
    err()
