import itertools
import random
import sys
import time
import math
import numpy as np


def held_karp(task_start, task_end, velocity):
    """
    Implementation of Held-Karp, an algorithm that solves the Traveling
    Salesman Problem using dynamic programming with memoization.

    Parameters:
        dists: distance matrix

    Returns:
        A tuple, (cost, path).
    """
    # start_position S(0,0)
    S = (0, 0)
    n = len(velocity)
    stoe_dists = [[] for i in range(n)]
    stos_dists = [[] for i in range(n)]
    etoe_dists = [[] for i in range(n)]
    Stos_dists = []
    Stoe_dists = []

    for i in range(n):
        for j in range(n):
            stoe_dists[i].append(dist_cal(task_start[i], task_end[j]))
            stos_dists[i].append(dist_cal(task_start[i], task_start[j]))
            etoe_dists[i].append(dist_cal(task_end[i], task_end[j]))

    for i in range(n):
        Stos_dists.append(dist_cal(S, task_start[i]))
        Stoe_dists.append(dist_cal(S, task_end[i]))

    C = {}

    for k in range(1, n):
        C[(1 << k * 2, k, 0)] = (Stos_dists[k], 0, 0)
        C[(1 << k * 2 | 1 << k * 2 + 1, k, 1)] = (Stoe_dists[k], 0, 0)

    # Iterate subsets of increasing length and store intermediate results
    # in classic dynamic programming manner
    for subset_size in range(2, n):
        for subset in itertools.combinations(range(1, n), subset_size):
            for binary in make_bin(subset):
                # Set bits for all nodes in this subset
                bits = 0
                for bit in subset:
                    bits |= 1 << bit * 2
                for bit in binary:
                    bits |= 1 << bit * 2 + 1

                # Find the lowest cost to get to this subset
                for k in subset:
                    prev = bits & ~(1 << k * 2)
                    if k in binary:
                        prev &= ~(1 << k * 2 + 1)

                    res_start = []
                    res_end = []
                    for m in subset:
                        print(subset)
                        print(m)
                        print(binary)
                        print(prev)
                        print(bin(prev))
                        if m == 0 or m == k:
                            continue
                        if m in binary:
                            res_start.append((C[(prev, m, 1)][0] + stos_dists[m][k], m, 0))
                            res_end.append((C[(prev, m, 1)][0] + stoe_dists[m][k], m, 1))
                        else:
                            res_start.append((C[(prev, m, 0)][0] + stoe_dists[k][m], m, 0))
                            res_end.append((C[(prev, m, 0)][0] + etoe_dists[m][k], m, 1))
                    C[(bits, k, 0)] = min(res_start)
                    C[(bits, k, 1)] = min(res_end)


    # We're interested in all bits but the least significant (the start state)
    bits = (2**(n * 2) - 1) - 1

    # Calculate optimal cost
    res = []
    for k in range(1, n):
        res.append((C[(bits, k, 0)][0] + Stos_dists[k][0], k, 0))
        res.append((C[(bits, k, 1)][0] + Stoe_dists[k][0], k, 1))
    opt, parent, s_or_e = min(res)

    # Backtrack to find full path
    path = []
    path_which = []
    for i in range(n - 1):
        path.append(parent)
        path_which.apppend(s_or_e)
        new_bits = bits & ~(1 << parent * 2)
        if s_or_e == 1:
            new_bits &= ~(1 << parent * 2 + 1)
        _, parent, s_or_e = C[(bits, parent, s_or_e)]
        bits = new_bits

    # Add implicit start state
    path.append(0)

    return opt, list(reversed(path)), list(reversed(path_which))

def generate_distances(n):
    task_start = []
    task_end = []
    velocity = []
    for i in range(n):
        task_start.append([random.randint(1, 999),random.randint(1,999)])
        task_end.append([random.randint(1, 999),random.randint(1,999)])
        velocity.append(random.randint(1,99))

    return task_start,task_end,velocity


def read_distances(filename):
    task_start = []
    task_end = []
    velocity = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comments
            if line[0] == '#':
                continue
            else:
                l = list(map(int, (line.strip()).split(',')))
                task_start.append(tuple(l[0:2]))
                task_end.append(tuple(l[2:4]))
                velocity.append(l[4])

    return task_start, task_end, velocity


def dist_cal(start, end):
    dist = math.sqrt((start[0] - end[0]) ** 2 + (start[1] - end[1]) ** 2)
    return dist

def make_bin(subset):
    binary = []
    subset_size = len(subset)
    for i in range(subset_size):
        for j in itertools.combinations(subset,i):
            binary.append(j)
    return binary

if __name__ == '__main__':
    arg = sys.argv[1]

    if arg.endswith('.csv'):
        task_start, task_end, velocity = read_distances(arg)
        n = len(velocity)
        print(task_start,task_end,velocity)

    else:
        task_start, task_end, velocity = generate_distances(int(arg))
        print(task_start)
        print(task_end)

    # Pretty-print the distance matrix
    print('')
    start = time.time()
    print(held_karp(task_start, task_end, velocity))
    process_time = time.time() - start
    print(process_time)