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
    # start_position R0(0,0)
    R_0 = (0, 0)
    n = len(velocity)
    stoe_dists = [[] for i in range(n)]
    stos_dists = [[] for i in range(n)]
    etoe_dists = [[] for i in range(n)]

    for i in range(n):
        for j in range(n):
            stoe_dists[i].append(dist_cal(task_start[i], task_end[j]))
            stos_dists[i].append(dist_cal(task_start[i], task_start[j]))
            etoe_dists[i].append(dist_cal(task_end[i], task_end[j]))

    C = {}
    # Iterate subsets of increasing length and store intermediate results
    # in classic dynamic programming manner
    for subset_size in range(2, n):
        for subset in itertools.combinations(range(1, n), subset_size):
            for bin in itertools
            # Set bits for all nodes in this subset
            bits = 0
            for bit in subset:
                bits |= 1 << bit * 2
                bits |=

            # Find the lowest cost to get to this subset
            for k in subset:
                prev = bits & ~(1 << k)

                res = []
                for m in subset:
                    if m == 0 or m == k:
                        continue
                    res.append((C[(prev, m)][0] + dists[m][k], m))
                C[(bits, k)] = min(res)

    # We're interested in all bits but the least significant (the start state)
    bits = (2**n - 1) - 1

    # Calculate optimal cost
    res = []
    for k in range(1, n):
        res.append((C[(bits, k)][0] + dists[k][0], k))
    opt, parent = min(res)

    # Backtrack to find full path
    path = []
    for i in range(n - 1):
        path.append(parent)
        new_bits = bits & ~(1 << parent)
        _, parent = C[(bits, parent)]
        bits = new_bits

    # Add implicit start state
    path.append(0)

    return opt, list(reversed(path))

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
