<<<<<<< HEAD
=======
import itertools
>>>>>>> ec550e7829356d401baee9f63536ed13a23081f9
import random
import sys
import time
import math
<<<<<<< HEAD
=======
import numpy as np
>>>>>>> ec550e7829356d401baee9f63536ed13a23081f9


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

    # C[total_dist,[connect to start or end],path]
    C = [[] for i in range(n)]

    for i in range(n):
        C[0].append([dist_cal(R_0, task_start[i]), [0], i])
        C[0].append([dist_cal(R_0, task_end[i]), [1], i])

    path_root = []
    path_which = []
    ans = []

    l = [i for i in range(n)]

    for i in range(n-1):
        for j in range(len(C[i])):
            m = list(iter(l))
            for s in range(2,len(C[i][j])):
                m.remove(C[i][j][s])

            for t in range(len(m)):
                if C[i][j][1][len(C[i][j][1])-1] == 0:
                    q = C[i][j][0] + stoe_dists[m[t]][C[i][j][len(C[i][j])-1]]
                    d = [q]
                    d.append(list(iter(C[i][j][1])))
                    d.extend(C[i][j][2:len(C[i][j])])
                    d.append(m[t])
                    d[1].append(0)
                    C[i + 1].append(d)

                    q = C[i][j][0] + etoe_dists[C[i][j][len(C[i][j])-1]][m[t]]
                    d  = [q]
                    d.append(list(iter(C[i][j][1])))
                    d.extend(C[i][j][2:len(C[i][j])])
                    d.append(m[t])
                    d[1].append(1)
                    C[i + 1].append(d)

                if C[i][j][1][len(C[i][j][1])-1] == 1:
                    q = C[i][j][0] + stos_dists[C[i][j][len(C[i][j])-1]][m[t]]
                    d = [q]
                    d.append(list(iter(C[i][j][1])))
                    d.extend(C[i][j][2:len(C[i][j])])
                    d.append(m[t])
                    d[1].append(0)
                    C[i + 1].append(d)

                    q = C[i][j][0] + stoe_dists[C[i][j][len(C[i][j])-1]][m[t]]
                    d = [q]
                    d.append(list(iter(C[i][j][1])))
                    d.extend(C[i][j][2:len(C[i][j])])
                    d.append(m[t])
                    d[1].append(1)
                    C[i + 1].append(d)

    for i in range(len(C[n-1])):
        if C[n-1][i][1][len(C[n-1][i][1]) - 1] == 0:
            q = C[n-1][i][0] + C[0][C[n-1][i][len(C[n-1][i]) - 1] * 2 + 1][0]
        else:
            q = C[n-1][i][0] + C[0][C[n-1][i][len(C[n-1][i]) - 1] * 2][0]
        ans.append(q)
        path_root.append(C[n-1][i][2:len(C[n-1][i])])
        path_which.append(C[n-1][i][1])

    opt = min(ans)
    idx = ans.index(opt)
    path = [path_root[idx], path_which[idx]]

    return opt, path


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
