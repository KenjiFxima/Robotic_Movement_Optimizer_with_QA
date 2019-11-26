import math
import minorminer
import numpy as np
from pyqubo import Array, Constraint, Placeholder, solve_qubo
import sys

def dist_cal(start, end):
    dist = math.sqrt((start[0] - end[0]) ** 2 + (start[1] - end[1]) ** 2)
    return dist

def read_distances(filename):
    task = []
    s = []
    velocity = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comments
            if line[0] == '#':
                continue
            else:
                l = list(map(int, (line.strip()).split(',')))
                s.append(list(l[0:2]))
                s.append(list(l[2:4]))
                velocity.append(l[4])
    set = np.array(s)
    n = int(len(set) / 2)
    for i in range(n):
        task.append(set)
    task = np.array(task)

    return task, velocity

#問題の読み込み
arg = sys.argv[1]
task, velocity = read_distances(arg)
n = len(velocity)

# 問題インスタンスの生成
x = Array.create('x',shape = (n,n * 2),vartype = 'BINARY')
start = [0,0]

w0 = []
w = [[] for i in range(n * 2)]

for i in range(n * 2):
    w0.append(dist_cal(start, task[0][i]))
for i in range(n * 2):
    for j in range(n * 2):
        w[j].append(dist_cal(task[0][i],task[0][j]))

w0 = np.array(w0)
w = np.array(w)
p = []
for i in range(n):
    p.append(max(np.amax(w[:n,i*2:(i+1)*2]),np.amax(w0[i*2:(i+1)*2]))+np.amax(w[i]))
Pt = max(p)

H_dists = sum(x[0] * w0)
for t in range(1,n-1):
    for j in range(n * 2):
        H_dists += sum(np.multiply(x[t + 1] * w[j],x[t][j]))
H_dists = sum(x[n-1] * w0)
H_tasks = 0
for i in range(n):
    H_tasks += p[i] * (sum((sum(x[:n,i*2:(i+1)*2])) - 1) ** 2)
H_time = 0
for i in range(n):
    H_time += Pt * ((sum(x[i]) - 1) ** 2)

H_cost = H_dists + H_tasks + H_time
model = H_cost.compile()
Q, offset = model.to_qubo()
print(x)