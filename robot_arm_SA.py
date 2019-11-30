#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ocean_template.py:    D-Wave Ocean SDKを用いた最適化用途のサンプリングテンプレートコード

#%%
# 本コードの実行のためには dwave-ocean-sdk モジュールがインストールされている
# 必要があります。以下の1行のコメントアウトを解除して、dwave-ocean-sdkを
# インストールすれば dimod, minorminerなどのD-Wave Ocean SDKに含まれる
# モジュールとともに numpy などの依存関係のあるパッケージもインストールされる。
#!pip install dwave-ocean-sdk
#!pip show dwave-ocean-sdk
#!pip show numpy

import dimod
from dwave.embedding import MinimizeEnergy, embed_bqm
from dwave.system import DWaveSampler
import math
import minorminer
import numpy as np
from pyqubo import Array, Constraint, Placeholder, solve_qubo
import sys
from pyqubo import Binary

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
row = n
column = n * 2

# 問題インスタンスの生成
x = Array.create('x', shape = (1, row * column),vartype = 'BINARY')
start = [0,0]

w0 = []
w = [[] for i in range(column)]

for i in range(column):
    w0.append(dist_cal(start, task[0][i]))
for i in range(column):
    for j in range(column):
        w[j].append(dist_cal(task[0][i],task[0][j]))

w0 = np.array(w0)
w = np.array(w)
p = []
for i in range(row):
    p.append(max(np.amax(w[:n,i * 2:(i + 1) * 2]),max(w0[i * 2:(i + 1) * 2]))+np.amax(w[i]))
Pt = max(p)
H_dists = sum(x[0, :column] * w0)
for t in range(1,row - 1):
    for j in range(column):
        H_dists += sum(x[0, column * t: column * (t + 1)] * w[j]) * x[0, column * t + j]
H_dists += sum(x[0, column * (row - 1): column * row] * w0)
H_tasks = 0
for i in range(row):
    s = 0
    for j in range(row):
        s += sum(x[0, column * j + i * 2: column * j + (i + 1) * 2])
    H_tasks += p[i] * ((s - 1) ** 2)
H_tasks = Constraint(H_tasks, "tasks")
H_time = 0
for i in range(row):
    H_time += Pt * ((sum(x[0, i * column: (i + 1) * column]) - 1) ** 2)
H_time = Constraint(H_time, "time")
H_cost = H_dists + Placeholder("tasks") * H_tasks + Placeholder("time") * H_time
model = H_cost.compile()
feed_dict = {'tasks': 2.0, 'time': 2.0}
Q, offset = model.to_qubo(feed_dict=feed_dict)