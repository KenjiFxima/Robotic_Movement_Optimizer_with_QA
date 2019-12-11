#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!pip install dwave-ocean-sdk
#!pip show dwave-ocean-sdk
#!pip show numpy

import dimod
from dwave.embedding import MinimizeEnergy, embed_bqm
from dwave.system import DWaveSampler
import math
import minorminer
import numpy as np
from pyqubo import Array, Constraint, Placeholder, solve_qubo, Binary
import sys
import re

def dist_cal(start, end):
    dist = math.sqrt((start[0] - end[0]) ** 2 + (start[1] - end[1]) ** 2)
    return dist

def read_distances(filename):
    tasks = []
    velocity = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comments
            if line[0] == '#':
                continue
            else:
                l = list(map(int, (line.strip()).split(',')))
                tasks.append(list(l[0:2]))
                tasks.append(list(l[2:4]))
                velocity.append(l[4])
    tasks = np.array(tasks)

    return tasks, velocity

#問題の読み込み
arg = sys.argv[1]
tasks, velocity = read_distances(arg)
n = len(velocity)
row = n
column = n * 2

#問題インスタンスの作成
x = []
for i in range(row * column):
    x.append(Binary('x[{}]'.format(i)))
x = Array(np.reshape(np.array(x),(1,row * column)))

#アームの開始位置
start = [0,0]

#重み(移動距離)の計算
w0 = []
w = [[] for i in range(column)]
wf = []

for i in range(column):
    w0.append(dist_cal(start, tasks[i]))
    for j in range(column):
        if i % 2 == 0:
            w[i].append(dist_cal(tasks[i + 1], tasks[j]))
        else:
            w[i].append(dist_cal(tasks[i - 1], tasks[j]))
    if i % 2 == 0:
        wf.append(dist_cal(tasks[i + 1], start))
    else:
        wf.append(dist_cal(tasks[i - 1], start))
w0 = np.array(w0)
w = np.array(w)
wf = np.array(wf)

p = []

for i in range(row):
    p.append(max(np.amax(w[:n,i * 2:(i + 1) * 2]),max(w0[i * 2:(i + 1) * 2]))+np.amax(w[i * 2 : (i + 1) * 2]))

Pt = max(p)

#コスト関数
H_dists = sum(x[0, :column] * w0)
for t in range(1, row):
    for i in range(column):
        H_dists += sum(np.multiply(w[i] * x[0, column * t: column * (t + 1)], x[0, column * (t - 1) + i]))
H_dists += sum(x[0, column * (row - 1): column * row] * wf)

#制約項-1
H_tasks = 0
for i in range(row):
    s = 0
    for j in range(row):
        s += x[0, column * j + (i * 2)] + x[0, column * j + (i * 2) + 1]
    H_tasks += p[i] * ((s - 1) ** 2)
H_tasks = Constraint(H_tasks, "tasks")

#制約項-2
H_time = 0
for i in range(row):
    H_time += Pt * ((sum(x[0, i * column: (i + 1) * column]) - 1) ** 2)
H_time = Constraint(H_time, "time")

H_cost = H_dists + Placeholder("tasks") * H_tasks + Placeholder("time") * H_time
model = H_cost.compile()

#制約項の重み
feed_dict = {'tasks': 1.0, 'time': 1.0}

#QUBOの作成
qubo,offset = model.to_qubo(feed_dict=feed_dict)

Q = {(int(re.search(r"x\[([0-9]+)\]", i)[1]),
       int(re.search(r"x\[([0-9]+)\]", j)[1])): v for (i, j), v in qubo.items()}

#埋め込み用のグラフ
S = list(Q.keys())

# この時点でIsing形式用のJ, h, BINARY形式用のQが生成済みである。
# ISING形式の場合
#bqm = dimod.BinaryQuadraticModel.from_ising(h, J)
# BINARY形式の場合
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

url = "https://cloud.dwavesys.com/sapi"
token = "sigU-299eaf05a65136c85cdfe87be0618aadec5edc91"
solver_name = "DW_2000Q_5"

sampler = DWaveSampler(endpoint=url, token=token, solver=solver_name)

# minorminerでエンベディング
embedding = minorminer.find_embedding(S, sampler.edgelist)
bqm_embed = embed_bqm(bqm, embedding, sampler.adjacency)

# D-Waveによるサンプリング
result = sampler.sample(bqm_embed, num_reads=1000, postprocess="optimization", beta=3.0)

# minimize energyによる後処理
cbm = MinimizeEnergy(bqm, embedding)
unembedded, idx = cbm(result, list(embedding.values()))

# アンエンベッドされた解に関して、エネルギーの再計算や重複解などの整理をする
# 出力結果はdimod.SampleSetの形式 (Ocean SDKによる他のサンプリング結果と同じデータ形式)
sample = dimod.SampleSet.from_samples_bqm(unembedded, bqm, num_occurrences=result.record['num_occurrences']).aggregate()

# 出力のテンプレート
print(sample)
print(sample.record['sample'])
print(sample.record['energy'])
print(sample.record['num_occurrences'])

# 最低エネルギー状態の取り出し
print(sample.lowest())
print(sample.lowest().record['sample'])
print(sample.lowest().record['energy'])
print(sample.lowest().record['num_occurrences'])