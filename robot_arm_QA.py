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
x = Array.create('x', shape = (row * column),vartype = 'BINARY')
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
H_dists = sum(x[0: column] * w0)
for t in range(1,row - 1):
    for j in range(column):
        H_dists += sum(np.multiply(x[column * t: column * (t + 1)] * w[j],x[column * t + j]))
H_dists += sum(x[column * (row - 1): column * row] * w0)
H_tasks = 0
for i in range(row):
    sum = 0
    for j in range(row):
        sum += sum(x[column * j + i: column * j + i + 2])
    H_tasks += p[i] * ((sum(sum(x[:n,i * 2:(i + 1) * 2])) - 1) ** 2)
H_tasks = Constraint(H_tasks, "tasks")
H_time = 0
for i in range(n):
    H_time += Pt * ((sum(x[i]) - 1) ** 2)
H_time = Constraint(H_time, "time")
H_cost = H_dists + Placeholder("tasks") * H_tasks + H_time * Placeholder("time")
model = H_cost.compile()
feed_dict = {'tasks': 2.0, 'time': 2.0}
Q, offset = model.to_qubo(feed_dict=feed_dict)

S = []

for i in range(n):
    for j in range(n * 2):
        S.append(('x[{0}][{1}]'.format(i, j), 'x[{0}][{1}]'.format(i, j)))
        for k in range(i + 1, n):
            S.append(('x[{0}][{1}]'.format(i, j), 'x[{0}][{1}]'.format(k, j)))
            if j % 2 == 0:
                S.append(('x[{0}][{1}]'.format(i, j), 'x[{0}][{1}]'.format(k, j + 1)))
            else:
                S.append(('x[{0}][{1}]'.format(i, j), 'x[{0}][{1}]'.format(k, j - 1)))
        for k in range(j + 1, n * 2):
            S.append(('x[{0}][{1}]'.format(i, j), 'x[{0}][{1}]'.format(i, k)))

print(S)

# この時点でIsing形式用のJ, h, BINARY形式用のQが生成済みである。
# ISING形式の場合
#bqm = dimod.BinaryQuadraticModel.from_ising(h, J)
# BINARY形式の場合
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
print(bqm)

# %%
url = "https://cloud.dwavesys.com/sapi"
token = "sigU-299eaf05a65136c85cdfe87be0618aadec5edc91"
solver_name = "DW_2000Q_5"

sampler = DWaveSampler(endpoint=url, token=token, solver=solver_name)

# minorminerでエンベディング
embedding = minorminer.find_embedding(S, sampler.edgelist)
#print('bqm : {0}'.format(bqm))
#print('embedding : {0}'.format(embedding))
bqm_embed = embed_bqm(bqm, embedding, sampler.adjacency)

# D-Waveによるサンプリング
result = sampler.sample(bqm_embed, num_reads=1000, postprocess="optimization", beta=3.0)

# minimize energyによる後処理
cbm = MinimizeEnergy(bqm, embedding)
unembedded, idx = cbm(result, list(embedding.values()))

# アンエンベッドされた解に関して、エネルギーの再計算や重複解などの整理をする
# 出力結果はdimod.SampleSetの形式 (Ocean SDKによる他のサンプリング結果と同じデータ形式)
sample = dimod.SampleSet.from_samples_bqm(unembedded, bqm, num_occurrences=result.record['num_occurrences']).aggregate()

#%%
# 出力のテンプレート
print(sample)
print(sample.record['sample'])
print(sample.record['energy'])
print(sample.record['num_occurrences'])

# %%
# 最低エネルギー状態の取り出し
print(sample.lowest())
print(sample.lowest().record['sample'])
print(sample.lowest().record['energy'])
print(sample.lowest().record['num_occurrences'])
