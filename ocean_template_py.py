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
x = Array.create('x',shape = (n * 2,5),vartype = 'BINARY')
start = [0,0]
print(task,velocity)

w0 = []
w = [[] for i in range(n * 2)]

for i in range(n * 2):
    w0.append(dist_cal(start, task[0][i]))
for i in range(n * 2):
    for j in range(n * 2):
        w[j].append(dist_cal(task[0][i]),task[0][j])

w0 = np.array(w0)
w = np.array(w)
p = []
for i in range(n):
    p.append(max(np.amax(w[:n,i*2:(i+1)*2]),max(w0[i*2:(i+1)*2]))+max(np.amax(w[i])))

h = list(np.diag(q))
S = {}
J = {}
Q = {}

H_dists = sum(x[0] * w0)
for i in range(1,n):
    for j in range(n * 2):
        H_dists += sum(np.multiply(x[i + 1] * w[j],x[i][j]))
H_dists = sum(x[n] * w0)

H_tasks = sum((sum(x) ** 2))

# この時点でIsing形式用のJ, h, BINARY形式用のQが生成済みである。
#%%
# ISING形式の場合
#bqm = dimod.BinaryQuadraticModel.from_ising(h, J)
# BINARY形式の場合
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

# %%
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
