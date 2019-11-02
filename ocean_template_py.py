#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ocean_template.py:    D-Wave Ocean SDKを用いた最適化用途のサンプリングテンプレートコード

#%%
# 本コードの実行のためには dwave-ocean-sdk モジュールがインストールされている
# 必要があります。以下の1行のコメントアウトを解除して、dwave-ocean-sdkを
# インストールすれば dimod, minorminerなどのD-Wave Ocean SDKに含まれる
# モジュールとともに numpy などの依存関係のあるパッケージもインストールされる。
#!pip install dwave-ocean-sdk
!pip show dwave-ocean-sdk
!pip show numpy

#%%
import dimod
from dwave.embedding import MinimizeEnergy, embed_bqm
from dwave.system import DWaveSampler
import math
import minorminer
import numpy as np
from pyqubo import Array, Constraint, Placeholder, solve_qubo

def dist_cal(x_1,y_1,x_2,y_2):
    dist = abs(math.sqrt(x_1 ** 2 + y_1 ** 2) - math.sqrt(x_2 ** 2 + y_2 ** 2))
    return dist

#%%
N = 64      # 変数の数
rho = 0.5   # QUBO行列のスパースネス（非ゼロ変数の数の割合）

# 問題インスタンスの生成
q = np.triu(np.random.normal(0, 1, [N, N]).astype(np.float64))
Nzero = int(math.floor(N * (N - 1) / 2 * (1.0 - rho)))
while Nzero > 0:
    (y, x) = sorted(np.random.randint(0, N, 2))
    if abs(q[y][x]) > 0:
        q[y][x] = 0
        Nzero -= 1

h = list(np.diag(q))
S = {}
J = {}
Q = {}

for i in range(N):
    for j in range(i, N):
        if abs(q[i][j]) > 0:
            if i == j:
                h[(i)] = q[i][i]
            else:
                S[(i, j)] = 1
                J[(i, j)] = q[i][j]

            Q[(i, j)] = q[i][j]
num = 5

start_x = np.round(1000 * np.random.rand())
start_y = np.round(1000 * np.random.rand())
task_x_start = np.round(1000 * np.random.rand(num))
task_y_start = np.round(1000 * np.random.rand(num))
task_x_end = np.round(1000 * np.random.rand(num))
task_y_end = np.round(1000 * np.random.rand(num))
velocity = np.round(100 * np.random.rand(num))

dists = []

for i in range(num):
    x = []
    y = []
    for j in range(i,num):
        dists.append(dist_cal(task_x_end[i],task_y_end[i],task_x_start[j],task_y_start[j]))

H_dists = sum()

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