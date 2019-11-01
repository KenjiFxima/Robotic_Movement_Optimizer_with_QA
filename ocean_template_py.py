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

# この時点でIsing形式用のJ, h, BINARY形式用のQが生成済みである。
#%%
# ISING形式の場合
#bqm = dimod.BinaryQuadraticModel.from_ising(h, J)
# BINARY形式の場合
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

# %%
url = "https://cloud.dwavesys.com/sapi"
token = "sigU-25f9897d301781a321aeca3e7d18d59c599234c3"
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