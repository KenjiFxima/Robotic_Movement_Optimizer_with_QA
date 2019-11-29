# 問題設定
NUM_VER = 6
vertices = list(range(NUM_VER))
edges = [(0,1), (0,4), (0,5), (1,2), (1,3), (3,4), (4,5)]
from pyqubo import Array, Constraint, Placeholder, solve_qubo

# BINARY変数
x = Array.create('x', shape=NUM_VER, vartype='BINARY')

# プレースホルダー
param_cover = Placeholder("cover")

# QUBO形式で定式化
H_cover = Constraint(sum((1-x[u])*(1-x[v]) for (u,v) in edges), "cover")
H_vertices = sum(x)
H = H_vertices + param_cover * H_cover

# モデルをコンパイル
model = H.compile()

# プレースホルダーと合わせてQUBOを作成
feed_dict = {"cover": param_cover}
qubo, offset = model.to_qubo(feed_dict=feed_dict)
print(qubo)