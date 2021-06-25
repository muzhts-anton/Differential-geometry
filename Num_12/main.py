from sympy import *
a = Matrix([0, 1])
T = Matrix([
    [1, 1],
    [0, 0]
])
Q = Matrix([
    [1, 1],
    [0, 1]
])

# контрвариантные компоненты вектора b в декартовом ортонорм. базисе e
# b^i = T^ij * aj
b_ = Matrix([0, 0])
for i in range(2):
    for j in range(2):
        b_[i] += T[i, j] * a[j]

P = Q**(-1)

# компоненты b в базисе e
# b^j = b_^i * Pi^j
b = Matrix([0, 0])
for j in range(2):
    for i in range(2):
        b[j] += b_[i] * P[j, i]

# метрическая матрица g[i,j] = ei * ej
g = Matrix([
    [0, 0],
    [0, 0]
])
for i in range(2):
    for j in range(2):
        for m in range(2):
            g[i, j] += Q[m, i] * Q[m, j]

# ковариантные компоненты вектора b в базисе e
# bi = gij * b^j
B = Matrix([0, 0])
for i in range(2):
    for j in range(2):
        B[i] += g[i, j] * b[j]

# контрвариантные компоненты вектора c в декартовом ортонорм. базисе e
# c^j = ai * T^ij
c_ = Matrix([0, 0])
for j in range(2):
    for i in range(2):
        c_[j] += a[i] * T[i, j]

# компоненты c в базисе e
# c^j = c_^i * Pi^j
c = Matrix([0, 0])
for j in range(2):
    for i in range(2):
        c[j] += c_[i] * P[j, i]

# ковариантные компоненты вектора c в базисе e
# ci = gij * c^j
C = Matrix([0, 0])
for i in range(2):
    for j in range(2):
        C[i] += g[i, j] * c[j]

# print(B, C)