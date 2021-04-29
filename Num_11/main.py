from sympy import *
T = Matrix([
    [-4, -2],
    [3, -1]
])
B = Matrix([
    [2, 3],
    [-1, -2]
])
Q = Matrix([
    [-2, 0],
    [-5, -3]
])

# метрическая матрица g[i,j] = ei * ej
g = Matrix([
    [0, 0],
    [0, 0]
])
for i in range(2):
    for j in range(2):
        for m in range(2):
            g[i, j] += Q[m, i] * Q[m, j]

P = Q**(-1)

# вычилим контравариантные компоненты тензора С в декартовом ортонормированном базисе
# C_^ik = T^ij * Bj^k
C_ = Matrix([
    [0, 0],
    [0, 0]
])
for i in range(2):
    for k in range(2):
        for j in range(2):
            C_[i, k] += T[i, j] * B[k, j]

# вычислим контраваринатные компонены С в в базисе е
# C^ml = C^ik * Pi^m * Pk^l
C = Matrix([
    [0, 0],
    [0, 0]
])
for m in range(2):
    for l in range(2):
        for i in range(2):
            for k in range(2):
                C[m , l] += C_[i, k] * P[m, i] * P[l, k]

# вычислим ковариантные компоненты С в базисе е
#Cst = C^ml * gms * glt 
c = Matrix([
    [0, 0],
    [0, 0]
])
for s in range(2):
    for t in range(2):
        for m in range(2):
            for l in range(2):
                c[s, t] += C[m, l] * g[m, s] * g[l, t]
print(c)
