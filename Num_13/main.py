# T**B = T^ij * Bji
from sympy import *
T = Matrix([
    [1, -4],
    [-5, -3]
])
B = Matrix([
    [-4, 2],
    [-5, -1]
])

rez = 0
for i in range(2):
    for j in range(2):
        rez += T[i, j] * B[j, i]

print(rez)
