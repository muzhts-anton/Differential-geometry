# T**B = T^ij * Bji
from sympy import *
T = Matrix([
    [1, 2],
    [0, 7]
])
B = Matrix([
    [3, -1],
    [2, 1]
])
rez = 0
for i in range(2):
    for j in range(2):
        rez += T[i, j] * B[j, i]
print(rez)
