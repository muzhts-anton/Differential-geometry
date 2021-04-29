from sympy import *
a = Matrix([0, 1, 1])
T = Matrix([
    [-3, 2, 2],
    [2, 4, -5],
    [-1, -4, -5]
])

def e(i, j, k):
    if i == j or i == k or j == k:
        return 0
    elif i == 0 and j == 1 and k == 2:
        return 1
    elif i == 2 and j == 0 and k == 1:
        return 1
    elif i == 1 and j == 2 and k == 0:
        return 1
    else:
        return -1
 
# Am^l = A^i * T^kl * eikm
A = Matrix([
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
])
for l in range(3):
    for m in range(3):
        for i in range(3):
            for k in range(3):
                A[m, l] += a[i] * T[k, l] * e(i, k, m)
print(A)
