# problem 3 diffGem HW1

from sympy import symbols, diff, sin, cos, sqrt, simplify, Matrix

X_1, X_2, X_3 = symbols('X_1 X_2 X_3', positive=True)

A = Matrix([[1, 0, 0],
            [0, sqrt(3)/2, -1/2],
            [0, 1/2, sqrt(3)/2]])

X = [X_1,
     X_2,
     X_3]

x_sphere = [X_1 * sin(X_2) * cos(X_3),
            X_1 * sin(X_2) * sin(X_3),
            X_1 * cos(X_2)]


Q_mixed = Matrix([[0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0]])

Jacobian = Matrix([[0.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0]])

for i in range(3):
    for k in range(3):
        Jacobian[i, k] = diff(x_sphere[i], X[k])

# (1)
for i in range(3):
    for k in range(3):
        for j in range(3):
            Q_mixed[i, k] += A[i, j] * Jacobian[j, k]


# The latex output is too ugly. TODO(Tony): rewrite it

# (2)
r_covar_1 = [Q_mixed[0, 0], Q_mixed[1, 0], Q_mixed[2, 0]]
r_covar_2 = [Q_mixed[0, 1], Q_mixed[1, 1], Q_mixed[2, 1]]
r_covar_3 = [Q_mixed[0, 2], Q_mixed[1, 2], Q_mixed[2, 2]]

r_covar = [r_covar_1, r_covar_2, r_covar_3]

# (3)
g_covar = Matrix([[1, 0, 0],
                  [0, X_1**2, 0],
                  [0, 0, X_1**2 * (sin(X_2))**2]])

g_contra = g_covar**(-1)

# (4)
Q_contra = Matrix([[0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0]])

W = Matrix([[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]])

for i in range(3):
    for m in range(3):
        for j in range(3):
            Q_contra[i, m] += g_contra[i, j] * Jacobian[m, j]

for i in range(3):
    for j in range(3):
        for m in range(3):
            W[i, j] += A[j, m] * Q_contra[i, m]

# (5)
Christoffel_mixed = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                     Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                     Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

for i in range(3):
    for j in range(3):
        for m in range(3):
            for k in range(3):
                Christoffel_mixed[m][i, j] += 1/2 * g_contra[k, m] * (diff(g_covar[k, j], X[i]) + diff(g_covar[i, k], X[j]) - diff(g_covar[i, j], X[k]))

Christoffel_covar = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                     Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                     Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

for k in range(3):
    for i in range(3):
        for j in range(3):
            for m in range(3):
                if (i == 0) & (j == 0):
                    print(Christoffel_mixed[m][i, j] * g_covar[m, k])
                Christoffel_covar[i][j, k] += Christoffel_mixed[m][i, j] * g_covar[m, k]
print(Christoffel_covar)
# (6)
H = [0, 0, 0]
for i in range(3):
    H[i] = sqrt(g_covar[i, i])
    simplify(H[i])
