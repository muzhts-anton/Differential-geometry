# problem 4 diffGem HW1

from sympy import symbols, diff, sin, cos, Matrix

X_1, X_2, X_3 = symbols('X_1 X_2 X_3')

X = [X_1, X_2, X_3]

x = [X_1 * cos(X_2), X_1 * sin(X_2), X_3]

g_covar = Matrix([[1, 0, 0],
                  [0, X_1**2, 0],
                  [0, 0, 1]])

g_contra = Matrix([[1, 0, 0],
                  [0, (1 / X_1)**2, 0],
                  [0, 0, 1]])

T_contra = Matrix([[0, X_1, 0],
                   [X_3, 0, 0],
                   [2 * X_2, 0, 0]])

T_covar = Matrix([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0]])

# (1)
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                T_covar[i, j] += T_contra[k, l] * g_covar[k, i] * g_covar[l, j]


T_mixed_up = Matrix([[0, 0, 0],
                     [0, 0, 0],
                     [0, 0, 0]])

T_mixed_down = Matrix([[0, 0, 0],
                       [0, 0, 0],
                       [0, 0, 0]])

for i in range(3):
    for j in range(3):
        for k in range(3):
            T_mixed_up[i, j] += T_contra[i, k] * g_covar[k, j]

for i in range(3):
    for j in range(3):
        for k in range(3):
            T_mixed_down[i, j] += g_covar[k, i] * T_contra[k, j]

# (2)
Chr = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
       Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
       Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

Chr[0][1, 1] = - X_1
Chr[1][0, 1] = 1 / X_1
Chr[1][1, 0] = 1 / X_1


nabla_covar_T_covar = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                       Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                       Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

nabla_covar_T_mixed_up = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                          Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                          Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

nabla_covar_T_mixed_down = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                            Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                            Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

nabla_covar_T_contra = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                        Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                        Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]


nabla_contra_T_covar = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                        Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                        Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

nabla_contra_T_mixed_up = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                           Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                           Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

nabla_contra_T_mixed_down = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                             Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                             Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

nabla_contra_T_contra = [Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                         Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
                         Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]


for k in range(3):
    for i in range(3):
        for j in range(3):
            nabla_covar_T_contra[k][i, j] = diff(T_contra[i, j], X[k])
            for m in range(3):
                nabla_covar_T_contra[k][i, j] += T_contra[m, j] * Chr[i][m, k] + T_contra[i, m] * Chr[j][m, k]

for k in range(3):
    for i in range(3):
        for j in range(3):
            nabla_covar_T_covar[k][i, j] = diff(T_covar[i, j], X[k])
            for m in range(3):
                nabla_covar_T_covar[k][i, j] -= T_covar[m, j] * Chr[m][i, k] + T_covar[i, m] * Chr[m][j, k]

for k in range(3):
    for i in range(3):
        for j in range(3):
            nabla_covar_T_mixed_up[k][i, j] = diff(T_mixed_up[i, j], X[k])
            for m in range(3):
                nabla_covar_T_mixed_up[k][i, j] += T_mixed_up[m, j] * Chr[i][m, k] - T_mixed_up[i, m] * Chr[m][j, k]

for k in range(3):
    for i in range(3):
        for j in range(3):
            nabla_covar_T_mixed_down[k][i, j] = diff(T_mixed_down[i, j], X[k])
            for m in range(3):
                nabla_covar_T_mixed_down[k][i, j] += - T_mixed_down[m, j] * Chr[m][i, k] + T_mixed_down[i, m] * Chr[j][m, k]


for m in range(3):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                nabla_contra_T_contra[m][i, j] += g_contra[m, k] * nabla_covar_T_contra[k][i, j]

for m in range(3):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                nabla_contra_T_covar[m][i, j] += g_contra[m, k] * nabla_covar_T_covar[k][i, j]

for m in range(3):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                nabla_contra_T_mixed_up[m][i, j] += g_contra[m, k] * nabla_covar_T_mixed_up[k][i, j]

for m in range(3):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                nabla_contra_T_mixed_down[m][i, j] += g_contra[m, k] * nabla_covar_T_mixed_down[k][i, j]
