# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/



from sympy import *
T = Matrix([
    [1, 0],
    [6, 2]
])
B = Matrix([
    [8, 1],
    [3, 7]
])
Q = Matrix([
    [1, 0],
    [1, 2]
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
Cst = C^ml * gms * glt 
c = Matrix([
    [0, 0],
    [0, 0]
])
for s in range(2):
    for t in range(2):
        for m in range(2):
            for l in range(2):
                c[s, t] += C[m, l] * g[m, s] * g[l, t]
