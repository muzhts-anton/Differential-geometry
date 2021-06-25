from sympy import *
from copy import deepcopy
from numpy import array, linalg, dot
init_printing()

def matrixSimplify(matr):
	for i in range(len(matr)):
		for j in range(len(matr[i])):
			matr[i][j] = trigsimp(simplify(together(matr[i][j])))


def printRelations(x):
	for i in range(len(x)):
		print("x%d =" % (i + 1), x[i])
	print("")


def getJacobian(x, X):
	Q = [[0 for j in range(len(x))] for i in range(len(x))]
	for i in range(0, len(x)):
		for j in range(0, len(x)):
			Q[i][j] = diff(x[i], X[j])
	Q = array(Q)
	return Q


def printMatr(matr):
	lengths = [0 for i in range(len(matr[0]))]
	for i in range(len(matr)):
		for j in range(len(matr[i])):
			l = len(str(matr[i][j]))
			if (l > lengths[j]):
				lengths[j] = l
	maxl = max(lengths)
	for i in range(len(matr)):
		for j in range(len(matr[i])):
			l = len(str(matr[i][j]))
			print(" " * (maxl - l + 2) + str(matr[i][j]), sep="", end="")
		print("")
	print("")



if __name__ == '__main__':
	X1 = symbols('X1', positive=True)
	X2, X3 = symbols('X2 X3')

	x1 = X1*cos(X2)
	x2 = X1*sin(X2)
	x3 = X3 

	x = [x1, x2, x3]
	X = [X1, X2, X3]

	T = array([  # вводим тут матрицу, используйте sin(0) вместо 0
	[1, sin(0), sin(0)],
	[sin(0), -sqrt(2)/2, -sqrt(2)/2],
	[sin(0), sqrt(2)/2, -sqrt(2)/2]
	])


	print("Relations between xi and Xi:")
	printRelations(x)


	print("Jacobian matrix Q: ")
	Q = getJacobian(x, X)
	matrixSimplify(Q)
	print(Q, "\n")



    #Q^I_K = A^I_J * якобиеву матрицу 
	Ans = T @ Q # матричное умножение
	printMatr(Ans)

	