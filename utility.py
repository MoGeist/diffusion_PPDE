#!/usr/bin/env python

from ast import literal_eval
from fenics import *
from itertools import combinations_with_replacement
from math import floor, ceil
import numpy as np
import matplotlib.pyplot as plt


def update_params(argv, params):
	'''	Update default dictionary with system arguments'''
	
	if len(argv) > 1:
		params.update(literal_eval(argv[1]))
	
	return params


def matrix_sqrt(mat):
	'''Caculate matrix square root of SPD matrix via eigendecomposition'''
	
	w, V = np.linalg.eigh(mat)
	
	if np.min(w) < 0:
		print('Smallest eigenvalue: ', np.min(w))
		print('Warning negative eigenvalues set to Zero')
		diag_mat = np.diag(np.sqrt(np.maximum(0, w))) 
	else:
		diag_mat = np.diag(np.sqrt(w))
	
	mat_sqrt = np.linalg.multi_dot([V, diag_mat, V.T])

	return mat_sqrt


def sqrt_gram_matrix(V):
	'''Calculate square root of Gram matrix'''
	
	u = TrialFunction(V)
	v = TestFunction(V)
	G1 = assemble(u*v*dx).array()	
	G2 = assemble(dot(grad(u), grad(v))*dx).array()
	
	G = G1 + G2

	return matrix_sqrt(G.astype(np.float32))

def get_polynomial(num_var, max_deg):
	'''Generate list of all possible monomials 

	num_var:	number of variables 
	max_deg:	maximal degree of monomials
	'''		
	
	var = ['c'] # add dummy variable
	for i in range(num_var):
		var.append('x[' + str(i) + ']')
	
	# compute all combinations of variables
	terms = []
	for x in combinations_with_replacement(var, max_deg):
		terms.append(list(x))
	
	for el in terms:
		while "c" in el:
			el.remove("c")
	
	monomial_list = ['*'.join(x) for x in terms]
	monomial_list[0] = '1' # fix constant term
	
	return monomial_list


def trigonometric_basis(n, p, mu=1, sigma=-1):
	'''Generate list of expressions and array of coefficients in polynomial basis (T1)

	n: 		number of samples
	p: 		number of basis fcts.
	mu: 	shift
	sigma: 	decay rate
	'''	

	num_terms = p
	
	coeff = np.random.uniform(0, 1., (n, num_terms))
	expr_list = []
	
	template_a = "+(1+{}*{}*(1+{}))"
	template_sin = "sin({}*3.14159*x[0])*sin({}*3.14159*x[1])"	

	for i in range(n):
		temp = str(mu)  
		for j in range(1,num_terms+1):
			sin = template_sin.format(floor((j+2)/2), ceil((j+2)/2))
			a = template_a.format(j**(sigma), str(coeff[i,j-1]), sin) 
			temp += a 
		expr_list.append(temp)

	return expr_list, coeff				


def squares_basis(n, s, mu=0.1):
	'''Generate list of expressions and array of coefficients in chessboard basis (T2)

	n: 	number of samples
	s: 	number of squares per row -> p = s^2
	mu:	shift

	Warning: Currently only works on unit square (in other cases adjust scaling)
	'''	
	
	num_terms = s**2
	
	coeff = np.random.uniform(0.0, 1., (n, num_terms))
	expr_list = []

	step = 1 / s

	# indicator fct. for x[0] in interval [low, up]
	template_x0 = "*ceil(fmax(x[0]-{low},0))*ceil(fmax({up}-x[0],0))"
	# indicator fct. for x[1] in interval [low, up]
	template_x1 = "*ceil(fmax(x[1]-{low},0))*ceil(fmax({up}-x[1],0))"

	for i in range(n):
		temp = str(mu)
		count = 0
		for j in range(s):
			for k in range(s):				
				indicator_x0 = template_x0.format(low=j*step, up=(j+1)*step)
				indicator_x1 = template_x1.format(low=k*step, up=(k+1)*step)
				temp += '+' + str(coeff[i,count]) + indicator_x0 + indicator_x1
				count += 1		
		expr_list.append(temp)	
	
	return expr_list, coeff				


def cookies_basis(n, s, mu=0.1, mode='var'):
	'''Generate list of expressions and array of coefficients in cookie basis (T3)

	n: 		number of samples
	s: 		number of squares per row -> p = s^2
	mu:		shift
	mode:	'var' for variable radius, else: fixed radius 

	Warning: Currently only works on unit square (in other cases adjust scaling)
	'''	
	
	num_terms = s**2
	coeff = np.random.uniform(0.0, 1., (n, num_terms))

	if mode == "var":
		coeff_radii = np.random.uniform(0.5, 0.9, (n, num_terms))
		coeff = np.concatenate((coeff, coeff_radii), axis=1) 	

	expr_list = []

	step = 1 / s

	template_dist = 'sqrt(pow(x[0]-{c_x0},2)+pow(x[1]-{c_x1},2))'
	template_cookie = '+{}*ceil(fmax({radius}-{dist},0))' 

	for i in range(n):
		temp = str(mu) 
		count = 0
		for j in range(s):
			for k in range(s):
				if mode == "var":
					r = coeff[i, count+num_terms] / (2*s)
				else:
					r = 0.8 / (2*s)	 
				d = template_dist.format(c_x0=1/(2*s)+j*step, c_x1=1/(2*s)+k*step)
				cookie = template_cookie.format(coeff[i, count], radius=r, dist=d)
				temp += cookie
				count += 1			
		expr_list.append(temp)	
	
	return expr_list, coeff	


def polynomial_basis(n, k, mu=0.1):
	''' Generate list of expressions and coefficients in polynomial basis (T4)

	n: 	number of samples
	k: 	maximal degree of polynomial -> p = (k+2) choose 2 
	mu:	shift
	'''	

	poly = get_polynomial(2, k)
	num_terms = len(poly)	
	
	coeff = np.random.uniform(-1., 1., (n, num_terms))
	expr_list = []

	template = "fmax({},{})"

	for i in range(n):	
		temp = '' 
		for j in range(num_terms):
			temp += str(coeff[i, j]) + '*' + poly[j] + '+'
		expr_list.append(template.format(mu, temp[:-1]))

	return expr_list, coeff


if __name__ == '__main__':
	# plot some examples
	
	expr, coeff = trigonometric_basis(1,10, sigma=1)
	#expr, coeff = squares_basis(1,4)
	#expr, coeff = cookies_basis(1,3,mode='var')
	#expr, coeff = polynomial_basis(1,9)	

	mesh = UnitSquareMesh(100, 100)
	V = FunctionSpace(mesh, 'P', 1)
	u = interpolate(Expression(expr[0], degree=1), V)
	fig = plot(u)
	plt.colorbar(fig)
	plt.show()
