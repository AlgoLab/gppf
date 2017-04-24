#!/usr/bin/python

from gurobipy import *

import sys, os
from datetime import datetime

#user-created libraries
sys.path.append(os.path.abspath(os.getcwd() + '/utils'))
from matrix_utils import *


#==================================================================#
#========================= PREREQUISITES ==========================#
#==================================================================#

#FIXED
#Big M number for penalize conflicts
bigM = 10000


#FIXED
#input frequency matrix in examples

def print_help():
	print('Usage:')
	print('./e_dollo.py -h (--help)\n\t Prints help.')
	print('./e_dollo.py -f (--file) [path] \n\t Input matrix is generated from file in [path].')
	print('\t File must be Tab separated values with the first row and first column of headers.')
	print('\t Rows represent mutation reads and Columns represent samples.')
	print('./e_dollo.py -k [K]\n\t [K] for the appropriate k-Dollo model.')
	sys.exit(1)

interrupted = False

def callback_stop(model, where):
	if where == GRB.Callback.MIPSOL:
		if GRB.Callback.MIPSOL_SOLCNT > 0:
			if model.cbGet(GRB.Callback.MIPSOL_OBJ) <= 0.01:
				global interrupted 
				interrupted = True
				model.terminate()
			

input_file = ''
k_dollo = 0
import_type = ''

if len(sys.argv) > 1:
	ohelp = False
	ofile = ''
	oK = -1
	oi = 1
	for o in sys.argv[1:]:
		try:
			if o in ['-h', '--help']:
				ohelp = True
			if o in ['-f', '--file']:
				ofile = sys.argv[oi + 1]
			if o == '-k':
				ok = int(sys.argv [oi + 1])
		except:
			print_help()
		oi += 1
	if ohelp:
		print_help()
	elif ofile != '' and ok >= 0:
		input_file = ofile
		solution_name = input_file.split('/')[-1][:-4]
		k_dollo = ok
	else:
		print_help()
else:
	print_help()

# Start program
input_matrix, mutation_names = import_matrix_tab(input_file)


#Fixed parameters
num_samples = len(input_matrix)
num_mutations = len(input_matrix[0])
num_clones = num_mutations
max_error = 1
# num_clones = 5
# columns_constraints = unique_extend_enum(num_mutations)

print('Num samples: %d' % num_samples)
print('Num mutations: %d' % num_mutations)
print('Num clones: %d' % num_clones)

#==================================================================#
#========================== GUROBI MODEL ==========================#
#==================================================================#
start_model_time = datetime.now()
model = Model('Persistent Phylogeny Model')

#---------------------------------------------------#
#------------------- VARIABLES ---------------------#
#---------------------------------------------------#

#------------------VARIABLE----------------
#u[s,c] = continuous
print('Generating variables u.. (iter: n * m)')
u = {}
s = 0
while s < num_samples:
	c = 0
	while c < num_clones:
		u[s,c] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, obj=0, name='u[{0},{1}]'.format(s,c))
		c +=1
	s +=1

#------------------VARIABLE----------------
# error[i,j]
error = {}
s = 0
while s < num_samples:
	m = 0
	while m < num_mutations:
		error[s,m] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, obj=1, name='error[{0},{1}]'.format(s,m))
		m += 1
	s += 1  

#------------------VARIABLE----------------
#x[c,m] = 1 iff clone c has mutation m
print('Generating variables x.. (iter: m^2)')
x = {}
c = 0
while c < num_clones:
	m = 0
	while m < num_mutations:
		x[c,m] = model.addVar(vtype=GRB.BINARY, obj=0, name='c[{0}, {1}]'.format(c,m))
		m += 1
	c += 1

print('Generating variables xp.. (iter: n * m^2)')
xp = {}
s = 0
while s < num_samples:
	c = 0
	while c < num_clones:
		m = 0
		while m < num_mutations:
			xp[s,c,m] = model.addVar(vtype=GRB.CONTINUOUS, obj=0, name='xp[{0},{1},{2}]'.format(s,c, m))
			m += 1
		c += 1
	s += 1

print('Generating variables y.. (iter: 2m^2)')
y = {}
c = 0
while c < num_clones:
	m = 0
	while m < num_mutations:
		k = 0
		while k <= k_dollo:
			# print('(%d, %d)' % (c, k_dollo*m+k))
			y[c, (k_dollo+1)*m + k] = model.addVar(vtype=GRB.BINARY,obj=0, name='y[{0},{1}]'.format(c, m))
			k += 1
		m += 1
	c += 1
print('Generating variables B.. (iter: 4m^2)')
B = {}
p = 0
while p < (k_dollo+1)*num_mutations:
	q = p + 1
	while q < (k_dollo+1)*num_mutations:
		B[p,q,1,1] = model.addVar(vtype=GRB.CONTINUOUS, obj=0, name='B[{0},{1},1,1]'.format(p,q))
		B[p,q,1,0] = model.addVar(vtype=GRB.CONTINUOUS, obj=0, name='B[{0},{1},1,0]'.format(p,q))
		B[p,q,0,1] = model.addVar(vtype=GRB.CONTINUOUS, obj=0, name='B[{0},{1},0,1]'.format(p,q))
		B[p,q,0,0] = 1
		q += 1
	p += 1

#-------------OBJECTIVE FUNCTION-----------
model.modelSense = GRB.MINIMIZE
model.update()

#---------------------------------------------------#
#------------------ CONSTRAINTS --------------------#
#---------------------------------------------------#
print('Generating constraint (1).. (iter: n * m)')
# split-row with errors
s = 0
while s < num_samples:
	m = 0
	while m < num_mutations:
		model.addConstr(
			quicksum(xp[s,c,m] for c in range(num_clones)) - input_matrix[s][m] <= error[s,m], '(sr_e+)[{0},{1},{2}]'.format(s,m,c))
		model.addConstr(
			quicksum(xp[s,c,m] for c in range(num_clones)) - input_matrix[s][m] >= -error[s,m], '(sr_e-)[{0},{1},{2}]'.format(s,m,c))
		model.addConstr(
			error[s,m] <= max_error * input_matrix[s][m], '(max_e)[{0},{1}]'.format(s,m))
		m += 1
	s += 1  


print('Generating constraint (2),(3),(4),(5).. (iter: n * m^2)')
s = 0
while s < num_samples:
	m = 0
	while m < num_mutations:
		c = 0
		while c < num_clones:
			model.addConstr(xp[s,c,m] >= 0, '(2)[{0},{1},{2}]'.format(s,m,c))
			model.addConstr(xp[s,c,m] <= x[c,m], '(3)[{0},{1},{2}]'.format(s,m,c))
			model.addConstr(xp[s,c,m] <= u[s,c], '(4)[{0},{1},{2}]'.format(s,m,c))
			model.addConstr(xp[s,c,m] >= u[s,c] + x[c,m] - 1, '(5)[{0},{1},{2}]'.format(s,m,c))
			c += 1
		m += 1
	s += 1  

# sum <= 1
s = 0
while s < num_samples:
	model.addConstr(
		quicksum(u[s,c] for c in range(num_clones)) <= 1, 'sum[s]'.format(s))
	s += 1



# (1') (2')
c = 0
while c < num_clones:
	m = 0
	while m < num_mutations:
		model.addConstr(quicksum(y[c, (k_dollo+1)*m + k] for k in range(1, k_dollo+1)) <= y[c, (k_dollo+1)*m], '(1\')[{0}, {1}]'.format(c,m))
		m += 1
	c += 1

c = 0
while c < num_clones:
	m = 0
	while m < num_mutations:
		model.addConstr(y[c, (k_dollo+1)*m] - quicksum(y[c, (k_dollo+1)*m + k] for k in range(1, k_dollo+1)) == x[c,m], '(1\')[{0}, {1}]'.format(c,m))
		m += 1
	c += 1

# (3')
c = 0
while c < num_clones:
	p = 0
	while p < (k_dollo+1)*num_mutations:
		if p % (k_dollo+1) == 0:
			q = p + k_dollo + 1
		else:
			q = p + 1
		while q < (k_dollo+1)*num_mutations:
			model.addConstr(
				y[c, p] + y[c, q] - B[p, q, 1, 1] <= 1, 'B[{0},{1},1,1]_{2}'.format(p,q,c))
			model.addConstr(
				- y[c, p] + y[c, q] - B[p, q, 0, 1] <= 0, 'B[{0},{1},0,1]_{2}'.format(p,q,c))
			model.addConstr(
				y[c, p] - y[c, q] - B[p, q, 1, 0] <= 0, 'B[{0},{1},1,0]_{2}'.format(p,q,c))
			q += 1
		p += 1	
	c += 1

p = 0
while p < (k_dollo+1)*num_mutations:
	q = p + 1
	while q < (k_dollo+1)*num_mutations:
		model.addConstr(
			B[p,q,0,1] + B[p,q,1,0] +B[p,q,1,1] <= 2,
                    'Conf[{0},{1}]'.format(p,q))
		q += 1
	p += 1	



print('Time to complete model creation: {0}'.format(datetime.now() - start_model_time))
#OPTIMIZE MODEL AND PRINT SOLUTION
print('#----- GUROBI OPTIMIZATION ----#')
print('#~~~~~ Total K Dollo characters: {0}#'.format(k_dollo))
start_optimize = datetime.now()
model.optimize(callback_stop)


#==================================================================#
#======================= POST OPTIMIZATION ========================#
#==================================================================#

def print_solution():
	print('#----- SOLUTIONS -------#')
	global f_out
	if model.status == GRB.Status.OPTIMAL or interrupted:
		if not interrupted:
			f_out.write('Time to optimize: {0}\n'.format(datetime.now() - start_optimize))
			print('Optimal objective value: {0}'.format(model.getAttr(GRB.Attr.ObjVal)))
		else:
			print('Interrupted')
		clone_used = []
		c = 0
		while c < num_clones:
			s = 0
			us = False
			while s < num_samples:
				if u[s,c].X > 0:
					us = True
				s +=1
			if us:
				clone_used.append(c)
			c +=1

		print('~~~~~~~~~~~~Total clone used: {0}'.format(len(clone_used)))

		#generate clones
		uclone_matrix = []
		c = 0
		while c < num_clones:
			if c in clone_used:
				clone = []
				m = 0
				while m < num_mutations:
					clone.append(int(x[c,m].X))
					m += 1
				uclone_matrix.append(clone)
			c += 1

		clone_matrix = []
		c = 0
		while c < num_clones:
			clone = []
			m = 0
			while m < num_mutations:
				clone.append(int(x[c,m].X))
				m += 1
			clone_matrix.append(clone)
			c += 1

		# print proportions
		print('~~~~~~ Split-Row\nf[sample, n. clone]\t\t[clone]')
		f_out.write('~~~~~~ Split-Row\nf[sample, n. clone]\t\t[clone]\n')
		s = 0
		while s < num_samples:
			c = 0
			while c < num_clones:
				if u[s,c].X > 0:
					print('u[{0},{1}] = {2:.4f}\t\t\t{3}'.format(s,c, u[s,c].X, clone_matrix[c]))
					f_out.write('u[{0},{1}] = {2:.4f}\t\t\t{3}\n'.format(s,c, u[s,c].X, clone_matrix[c]))
				c+=1
			s+=1

		print('~~~~~~ Clone matrix')
		f_out.write('~~~~~~ Clone matrix\n')
		print_lmatrix(uclone_matrix, f_out)

		print('~~~~~~ Usage matrix')
		f_out.write('~~~~~~ Usage matrix\n')
		s = 0
		while s < num_samples:
			c = 0
			row = ''
			while c < num_clones:
				row += '{0} '.format(u[s,c].X)
				c +=1
			print(row[:-1])
			f_out.write(row[:-1] + '\n')
			s +=1

		# print y
		print('~~~~~~ Extended matrix')
		f_out.write('~~~~~~ Extended matrix\n')
		c = 0
		while c < num_clones:
			m = 0
			row = ''
			while m < (k_dollo+1)*num_mutations:
				row += str(int(y[c,m].X)) +' '
				m+=1
			print(row)
			f_out.write(row + '\n')
			c += 1

		# print('~~~~~~ Usage matrix')
		# c = 0
		# while c < num_clones:
		# 	s = 0
		# 	row = ''
		# 	while s < num_samples:
		# 		for i in range(num_mutations):
		# 			row += str(u[s,c].X * clone_matrix[c][i]) + ' '
		# 		s += 1
		# 	print(row)
		# 	c += 1
		print('~~~~~~ Error matrix')
		f_out.write('~~~~~~ Error matrix\n')
		s = 0
		while s < num_samples:
			m = 0
			row = ''
			while m < num_mutations:
				row += '{0} '.format(error[s,m].X)
				m += 1
			print(row[:-1])
			f_out.write(row + '\n')
			s += 1  

	else:
		model.computeIIS()
		for c in model.getConstrs():
			if c.getAttr(GRB.Attr.IISConstr) > 0:
				print c.getAttr(GRB.Attr.ConstrName)
				f_out.write(str(c.getAttr(GRB.Attr.ConstrName)) + '\n')
		print('Solution not found')
		# k_pers += 1
		# print('Starting with {0}'.format(k_pers))

# generate_model()
# generate_constraint(0)

print('\n---------INPUT MATRIX---------')
print('Execution time: {0}'.format(datetime.now()-start_model_time))
f_out = open('results.txt', 'w+')
f_out.write('---------INPUT MATRIX---------\n')
print_lmatrix(input_matrix,f_out)
print_solution()


