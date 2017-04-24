#!/usr/bin/python
from gurobipy import *

import sys, os
from datetime import datetime
import numpy as np
from numpy import linalg as LA

#user-created libraries
sys.path.append(os.path.abspath(os.getcwd() + '/utils'))
from matrix_utils import *
from tree import *



#==================================================================#
#========================= PREREQUISITES ==========================#
#==================================================================#

def print_help():
	print('Usage:')
	print('./e_persistent.py -h (--help)\n\t Prints help.')
	print('./e_persistent.py -f (--file) [path] \n\t Input matrix is generated from file in [path].')
	print('\t File must be Tab separated values with the first row and first column of headers.')
	print('\t Rows represent mutation reads and Columns represent samples.')
	print('./e_persistent.py -p (--persistent) [P]\n\t [P] is the maximum number of persistent character allowed.')
	print('./e_persistent.py -t (--time) [P]\n\t [time] is the maximum running time.')
	print('./e_persistent.py -c (--clones) [P]\n\t [clones] is the percentage of clones allowed.')
	print('\t The actual amount of maximum clones used is calculated by : [clones] * #Mutations.')
	sys.exit(1)

interrupted = False

input_file = ''
k_pers = 0
mutation_names = None
solution_name = ''
max_tree_depth = 10

#time to stop:
stop_time = 7200
mut_mod = 1

if len(sys.argv) > 1:
	ohelp = False
	ofile = ''
	opers = -1
	otime = -1
	oi = 1
	oclones = -1
	for o in sys.argv[1:]:
		try:
			if o in ['-h', '--help']:
				ohelp = True
			if o in ['-f', '--file']:
				ofile = sys.argv[oi + 1]
			if o in ['-p', '--persistent']:
				opers = int(sys.argv [oi + 1])
			if o in ['-t', '--time']:
				otime = int(sys.argv [oi + 1])
			if o in ['-c', '--clones']:
				oclones = float(sys.argv [oi + 1])
		except:
			print_help()
		oi += 1
	if ohelp:
		print_help()
	elif ofile != '' and opers >= 0 and otime >= 0 and oclones >= 0:
		input_file = ofile
		solution_name = input_file.split('/')[-1][:-4]
		k_pers = opers
		stop_time = otime
		mut_mod = oclones
	else:
		print_help()
else:
	print_help()

# Start program
input_matrix, mutation_names = import_matrix_tab(input_file)
matrix_name = os.path.basename(input_file).split('.')[0]

#Fixed parameters
num_samples = len(input_matrix)
num_mutations = len(input_matrix[0])
num_clones = int(num_mutations * mut_mod)
max_error = 1

print('Num samples: %d' % num_samples)
print('Num mutations: %d' % num_mutations)
print('Num clones: %d' % num_clones)

#==================================================================#
#========================== GUROBI MODEL ==========================#
#==================================================================#
start_model_time = datetime.now()
model = Model('Persistent Phylogeny Model')
model.setParam('Threads', 3)
model.setParam('TimeLimit', stop_time)

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
		y[c, 2*m] = model.addVar(vtype=GRB.BINARY,obj=0, name='y[{0},{1}]'.format(c, m))
		y[c, 2*m+1] = model.addVar(vtype=GRB.BINARY,obj=0, name='y[{0},{1}]'.format(c, m))
		m += 1
	c += 1

print('Generating variables B.. (iter: 4m^2)')
B = {}
p = 0
while p < 2*num_mutations:
	q = p + 1
	while q < 2*num_mutations:
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
# def generate_constraint(k_pers):
print('Generating constraint (1).. (iter: n * m)')
# split-row with errors
s = 0
while s < num_samples:
	m = 0
	while m < num_mutations:
		model.addConstr(
			quicksum(xp[s,c,m] for c in range(num_clones)) - input_matrix[s][m] <= error[s,m], 
			'(sr_e+)[{0},{1},{2}]'.format(s,m,c))
		model.addConstr(
			quicksum(xp[s,c,m] for c in range(num_clones)) - input_matrix[s][m] >= -error[s,m], 
			'(sr_e-)[{0},{1},{2}]'.format(s,m,c))
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
print('Generating constraint (1\'),(2\').. (iter: m^2)')
c = 0
while c < num_clones:
	m = 0
	while m < num_mutations:
		model.addConstr(y[c, 2*m] == y[c, 2*m + 1] + x[c, m], '(1\')[{0}, {1}]'.format(c,m))
		m += 1
	c += 1

# print('Generating constraint for max_tree_depth (iter: c)')
# c = 0
# while c < num_clones:
# 	model.addConstr(quicksum(y[c,m] for m in range(num_mutations*2)) <= max_tree_depth, 'max tree depth')
# 	c += 1

# (3')
print('Generating constraint for B')	
c = 0
while c < num_clones:
	p = 0
	while p < num_mutations*2:
		q = p + 1
		while q < num_mutations*2:
			model.addConstr(
				y[c, p] + y[c, q] - B[p, q, 1, 1] <= 1, 'B[{0},{1},1,1]_{2}'.format(p,q,c))
			model.addConstr(
				- y[c, p] + y[c, q] - B[p, q, 0, 1] <= 0, 'B[{0},{1},0,1]_{2}'.format(p,q,c))
			model.addConstr(
				y[c, p] - y[c, q] - B[p, q, 1, 0] <= 0, 'B[{0},{1},1,0]_{2}'.format(p,q,c))
			q += 1
		p += 1	
	c += 1

# (4')
print('Generating constraint for no conflicts')
p = 0
while p < num_mutations*2:
	q = p + 1
	if p %2 == 0:
		q +=1
	while q < num_mutations*2:
		model.addConstr(
			B[p,q,0,1] + B[p,q,1,0] +B[p,q,1,1] <= 2,
                    'Conf[{0},{1}]'.format(p,q))
		q += 1
	p += 1

# < k persistent characters
print('Generating constraint for persistent characters')
model.addConstr(
	quicksum(B[2*p,2*p+1,1,1] for p in range(num_mutations)) <= k_pers, 'Number of persistent characters')


print('Time to complete model creation: {0}'.format(datetime.now() - start_model_time))
#OPTIMIZE MODEL AND PRINT SOLUTION
print('#----- GUROBI OPTIMIZATION ----#')
print('#~~~~~ Total Persistent characters: {0}#'.format(k_pers))
start_optimize = datetime.now()
# model.optimize(callback_time)
model.optimize()


#==================================================================#
#======================= POST OPTIMIZATION ========================#
#==================================================================#


def print_solution():
	global f_out
	global final_error_interr
	if model.status == GRB.Status.OPTIMAL or model.status == GRB.Status.TIME_LIMIT:
		final_error = 0

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


		#print('~~~~~~ Clone matrix')
		# f_out.write('~~~~~~ Clone matrix\n')
		# print_lmatrix(uclone_matrix, f_out)

		#print('~~~~~~ Usage matrix')
		#f_out.write('~~~~~~ Usage matrix\n')
		usage_matrix = []
		s = 0
		while s < num_samples:
			c = 0
			row = []
			while c < num_clones:
				row.append(u[s,c].X)
				c +=1
			usage_matrix.append(row)
			s +=1
		# print_lmatrix(usage_matrix, f_out)


		test_matrix = np.dot(usage_matrix, clone_matrix)
		diff = input_matrix - test_matrix
		for i in range(num_samples):
			for j in range(num_mutations):
				final_error += abs(diff[i][j])

		test_accuracy = LA.norm(diff) / LA.norm(input_matrix)

		f_out.write('{0},{1},{2},{3},{4},{5},{6},{7},{8}\n'.format(
			matrix_name,
			num_samples,
			num_mutations,
			mut_mod,
			len(clone_used),
			k_pers,
			int((datetime.now() - start_optimize).total_seconds()),
			final_error,
			test_accuracy))
	else:
		f_out.write('{0},{1},{2},{3},{4},{5},{6},{7},{8}\n'.format(
			matrix_name,
			num_samples,
			num_mutations,
			mut_mod,
			num_clones,
			k_pers,
			int((datetime.now() - start_optimize).total_seconds()),
			'unfeasible',
			'unfeasible'))

		


f_out = open('test_results/pers.txt', 'a+')
if k_pers == num_mutations:
	f_out.close()
	f_out = open('test_results/kpers.txt', 'a+')
print_solution()
f_out.close()
