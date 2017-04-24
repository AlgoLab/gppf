#!/usr/bin/python
from gurobipy import *

import sys, os
from datetime import datetime
import numpy as np

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
	sys.exit(1)

interrupted = False
final_error_interr = 0

def callback_stop(model, where):
	if where == GRB.Callback.MIPSOL:
		if GRB.Callback.MIPSOL_SOLCNT > 0:
			global num_mutations
			global num_clones
			treshold = num_mutations * num_clones * 0.01
			if model.cbGet(GRB.Callback.MIPSOL_OBJ) <= 0.6: #<----------Change according to data (0.03 * M)
				global interrupted 
				interrupted = True
				global final_error_interr
				final_error_interr = GRB.Callback.MIPSOL_OBJ
				model.terminate()
			

input_file = ''
k_pers = 0
mutation_names = None
solution_name = ''
max_tree_depth = 10

if len(sys.argv) > 1:
	ohelp = False
	ofile = ''
	opers = -1
	oi = 1
	for o in sys.argv[1:]:
		try:
			if o in ['-h', '--help']:
				ohelp = True
			if o in ['-f', '--file']:
				ofile = sys.argv[oi + 1]
			if o in ['-p', '--persistent']:
				opers = int(sys.argv [oi + 1])
		except:
			print_help()
		oi += 1
	if ohelp:
		print_help()
	elif ofile != '' and opers >= 0:
		input_file = ofile
		solution_name = input_file.split('/')[-1][:-4]
		k_pers = opers
	else:
		print_help()
else:
	print_help()

# Start program
input_matrix, mutation_names = import_matrix_tab(input_file)

#Fixed parameters
num_samples = len(input_matrix)
num_mutations = len(input_matrix[0])
#mutation modifier <--------------------------------------Change this
mut_mod = 1
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
print('Generating constraint (1\'),(2\').. (iter: m^2)')
c = 0
while c < num_clones:
	m = 0
	while m < num_mutations:
		model.addConstr(y[c, 2*m] == y[c, 2*m + 1] + x[c, m], '(1\')[{0}, {1}]'.format(c,m))
		m += 1
	c += 1

print('Generating constraint for max_tree_depth (iter: c)')
c = 0
while c < num_clones:
	model.addConstr(quicksum(y[c,m] for m in range(num_mutations*2)) <= max_tree_depth, 'max tree depth')
	c += 1

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
model.optimize(callback_stop)


#==================================================================#
#======================= POST OPTIMIZATION ========================#
#==================================================================#


def print_solution():
	#print('#----- SOLUTIONS -------#')
	global f_out
	global final_error_interr
	f_out.write(str(datetime.now() - start_optimize) + '\n')
	if model.status == GRB.Status.OPTIMAL or interrupted:
		if not interrupted:
			f_out.write(str(model.getAttr(GRB.Attr.ObjVal)) + '\n')
		else:
			f_out.write(str(final_error_interr) + '\n')
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

		#print('~~~~~~~~~~~~Total clone used: {0}'.format(len(clone_used)))

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
		# print('~~~~~~ Split-Row\nf[sample, n. clone]\t\t[clone]')
		# f_out.write('~~~~~~ Split-Row\nf[sample, n. clone]\t\t[clone]\n')
		# s = 0
		# while s < num_samples:
		# 	c = 0
		# 	while c < num_clones:
		# 		if u[s,c].X > 0:
		# 			print('u[{0},{1}] = {2:.4f}\t\t\t{3}'.format(s,c, u[s,c].X, clone_matrix[c]))
		# 			f_out.write('u[{0},{1}] = {2:.4f}\t\t\t{3}\n'.format(s,c, u[s,c].X, clone_matrix[c]))
		# 		c+=1
		# 	s+=1

		#print('~~~~~~ Clone matrix')
		# f_out.write('~~~~~~ Clone matrix\n')
		print_lmatrix(uclone_matrix, f_out)

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
		print_lmatrix(usage_matrix, f_out)

		#print('~~~~~~ Extended matrix')
		#f_out.write('~~~~~~ Extended matrix\n')
		# extended_matrix = []
		# c = 0
		# while c < num_clones:
		# 	m = 0
		# 	row = []
		# 	while m < num_mutations*2:
		# 		row.append(int(y[c,m].X))
		# 		m+=1
		# 	extended_matrix.append(row)
		# 	c += 1
		# print_lmatrix(extended_matrix, f_out)

		# print('~~~~~~ Error matrix')
		# f_out.write('~~~~~~ Error matrix\n')
		# s = 0
		# while s < num_samples:
		# 	m = 0
		# 	row = ''
		# 	while m < num_mutations:
		# 		row += '{0} '.format(error[s,m].X)
		# 		m += 1
		# 	print(row[:-1])
		# 	f_out.write(row + '\n')
		# 	s += 1  

		#print('~~~~~~ Tree (in "%s_usage.tree" it is also included the usage matrix)' % solution_name)
		#build_tree(np.array(extended_matrix), mutation_names, np.array(usage_matrix), solution_name + '_usage.tree')

	else:
		model.computeIIS()
		for c in model.getConstrs():
			if c.getAttr(GRB.Attr.IISConstr) > 0:
				print c.getAttr(GRB.Attr.ConstrName)
				f_out.write(str(c.getAttr(GRB.Attr.ConstrName)) + '\n')
		#print('Solution not found')


#print('\n---------INPUT MATRIX---------')
#print('Execution time: {0}'.format(datetime.now()-start_model_time))
f_out = open(solution_name + '_results.txt', 'w+')
#f_out.write('---------INPUT MATRIX---------\n')
#print_lmatrix(input_matrix,f_out)
print_solution()
f_out.close()
