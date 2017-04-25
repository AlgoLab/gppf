from gurobipy import *
from numpy import linalg as LA
import numpy as np
from datetime import datetime

def check_file_existance(path):
	if not os.path.exists(path):
		fo = open(path, 'w+')
		fo.write('matrix_name,num_samples,num_mutations,mut_mod,' +
					'clone_used,k,time,total_error,accuracy\n')
		fo.close()
	return open(path, 'a')

def print_exp_solution(num_clones, num_samples, num_mutations, start_optimize, 
							args, model, u, x, input_matrix, matrix_name):
	if args.model == 'dollo' and args.k == 0:
		args.model == 'perfect'
	f_out = check_file_existance('res_exp_{}.csv'.format(args.model))
	
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
			args.clones,
			len(clone_used),
			args.k,
			int((datetime.now() - start_optimize).total_seconds()),
			final_error,
			test_accuracy))
	else:
		f_out.write('{0},{1},{2},{3},{4},{5},{6},{7},{8}\n'.format(
			matrix_name,
			num_samples,
			num_mutations,
			args.clones,
			num_clones,
			args.k,
			int((datetime.now() - start_optimize).total_seconds()),
			'unfeasible',
			'unfeasible'))

def print_solution():
	print('TODO')