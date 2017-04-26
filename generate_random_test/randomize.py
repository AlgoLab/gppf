#!/usr/bin/python

from matrix_utils import *
import os, argparse, sys
import numpy as np
import random

def import_hudson(path):
	fi = open(path, 'r')
	matrix = []
	for line in fi:
		stripped = line.replace(' ', '')
		row = []
		for c in stripped[:-1]:
			row.append(int(c))
		if row != []:
			matrix.append(row)

	# print(matrix)
	
	base = []
	i = 0
	while i < len(matrix):
		j = 0
		row = []
		while j < len(matrix[0]):
			# print('(%d, %d)' % (i,j))
			if matrix[i][j] == 0:
				row.append(0)
			else:
				if matrix[i][j+1] == 0:
					row.append(1)
				else:
					row.append(0)
			j += 2
		base.append(row)
		i += 1

	npbase = np.array(base)
	for i in range(len(matrix[0])*0.4):
		npbase[np.random.randint(0,20)][np.random.randint(0,10)] = 1


	return np.array(matrix), npbase, len(matrix), len(matrix[0])


def main(fpath):
	filepath = fpath
	filename = os.path.basename(filepath).split('.')[0]
	outpath = os.path.abspath(os.path.dirname(filepath))

	clones_mod = 1


	n_sample = 5  #<---------------------------------------CHANGE Nr. SAMPLEs
	extended, b, rows, cols = import_hudson(filepath)

	clones_tot = int(rows * clones_mod)

	u = []
	for i in range(n_sample):
		p = random.randint(0, 80) #prob che sia zero
		row = []
		resto = random.randint(0, 100)
		sum = float(resto)
		for j in range(clones_tot): #colons - mutations
			#x = random.randint(0, 100)
			x = np.random.dirichlet((1,1))[0] * 100
			if x > p:
				sum += x
				row.append(x)
			else:
				row.append(0)
		row = map(lambda x: float("{0:.4f}".format(x/sum)), row)
		u.append(row)

	np.array(u)
	f = np.dot(u, b)

	with open('%s/%s.sim' % (outpath, filename), 'w+') as fout:
		for i in range(cols/2):
			fout.write('Mut_%d\t' % i)
			fout.write('Mut_%d-\t' % i)
		fout.write('\n')


	with open('%s/%s.sim' % (outpath, filename), 'a') as fout:
		np.savetxt(fout, f, delimiter='\t', fmt='%.4f')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('file')
	args = parser.parse_args()
	if not '.hudson' in args.file:
		print('Please use a .hudson file, generated with ms by Hudson [cite].\n')
		sys.exit(1)
	else:
		main(args.file)

