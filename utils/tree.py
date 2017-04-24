#!/usr/bin/python
import sys, os, random
import numpy as np

f_out = None

hash_clone = {}

class Node:
	def __init__(self, name, parent, mut):
		self.name = name
		self.parent = parent
		self.children = []
		self.binary = []
		# print(mut)
		if parent:
			parent.children.append(self)
			self.binary = list(self.parent.binary)
			self.binary[mut] = 1
			key = ''.join(str(e) for e in self.binary)
			global hash_clone
			hash_clone[key] = self
		else:
			for i in range(mut):
				self.binary.append(0)
		# print(self.binary)

	def print_node(self):
		if not self.parent == None:
			print('"%s" -- "%s" [penwidth=3];' % (self.parent.name, self.name))
			global f_out
			if f_out:
				f_out.write('"%s" -- "%s" [penwidth=3];\n' % (self.parent.name, self.name))


def print_tree(node):
	if len(node.children) == 0:
		node.print_node()
	else:
		node.print_node()
		for child in node.children:
			print_tree(child)

def print_dot_tree(node, usage_matrix, names, clones, usage_tree=True):
	print('graph phylogeny {')
	f_out.write('graph phylogeny {\n')
	print_tree(node)

	if usage_tree:

		for i in range(len(usage_matrix)):
			# Randomize color and create sample node
			r = lambda: random.randint(25,255)
			color = '#%2X%2X%2X' % (r(),r(),r())
			f_out.write('{{ rank=sink; "Sample {0}" [shape=box color="{1}" fontcolor="{1}"]; }}\n'.format(i+1, color))

			for j in range(len(usage_matrix[0])):
				if usage_matrix[i][j] > 0:
					clone = clones[j]
					key = ''.join(str(e) for e in clone)
					global hash_clone
					c = hash_clone[key]
					f_out.write('"Sample {0}" -- "{1}" [label="{2:.4f}" color="{3}" fontcolor="{3}"];\n'
								.format(i+1, c.name, 
								usage_matrix[i][j], color))

	f_out.write('}\n')
	print('}')

# True if col1 contains col2
def contains(col1, col2):
	for i in range(len(col1)):
		if not col1[i] >= col2[i]:
			return False
	return True

# def identify_clones(matrix, names, mut_nod):
# 	clones = []
# 	for i in range(len(matrix)):
# 		c = []
# 		for j in range(len(matrix[0])):
# 			if matrix[i][j] == 1:
# 				c.append(mut_nod[names[j]])
# 		clones.append(c)
# 	return clones


def build_tree(matrix, names, usage, output_file):
	global f_out
	if not output_file == 'random':
		f_out = open(output_file, 'w+')
	rows = len(matrix)
	cols = len(matrix[0])

	dimensions = np.sum(matrix, axis=0)
	# ordered indeces
	indeces = np.argsort(dimensions)
	dimensions = np.sort(dimensions)

	mutations_name = []
	for i in range(cols):
		mutations_name.append(names[indeces[i]])
	
	# REMEMBER:
	# get the i-th columk with matrix[:,i]

	root = Node('root', None, cols)

	driver_mut = Node(mutations_name[-1], root, indeces[-1])

	mut_nod = {}

	mut_nod[mutations_name[cols-1]] = driver_mut

	i = cols - 2
	while i >=0:
		if dimensions[i] == 0:
			break

		attached = False
		for j in range(i+1, cols):
			if contains(matrix[:, indeces[j]], matrix[:, indeces[i]]):
				node = Node(mutations_name[i], mut_nod[mutations_name[j]], indeces[i])
				mut_nod[mutations_name[i]] = node
				attached = True
				break

		if not attached:
			node = Node(mutations_name[i], root, indeces[i])
			mut_nod[mutations_name[i]] = node				
		i -=1

	if not output_file == 'random':
		# identified_clones = identify_clones(matrix, names, mut_nod)
		print_dot_tree(root, usage, names, matrix)
	return root, mut_nod
