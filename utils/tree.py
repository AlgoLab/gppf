#!/usr/bin/python
import sys
import os
import random
import numpy as np

f_app = None

hash_clone = {}


class Node:
    def __init__(self, name, parent, mut):
        self.name = name
        self.parent = parent
        self.children = []
        self.binary = []
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

    def print_node(self):
        if not self.parent == None:
            global f_app
            f_app.write('\t"%s" -- "%s";\n' % (self.parent.name, self.name))


def print_tree(node):
    if len(node.children) == 0:
        node.print_node()
    else:
        node.print_node()
        for child in node.children:
            print_tree(child)


def print_dot_tree(node, usage_matrix, names, clones):
    f_app.write('graph phylogeny {\n')
    print_tree(node)
    f_app.write('}\n')

# True if col1 contains col2


def contains(col1, col2):
    for i in range(len(col1)):
        if not col1[i] >= col2[i]:
            return False
    return True


def build_tree(matrix, names, usage, append_file):
    global f_app
    f_app = append_file
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

    root = Node('germline', None, cols)

    driver_mut = Node(mutations_name[-1], root, indeces[-1])

    mut_nod = {}

    mut_nod[mutations_name[cols-1]] = driver_mut

    i = cols - 2
    while i >= 0:
        if dimensions[i] == 0:
            break

        attached = False
        for j in range(i+1, cols):
            if contains(matrix[:, indeces[j]], matrix[:, indeces[i]]):
                node = Node(
                    mutations_name[i], mut_nod[mutations_name[j]], indeces[i])
                mut_nod[mutations_name[i]] = node
                attached = True
                break

        if not attached:
            node = Node(mutations_name[i], root, indeces[i])
            mut_nod[mutations_name[i]] = node
        i -= 1

    print_dot_tree(root, usage, names, matrix)
    return root, mut_nod
