from gurobipy import *
from numpy import linalg as LA
import numpy as np
from datetime import datetime
from tree import *
from matrix_utils import *


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
        args.model = 'perfect'
    if args.model == 'persistent' and args.k == num_mutations:
        args.model = 'persistent_full'
    f_out = check_file_existance('res_exp_{}.csv'.format(args.model))

    if model.status == GRB.Status.OPTIMAL or model.status == GRB.Status.TIME_LIMIT:
        final_error = 0

        clone_used = []
        c = 0
        while c < num_clones:
            s = 0
            us = False
            while s < num_samples:
                if u[s, c].X > 0:
                    us = True
                s += 1
            if us:
                clone_used.append(c)
            c += 1

        uclone_matrix = []
        c = 0
        while c < num_clones:
            if c in clone_used:
                clone = []
                m = 0
                while m < num_mutations:
                    clone.append(int(x[c, m].X))
                    m += 1
                uclone_matrix.append(clone)
            c += 1

        clone_matrix = []
        c = 0
        while c < num_clones:
            clone = []
            m = 0
            while m < num_mutations:
                clone.append(int(x[c, m].X))
                m += 1
            clone_matrix.append(clone)
            c += 1

        usage_matrix = []
        s = 0
        while s < num_samples:
            c = 0
            row = []
            while c < num_clones:
                row.append(u[s, c].X)
                c += 1
            usage_matrix.append(row)
            s += 1

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
    f_out.close()


def print_solution(num_clones, num_samples, num_mutations, start_optimize,
                   args, model, u, x, y, error, input_matrix, matrix_name, mutation_names):
    if args.model == 'dollo' and args.k == 0:
        args.model = 'perfect'
        args.k = ''
    if args.model == 'persistent' and args.k == num_mutations:
        args.model = 'persistent_full'
        args.k = 'full'
    f_out = open('res_{}.txt'.format(matrix_name), 'w+')
    if model.status == GRB.Status.OPTIMAL or model.status == GRB.Status.TIME_LIMIT:
        final_error = 0

        clone_used = []
        c = 0
        while c < num_clones:
            s = 0
            us = False
            while s < num_samples:
                if u[s, c].X > 0:
                    us = True
                s += 1
            if us:
                clone_used.append(c)
            c += 1

        uclone_matrix = []
        c = 0
        while c < num_clones:
            if c in clone_used:
                clone = []
                m = 0
                while m < num_mutations:
                    clone.append(int(x[c, m].X))
                    m += 1
                uclone_matrix.append(clone)
            c += 1

        clone_matrix = []
        c = 0
        while c < num_clones:
            clone = []
            m = 0
            while m < num_mutations:
                clone.append(int(x[c, m].X))
                m += 1
            clone_matrix.append(clone)
            c += 1

        usage_matrix = []
        s = 0
        while s < num_samples:
            c = 0
            row = []
            while c < num_clones:
                row.append(u[s, c].X)
                c += 1
            usage_matrix.append(row)
            s += 1

        if args.model in ('persistent', 'persistent_full'):
            k_model = 2
        elif args.model == 'dollo':
            k_model = args.k + 1
        elif args.model == 'caminsokal':
            k_model = args.k
        else:
            k_model = 1

        extended_matrix = []
        c = 0
        while c < num_clones:
            m = 0
            row = []
            while m < num_mutations*k_model:
                row.append(int(y[c, m].X))
                m += 1
            extended_matrix.append(row)
            c += 1

        error_matrix = []
        s = 0
        while s < num_samples:
            m = 0
            row = []
            while m < num_mutations:
                row.append(error[s, m].X)
                m += 1
            error_matrix.append(row)
            s += 1

        inferred_matrix = np.dot(usage_matrix, clone_matrix)
        diff = input_matrix - inferred_matrix
        for i in range(num_samples):
            for j in range(num_mutations):
                final_error += abs(diff[i][j])

        test_accuracy = LA.norm(diff) / LA.norm(input_matrix)

        f_out.write('Elapsed time: {0} sec\n'.format(
            (datetime.now() - start_optimize).total_seconds()))
        f_out.write('Input name: {0}\nNumber of samples: {1}\n'.format(
            matrix_name, num_samples))
        f_out.write('Number of mutations: {0}\nClone limit: {1}\n'.format(
            num_mutations, args.clones))
        f_out.write('Total clone used: {0}\nModel (k): {1} ({2})\n'.format(
            len(clone_used), args.model, args.k))
        f_out.write('Total error: {0}\nSolution accuracy: {1}\n'.format(
            final_error, test_accuracy))
        f_out.write('Clonal matrix:\n')
        print_lmatrix(uclone_matrix, f_out)
        f_out.write('Usage matrix:\n')
        print_lmatrix(usage_matrix, f_out)
        if not args.model == 'perfect':
            f_out.write('Extended matrix:\n')
            print_lmatrix(extended_matrix, f_out)
        f_out.write('Error matrix:\n')
        print_lmatrix(error_matrix, f_out)
        f_out.write('Inferred matrix F:\n')
        print_lmatrix(inferred_matrix, f_out)
        f_out.write('Tree in DOT code:\n')
        build_tree(np.array(extended_matrix), mutation_names,
                   np.array(usage_matrix), f_out)
    else:
        f_out.write('Elapsed time: {0} sec\n'.format(
            (datetime.now() - start_optimize).total_seconds()))
        f_out.write('Input name: {0}\nNumber of samples: {1}\n'.format(
            matrix_name, num_samples))
        f_out.write('Number of mutations: {0}\nClone limit:{1}\n'.format(
            num_mutations, args.clones))
        f_out.write('Model (k): {1} ({2})\n\nModel unfeasible.'.format(
            args.model, args.k))
    f_out.close()
