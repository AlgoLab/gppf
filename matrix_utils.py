import re

def import_matrix_tab(file_path, type, k_model):
    if '.sim' in file_path:
        with open(file_path, 'r') as f_in:
            first_line = f_in.readline()[:-2]
            names = list(first_line.split('\t'))

            m = []
            for line in f_in:
                row = list(line[:-1].split('\t'))
                m.append(map(lambda x: float(x), row))

            return m, names

    else:
        fo = open(file_path, 'r')
        first_line = fo.readline()
        names = []
        num_s = (len(first_line.split('\t')) - 1) / 2
        m = [[] for k in range(num_s)]

        for line in fo:
            row = line[:-1].split('\t')[1:]
            mut = line[:-1].split('\t')[0]
            if ',' in mut:
                mut = mut.split(',')[0]
            names.append(mut)
            if type == 'persistent':
                names.append(mut + '-')
            elif type == 'dollo':
                for i in range(k_model):
                    names.append(mut + '-{0}'.format(i+1))
            elif type == 'caminsokal':
                for i in range(k_model-1):
                    names.append(mut + '+{0}'.format(i+1))

            i = 0
            while i < len(row):
                if not row[i+1] == '':
                    x = float(row[i+1]) / (int(row[i+1]) + int(row[i]))
                    m[i/2].append(x)
                else:
                    m[i/2].append(0)
                i += 2
        fo.close()
        return m, names

def read_matrix_tab(file_path):
    def parse_line(array):
        return [int(x) for x in array]

    if '.sim' in file_path:
        with open(file_path, 'r') as f_in:
            first_line = f_in.readline()[:-2]
            names = list(first_line.split('\t'))

            m = []
            for line in f_in:
                row = list(line[:-1].split('\t'))
                m.append(map(lambda x: float(x), row))

            return m, names

    else:
        with open(file_path, 'r') as fo:
            lines = fo.readlines()
            matrix = [parse_line(re.split("\s", line.rstrip("\n"))) for line in lines]

            return matrix

def expand_name(s, max_gains, max_losses):
    positive_names = [s + '+' + str(i) for i in range(max_gains)]
    negative_names = [s + '-' + str(i) for i in range(max_losses)]
    return positive_names + negative_names

def compute_names(matrix, max_gains, max_losses):
    names = [str(i) for i in range(len(matrix[0]))]
    expanded = [expand_name(name) for name in names]

    return expanded


def print_lmatrix(lmatrix, file_out=None):
    for i in lmatrix:
        row = ''
        for c in i:
            row += str(c) + ' '
        if file_out:
            file_out.write(row + '\n')
        else:
            print(row)
