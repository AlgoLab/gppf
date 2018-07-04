def import_ternary_matrix_tab(file_path):
	with open(file_path, 'r') as fin:
		m = []

		for line in fin:
			row = [int(x) for x in line.strip().split()]
			m.append(row)
			
		return m

def import_matrix_tab(file_path, type, k_model):
	if '.sim' in file_path:
		with open(file_path, 'r') as f_in:
			first_line = f_in.readline()[:-2]
			names = list(first_line.split('\t'))

			m = []
			for line in f_in:
				row = list(line[:-1].split('\t'))
				m.append(map(lambda x: float(x),row))

			return m, names

	else:
		m = []
		fo = open(file_path, 'r')
		first_line = fo.readline()
		names = []
		num_s = (len(first_line.split('\t')) - 1) / 2
		for k in range(num_s):
			m.append([])

		for line in fo:
			s = []
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
					x = float(row[i+1])/ (int(row[i+1]) + int(row[i]))
					m[i/2].append(x)
				else:
					m[i/2].append(0)
				i +=2
		fo.close()
		return m, names

def print_lmatrix(lmatrix, file_out=None):
	for i in lmatrix:
		row = ''
		for c in i:
			row += str(c) + ' '
		if file_out:
			file_out.write(row + '\n')
		else:
			print(row)
	
