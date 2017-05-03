#!/usr/bin/python
import sys

# The exp argument is used to properly calculate and display output
# Use:
# ./make_table.py 1 			for exp1
# ./make_table.py 2 			for exp2
# All the other configurations will cause the program to exit.


tot_exp = '0'
exp = sys.argv[1]
if exp == '1':
	tot_exp = '100'
elif exp == '2':
	tot_exp = '10'
else:
	sys.exit(0)


def import_csv(path):
	csv = []
	for line in open(path, 'r'):
		csv.append(line[:-1].split(','))
	return(csv)

def get_elem(data, item):
	for i in data:
		if i[0] == item:
			return i
def get_mut_mod(data, mod):
	out = []
	for item in data:
		if item[3] == mod:
			out.append(item)
	return out

perf_csv = import_csv('res_exp_perfect.csv')
pers_csv = import_csv('res_exp_persistent_full.csv')
dollo_csv = import_csv('res_exp_dollo.csv')
caminsokal_csv = import_csv('res_exp_caminsokal.csv')

for mod in ('1.0', '0.8', '0.6', '0.4'):

	perf_m = get_mut_mod(perf_csv, mod)
	pers_m = get_mut_mod(pers_csv, mod)
	dollo_m = get_mut_mod(dollo_csv, mod)
	cs_m = get_mut_mod(caminsokal_csv, mod)

	print('Clone limit: {0}%'.format(int(float(mod)*100)))
	print('\t\t Persistent\t Dollo\t\t Camin-Sokal')
	for perc in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
		countpers = 0
		countdollo = 0
		countcs = 0
		for c in range(int(tot_exp)):
			pf = get_elem(perf_m, str(c))
			ps = get_elem(pers_m, str(c))
			do = get_elem(dollo_m, str(c))
			cs = get_elem(cs_m, str(c))
			if float(ps[8]) <= (1 - perc) * float(pf[8]):
				countpers += 1
			if float(do[8]) <= (1 - perc) * float(pf[8]):
				countdollo += 1
			if float(cs[8]) <= (1 - perc) * float(pf[8]):
				countcs += 1
		if(exp == '1'):
			print('<= {0}\\% PPE \t {1:02d}/{4} \t {2:02d}/{4} \t {3:02d}/{4}'.format(
				int(100 - perc*100),countpers, countdollo, countcs, tot_exp))
		else:
			print('<= {0}\\% PPE \t {1:02d}/{4} \t\t {2:02d}/{4} \t\t {3:02d}/{4}'.format(
				int(100 - perc*100),countpers, countdollo, countcs, tot_exp))
	print('\n')