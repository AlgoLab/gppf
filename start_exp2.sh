#/bin/bash
# Start Exp.2 described in the paper
# This experiment could, obviously, be parallelized

clones=( 1 0.8 0.6 0.4 )
for c in "${clones[@]}"
do
	for i in {0..9}
	do
		./gppf -m perfect -f data/simulated/n25_m50/$i.sim -c $c -t 7200 -e
		./gppf -m persistent -f data/simulated/n25_m50/$i.sim -k full -c $c -t 7200 -e
		./gppf -m persistent -f data/simulated/n25_m50/$i.sim -k 4 -c $c -t 7200 -e
		./gppf -m dollo -f data/simulated/n25_m50/$i.sim -k 4 -c $c -t 7200 -e
		./gppf -m caminsokal -f data/simulated/n25_m50/$i.sim -k 4 -c $c -t 7200 -e
	done
done

