#/bin/bash
# Start Exp.1 described in the paper
# This experiment could, obviously, be parallelized

clones=( 1 0.8 0.6 0.4 )
for c in "${clones[@]}"
do
	for i in {0..99}
	do
		./gppf -m perfect -f data/simulated/n20_m20/$i.sim -c $c -t 300 -e
		./gppf -m persistent -f data/simulated/n20_m20/$i.sim -k full -c $c -t 300 -e
		./gppf -m dollo -f data/simulated/n20_m20/$i.sim -k 2 -c $c -t 300 -e
		./gppf -m caminsokal -f data/simulated/n20_m20/$i.sim -k 2 -c $c -t 300 -e
	done
done

