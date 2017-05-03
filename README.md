General Parsimony Phylogeny models from Frequencies
===================================================

Requirements
-------------
To run `gppf` are required:
- **Python** 2.7.12 +
- **Gurobi** 6.5.2 +

To reproduce experiments described in the paper the following are also required:
- **R** 3.3.3 with packages:
	- *ggplot2*
	- *grid*
	- *plyr*
	- *readr*

For generating samples:
- **ms** by Hudson


gppf
-------
Start the program with:<br>
`gppf [-h] -m {perfect,persistent,dollo,caminsokal} -f FILE [-k K] -t
            TIME -c CLONES [-e]`<br>
Where: 

`-m(--model) {perfect,persistent,dollo,caminsokal}`
Is a required arguments to specify to the phylogent model


`-f(--file) FILE` Specifies the path of the input file.


`-k K` Specifies the **k**-value of the selected model. 
(Persistent(k), Dollo(k),Camin-Sokal(k))
- For a peristent model without restriction, type `-k full`. 
- Do not specify any **k** value for Perfect model, if specified it will be ignored.


`-t(--time) TIME` Specifies the maximum time allowed for the computation. 
Type `-t 0` to not impose a limit.

`-c(--clones) CLONES` Specifies the percentage of clones allowed. Eg `-c 0.8`.
The actual amount of maximum clones used is calculated by : `[clones] * #Mutations`.

`-e(--exp)` Set this parameter to output experimental-format results. See experimental results.

Example of commands:
- `./gppf -m persistent -f data/simulated/n10_m20/21.sim -t 0 -c 0.8 -k 2` 
starts a Persistent(2) model without time limit and with a clone limit of 80%
- `./gppf -m dollo -f data/simulated/n10_m20/21.sim -t 300 -c 1 -k 4 -e`
starts a Dollo(4) model with a time limit of 300 seconds (5 minutes), with a clone limit of 100%
and with a experimental-format output.

Input format
----------------
**gppf** accept to different input format:
- The [simulated data format](data/simulated/n10_m20/0.sim), for file of extension `.sim` (required), is a frequency matrix F of tab separated values. The first line is a dummy line of mutation names.

	| Mut_0	 |	Mut_0- | Mut_1  | Mut_1- |	Mut_2  | Mut_2- |
	|--------|---------|--------|--------|---------|--------|
	| 0.4584 |	0.0000 | 0.0970 | 0.0000 |	0.1630 | 0.1630 |
	| 0.6222 |	0.2443 | 0.1420	| 0.2443 |	0.1073 | 0.1073 |
	| 0.6450 |	0.1400 | 0.1539	| 0.1400 |	0.1270 | 0.1270 |
	| 0.6110 |	0.0805 | 0.0838	| 0.0805 |	0.0794 | 0.0794 |
	| 0.6930 |	0.1108 | 0.0920	| 0.1108 |	0.0792 | 0.0792 |

- The [real data format](data/real/CLL077_deep.tx) is a format were the input is a tab-separated text file. The first line contains the sample names. The first column contains mutations ids. The consecutive pair of columns contains read counts
for reference and mutated alleles.

Detailed-format Output
-----------------
When not running an experiment **gppf** outputs all the information in detail in a file called `res_INPUTNAME.txt` created in the
same folder where gppf is running. The output file contains:
- Elapsed time
- Input name
- Number of samples
- Number of mutations
- Clone limit
- Total clone used
- Model (k)
- Total error
- Solution accuracy
- Clonal matrix
- Usage matrix
- Extended matrix
- Error matrix
- Inferred matrix F
- Tree in DOT code

We can see as an example the output of [CLL077](experimental_results/real/res_cll077_deep.txt),
to run the execution of *CLL077* as described in the paper use:

`./gppf -m persistent -f data/real/cll077_deep.txt -t 0 -c 1 -k 2`


Experimental-format Output
------------------------------
In the experiments the programs that does not outputs all the
informations regarding the clonal matrix, the usage matrix and does not print
a tree of the reconstructed phylogeny.
It instead prints informations useful for testing. In specific it prints a 
**CSV** file with the following informations: <br>
`matrix_name, num_samples, num_mutations, mut_mod, 
clone_used, k, time, total_error, accuracy`

We can see as an example the Dollo(2) output of [Exp.2](experimental_results/simulated/exp2/dollo.csv)

Recreating experiments
-----------------------
To start the experiments described in the paper run the bash files `start_exp1.sh` and `start_exp2.sh`
note that this experiments can be easily parallelized by separating the `for` cycles in different files
or sessions. The bash files also recreate the plot and the table present in the paper by running programs
`plot_from_csv.R` and `make_table.py`, the last two programs require the output file of `gppf` to be in 
the root directory, (as default)



