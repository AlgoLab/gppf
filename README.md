General Parsimony Phylogeny models from Frequencies
===================================================

REFERENCE TO THE PAPER

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
- The [simulated format](data/simulated/n10_m20/0.sim), for file of extension `.sim` (required), is a frequency matrix F of tab separated values. The first line is a dummy line of mutation names.

	| Mut_0		Mut_0-	Mut_1	Mut_1-	Mut_2	Mut_2- |
	---------------------------------------------------
	| 0.3849	0.0000	0.3849	0.3271	0.3271	0.1661 |	
	| 0.4577	0.0000	0.4577	0.2598	0.2598	0.0000 |	
	| 0.3904	0.0000	0.3904	0.2582	0.2582	0.1561 |	
	| 0.1416	0.0000	0.1416	0.1416	0.1416	0.0000 |	
	| 0.5114	0.0000	0.5114	0.3046	0.3046	0.1853 |	
	

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

To start the experiments described if the paper run the bash files `start_exp1.sh` and `start_exp2.sh`
note that this experiments can be easily parallelized by separating the `for` cycles in different files
or sessions. 




