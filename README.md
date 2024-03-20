# chromsim

## Usage

### General-purpose parameters

A handful of arguments apply to all or most uses for this program. These are:

- `-o/--output-dir`: the directory in which any produced files are saved
- `-f/--filename`: the name (without a file ending) of the produced file in cases where only one output file is produced

### Running a simulation

To simulate a series of random inversions, run the script `main.py` with the argument `-S/--simulate`. This requires the additional arguments `-a/--asize A` (`A` is the number of genes in the first linkage group), `-b/--bsize B` (same thing for the second linkage group) in any case. This produces a `.inv` file containing all the information needed to reproduce the run (for instance for plotting it, see below). A result might look like this:

```
# timestamp 2024-03-18 15:54:46.171033
## params
# Asize 100
# Bsize 95
# window_size 1
# level_of_convergence 1.00
## results
# t100 146669
# t50 6489
# tS 137
# mu 0.99
# sigma 0.07
# AB_convergence 146669

6	90
114	155
34	106
59	142
45	118
.
.
.
```

The recorded results are:

- τ<sub>100</sub>: the number of inversion cycles until convergence (i.e. the total number of cycles in a "normal" run with no meddling options)
- τ<sub>50</sub>: the number of cycles by which 50% of all possible gene pairings have been recorded
- τ<sub>S</sub>: the number of cycles until the first entropy value within the 1.96σ of μ was recorded
- mu (μ): the arithmetic mean of the entropy values at each cycle (after a 25% burn-in)
- sigma (σ): the standard deviation of the entropy values at each cycle (after a 25% burn-in)
- AB_convergence: the number of cycles it took for each gene in group A to have interacted with each gene in group B and vice versa (though not necessarily with each gene in its own group)

followed by a list of positions of the cuts made for each inversion. For each inversion, two positions on the chromosomes are randomly chosen, and the sequence of genes between them gets inverted.

### Optional parameters
 
If you want to specify a window size (the distance along the chromosome to the left and right of a gene that are counted towards its new interactions after every iteration), you can do so by using `-w/--window-size W`. Otherwise, a default value of `W=1` will be used. If `-n/--cycle-number N` is used, only `N` inversion cycles will be simulated before exiting. Otherwise the simulation will run until all genes have converged (been within `W` of every other gene at least once over the course of the simulation). Similarly, by using `-l/--level-of-convergence L`, the simulation will stop once a fraction of all possible unique interactions between two genes equal to `L` has been recorded (`L` has to be within [0, 1], of course). `L=1` equates to omitting this option entirely.

### Using a `.rbh` file as a basis
 
 chromsim accepts a `.rbh` file as input using the option `-r/--rbh path/to/file.rbh`. Instead of passing the A and B sizes using the above options, you will need to specify the linkage groups' names using `-A/--group-a ALG_A` and `-B/--group-b ALG_B`, respectively. Additionally, the program expects the name of the chromosome or scaffold with `-c/--chromosome CHROM`. The optional parameters above still apply. In addition, you have the option to only run the simulation until the entropy value calculated from the gene positions in the `.rbh` file has been reached. To do so, add `-M/--find-m` to the execution command.
