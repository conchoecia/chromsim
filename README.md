# chromsim

## Running the program

Make sure that the necessary Python modules are available and execute the following command:

`inversion_sim/main.py [-h] <-C | -c CYCLE_NUMBER> [-w WINDOW_SIZE {1}] [-l LEVEL_OF_CONVERGENCE {1}] [-T] [-t TRANSLOCATIONS_PER_CYCLE {0}] [-o OUTPUT_DIR {"./"}] A_SIZE B_SIZE`

### Parameters

- `-h, --help`: Show detailed help.
- `-C, --converge`: Let the simulation run until all genes are converged (i.e. all possible interactions between any two genes have been explored).
- `-c, --cycle-number`: Let the simulation run for CYCLE_NUMBER cycles, regardless of when convergence is reached. Note that, if CYCLE_NUMBER is low, $\tau_{50\%}$ and the `converged_at` values will return as `-1`.
- `-w, --window-size`: Specify the distance over which two genes are said to be able to interact.
- `-l, --level-of-convergence`: Specify the fraction of possible interactions to be explored when `-C` is set.
- `-T, --plot-average-t50`: If this flag is set, the program will plot the average $\tau_{50\%}$ of previous runs (reading the data from log files), then exit (not running a simulation).
- `-t, --translocations-per-cycle`: Specify the number of genes to randomly move to a different position each cycle (in addition to inversion).
- `-o, --output-dir`: Set the directory to use as the basis for all file-based actions. If not present, `log/` and `diagrams/` subfolders will be created there.
- `A_SIZE, B_SIZE`: Set the number of genes in each ancestral group.
