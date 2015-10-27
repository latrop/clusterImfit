# clusterImfit
Run imfit in genetic mode on a cluster

# Algorithm
At first the searching of the best model is being performed via genetic algorithm in
parallel mode. When the GA is converged, then Levenberg-Marquardt optimisation is applied
to N best organisms of the last generation (where N is the number of nodes). Resulting model
is the one with best chi squared among results of LM-optimisation.

# GA stopping conditions
1) Maximum number of iterations is reached (**maxGenNumber** parameter in config).

2) Relative improvement in both best model and average generation fitness is less
then **fTol** for at least **fSpan** last generations.

# Usage
1) Place in the program directory a valid imfit input file

2) Place in the program directory fits files (object and (optional) psf and mask)

3) Configure program with config file. All parameter names are explained in the example of the config.

4) Run cluster_imfit.py without arguments.

Results will be in 'results' directory. 'results/generations' directory will contain one best
organism per generation, so one can see the progress of the optimisation.