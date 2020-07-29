# chromEvol_model_adequacy
Model adequacy for likelihood models of chromosome-number evolution.
This framework uses simulations to assess whether a given model of chromosome-number evolution provides a realistic description of the evolutionary process.

# Output:
Supplementary output files will be located in the output directory specified by the user. Model adequacy simulations together with final results will be located under */adequacy_test* sub-directory.
Output files include:
- In the output directory specified by the user:
  - tree_wo_counts.tree: phylogeny without coutns in the tips labels, based on the tree provided by the user. 
  - tree_with_counts.tree: phylogeny with coutns in the tips labels, based on the tree provided by the user. 
  - root_freq: root grequencies as extracted from the results file.
  - increasing_max_chr.txt: while simulating the upper bound was reached. Another iteration of simulations is running. 
  - orig_stats: empirical test statistics calculations.
  - stats_dist_sims: the distributions of each statistics, over *n* simulations.
  - true_percentiles: the percentiles in which the empirical statistics fall within the distributions.
  - percentiles_limits: the 2.5th and 97.5th percentiles per statistics.
  - adequacy_vec: final adequacy vector, with 0 denoting for inadequacy, and 1 denoting for adequacy, per statistic, in the following order: variance, entropy, parisomony score, parsimony-time score.
- A */simulations* directory is also created, which contains:
  - param_sim: simulations parameters (as extracted from the results file).
  - log.txt: log file for simulations.
  - subdirectories contatining the simulations.

# Prerequisites:
1. Python-3.6.5 or newer.
2. Python modules: os, argparse, re, ete3, scipy, numpy, subprocess.
3. R-3.5.1 or newer.
4. R packages: phytools, phangorn.
5. chromEvol: download from [here](https://www.tau.ac.il/~itaymay/cp/chromEvol/chromEvol) or extract from the downloaded chromEvol_model_adequacy repository (located inside the chromEvol folder). Make sure after download or extraction to add the file extension **.ext** to the *chromEvol* file.

# Installation
1. Download all scripts of [chromEvol model adequacy](https://github.com/MayroseLab/chromEvol_model_adequacy) to the same directory.
2. Verify that you have the following paths or input files (see examples in the [example directory](https://github.com/MayroseLab/chromEvol_model_adequacy/tree/master/example)):
- chromEvol executable path: chromevol_path = r"<YOUR_PATH>\chromevol.exe" (for exmplae: C:\Users\Downloads\chromevol.exe)
- R executable path: R_path = r"<YOUR_PATH>\Rscript.exe" (for example: C:\Program Files\R\R-3.6.2\bin\Rscript.exe)
- counts file in FASTA format
- phylogeny file in NEWICK format (the mlAncTree.tree inferred by chromEvol)
- chromEvol results file

# Running
## Modify paths to local executables (required)
In defs.py (definitions file) change the *R_path* and *chromevol_path* variables to **your executable paths**. Full paths are required. See **Installation** for examples.
## The -c <full_path_to_counts_file> parameter: (required)
The empirical data the model in question was fitted to, given in a FASTA format.
## The -t <full_path_to_tree_file> parameter: (required)
The inferred ancestral reconstruction phylogeny in a NEWICK format (by default, named *mlAncestors.tree* by chromEvol).
## The -r <full_path_to_results_file> parameter: (required)
The chromEvol results file from which the simulating parameters will be extracted.
## The -out <full_path_to_output_directory> parameter: (required)
The output directory to which all relevant output files, as well as the simulations and final results will be printed to. Note that in this directory supplementary files will be written, and the model adequacy test itself will be written to a subsequent */adequacy_test* sub-directory.
## The -n <number_of_simulations> parameter:
The number of simulations requested by the user. If not set by the user, the default is 1,000.

# Runnign example:
This can be downloaded from the [example](https://github.com/MayroseLab/chromEvol_model_adequacy/tree/master/example) folder.
*python ~/main.py -c ~/example/counts -t ~/example/tree -r ~/example/chromEvol.res -out ~/example/ -n 1000*
