# Data and code for the paper: Late Cretaceous origins for major nightshade lineages from total evidence timetree analysis
#### Ixchel González Ramírez, Rocío Deanna, and Stacey D. Smith

***

This repository contains the data and the code necessary to replicate the results and figures of this paper. The respitory contains both RevBayes and R code. Below are details on the organization og the files and instructions to run the analyses.


## Data


## Phylogenetic analyses
All the phylogenetic analyses are set up to be run from the root of this directory calling RevBayes from the command line, e.g. `rb path/to/analysis/header.rev`. We will indicate below the command to run each analysis. The code is written in a modular way, check the `modules/` folder to see all the different variants of each model component and these different model components are called by template files. If you want to replicate our analyses in your own computer, all you have to do is to modify the `header` files. In these, you can specify the directories, and name variants of your output path. 

```
|--headers          # contains the headers to run the analyses
    |--MCMC             # phylogenetic inference
    |--SummaryTrees     # summary trees if needed
    |
|--modules                  # contains all the submodels
    |--analysis                 # different MCMC settings
    |--cont_models              # different continuos traits evolution models
    |--diversification models   # different diversification models
    |--fossilization models     # different fossilization models
    |--mol_clock_models         # different molecular clock models
    |--mol_evol_models          # different molecular evolution models
    |--morph_clock_models       # different morphological clock models
    |--morph_evol_models        # different morphological evolution models
|--output          # by default the output will be written in this folder

```


## R code 


### Assessing convergence


### 'One-fossil at a time' posterior samples inspection




