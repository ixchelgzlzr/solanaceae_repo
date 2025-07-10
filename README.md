# Data and code for the paper: Late Cretaceous origins for major nightshade lineages from total evidence timetree analysis
#### Ixchel González Ramírez, Rocío Deanna, and Stacey D. Smith

***

This repository contains the data and the code necessary to replicate the results and figures of this paper. The respitory contains both RevBayes and R code. Below are details on the organization og the files and instructions to run the analyses.


## Data

All the data necessary for the analyses are in this folder. The folder `alignments/` has the 10 markers used. The file `Ages.txt` has a list of all the tips with their min and max age. For extant species, this are both 0, while for fossils these intervals represent the stratigraphic uncertainty of the fossils. The file `clade_constraints.Rev` contains the definition of the two clade constraints, and the file `R_backbone_from_sampled_3_tips.tree` contains the three tips backbone used for the TED analyses. The files `continuous.nexus` and `morphology.nexus` contain the continuous and categorical morphological data respectively. Finally, `epoch_timescale.csv` defines the time intervals in which rates are allowed to vary for the time-heterogeneous TED analysis.

**Important note:** The names of some tips were updated by the time we had concluded our analyses. We use the correspondence table `conversion_table.csv` to update the names in the final figures.

## Phylogenetic analyses
All the phylogenetic analyses are set up to be run from the root of this directory calling RevBayes from the command line, e.g. `rb path/to/analysis/header.rev`. We will indicate below the command to run each analysis. The code is written in a modular way, check the `modules/` folder to see all the different variants of each model component and these different model components are called by template files. If you want to replicate our analyses in your own computer, all you have to do is to modify the `headers/` files. In these, you can specify the director for your output path by specifying `output_dir` and `output_extra`.


```
# Folder structure

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
To run the preliminary analyses that only infer phylogenetic relationships between all the tips (fossils and extant) without inferring divergence times, you can run: 
* `rb headers/MCMC/only_morph_no_continuous.Rev`, to use only categorical morphological traits.
* `rb headers/MCMC/only_morph.Rev` , to use categorical and continuous morphological traits.
* `rb headers/MCMC/unrooted_all_data_mole_scaler.Rev`, to use categorical and continuous morphological traits, as well as molecular data.

To run the total evidence dating analyses to infer both divergence times and phylogenetic relationships, you can run:
* `rb headers/MCMC/model_TED_w_backbone_const.Rev`, to run the **time homogeneous** TED model.
* `rb headers/MCMC/model_TED_w_backbone_const_timehet.Rev`, to run the **time heterogeneous** TED model. 


**IMPORTANT NOTES**

* If you are running multiple MCMCs we recommend to manually run them and modify the `output_extra  =` argument of each header file. Not changing this value will overwrite your results.
* We discovered a back end bug in rb. When the first drawn tree has -inf probability, and rb tries to get a new starting tree it misrepresents the age of fossils. While hopefully this will be fixed soon, the workaround for now is **to stop the analysis `ctrl +C` and restart until the analysis finds a good starting tree in the first attempt**. To be more specific, if you see in screen something like the chunk below, you need to stop and re-start.

```
Processing of file "modules/morph_evol_models/Mk_G.Rev" completed
   Processing file "modules/cont_models/brownian_unlinked_global.Rev"
   Processing of file "modules/cont_models/brownian_unlinked_global.Rev" completed
   Processing file "modules/analysis/MCMC_long.Rev"
   Could not compute lnProb for node 'timetree': lnProb = -inf
   (((Brunfelsia_australis[&index=148]:9.906192,Anisod...
   
   
   Drawing new initial states ... 
   Could not compute lnProb for node 'timetree': lnProb = -inf
   (((((H4895_Hyoscyamus_undulatus[&index=148]:36.9771...
```  
If the first attempt to get a tree was succesful, you will see something like this:

```
   Processing of file "modules/morph_evol_models/Mk_G.Rev" completed
   Processing file "modules/cont_models/brownian_unlinked_global.Rev"
   Processing of file "modules/cont_models/brownian_unlinked_global.Rev" completed
   Processing file "modules/analysis/MCMC_long.Rev"
   
   Running burn-in phase of Monte Carlo sampler for 200 iterations.
   This simulation runs 1 independent replicate.
   The simulator uses 381 different moves in a random move schedule with 4358.6 moves per iteration
   

Progress:
0---------------25---------------50---------------75--------------100
```

## R code 


### Assessing convergence


### 'One-fossil at a time' posterior samples inspection




