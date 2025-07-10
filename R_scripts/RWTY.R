########
# RWTY #
########

# this script is for checking the convergence of chains


#install.packages("devtools")
#library(devtools)
#install_github("danlwarren/RWTY")
library(rwty)


# setting procesors number
rwty.processors <<- 8


# read trees for time homo
  # chain 1
my.trees.1 <- load.trees("output/final_runs/time_homogeneous/div_constant_foss_constant_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_1/tree.trees", type = "newick", gens.per.tree = 1)

 # chain 2
my.trees.2 <- load.trees("output/final_runs/time_homogeneous/div_constant_foss_constant_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_2/tree.trees", type = "newick", gens.per.tree = 1)


# running RWTY
solanaceae.rwty <- analyze.rwty(list(my.trees.1, my.trees.2), burnin = 1000)

#ape::unique.multiPhylo(my.trees.1[[1]][1:1000])
#-------------------------------------------------------------------------

# read trees for time hete
# chain 1
my.het.trees.1 <- load.trees("output/final_runs/time_heterogeneous/div_time_heterogeneous_foss_time_heterogeneous_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_1/tree.trees", type = "newick", gens.per.tree = 1)

# chain 2
my.het.trees.2 <- load.trees("output/final_runs/time_heterogeneous/div_time_heterogeneous_foss_time_heterogeneous_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_2/tree.trees", type = "newick", gens.per.tree = 1)


# running RWTY
het.solanaceae.rwty <- analyze.rwty(list(my.het.trees.1, my.het.trees.2), burnin = 1000)


#-------------------------------------------------------------------------
# EXTANT
#-------------------------------------------------------------------------

# read extant trees for time HOMOGENEOUS
# chain 1
ext.my.trees.1 <- load.trees("output/final_runs/time_homogeneous/div_constant_foss_constant_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_1/extant_tree.trees", type = "newick", gens.per.tree = 1)

# chain 2
ext.my.trees.2 <- load.trees("output/final_runs/time_homogeneous/div_constant_foss_constant_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_2/extant_tree.trees", type = "newick", gens.per.tree = 1)


# running RWTY
ext.solanaceae.rwty <- analyze.rwty(list(ext.my.trees.1, ext.my.trees.2), burnin = 1000)

ape::unique.multiPhylo(ext.my.trees.1[[1]][1:2000])



#-------------------------------------------------------------------------

# read extant trees for time heterogeneous
# chain 1
ext.het.my.trees.1 <- load.trees("output/final_runs/time_heterogeneous/div_time_heterogeneous_foss_time_heterogeneous_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_1/extant_tree.trees", type = "newick", gens.per.tree = 1)

# chain 2
ext.het.my.trees.2 <- load.trees("output/final_runs/time_heterogeneous/div_time_heterogeneous_foss_time_heterogeneous_moleclock_UCLN_moleQ_GTR_G_morphclock_linked_morphQ_Mk_G_MCMC_long_2/extant_tree.trees", type = "newick", gens.per.tree = 1)


# running RWTY
ext.het.solanaceae.rwty <- analyze.rwty(list(ext.het.my.trees.1, ext.het.my.trees.2), burnin = 1000)


