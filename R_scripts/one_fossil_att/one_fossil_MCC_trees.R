##############################################################################
# I am going to prune the fossils so we can make MCC trees for individual fossils
#################################################################

library(ape)
library(phangorn)
library(rlist)


# Read trace without burnin
tree.trace <- read.tree("output/final_runs/time_heterogeneous/trees_combined.trees")


# all the fossils
fossil_names <- c("2003_0073_Nicotianosemen_minuta", 
                  "CCN6065_second_seed_Capsicum_pliocenicum",
                  "f_Eophysaloides_inflata",
                  "f_Lycianthoides_calycina",
                  "f_Physalis_hunickenii",
                  "f_Physalis_infinemundi",
                  "H4895_Hyoscyamus_undulatus",
                  "K453A_Hyoscyamus_datureoides",
                  "K528_Solanum_miocenicus",
                  "K530_Solanoides_dorofeevii",
                  "K587B_Solanum_foveolatum",
                  "UF6500_Nephrosemen_reticulatus",
                  "V_40891_Solanispermum_reniforme",
                  "V_40898_Solanum_arnense")


# number of trees to use for summaries
N_trees = 13900

# make an empty list of MCC trees
MCC_trees_trimmed <- vector(mode = "list", length = length(fossil_names))


# for each fossil
for ( i in 1:length(fossil_names)){
  
  # get the name of the fossil
  this_fossil <- fossil_names[i]
  
  # Get the vector of fossils to drop from sample
  names_to_drop <- setdiff(fossil_names, this_fossil)
  
  # make empty list
  trimmed_trace <- vector(mode = "list", length = length(tree.trace))
  
  # drop tips from the trace 
  for ( t in 1:length(tree.trace)){
    
    trimmed_trace[[t]] <- drop.tip(tree.trace[[t]], names_to_drop)
    
  }
  

  # sample N trees
  sample_trees <- sample(trimmed_trace, N_trees)
  
  # make MCC tree
  #MCC_trees_trimmed[[i]] <- maxCladeCred(sample_trees, part = prop.part(sample_trees))
  MCC_trees_trimmed[[i]] <- maxCladeCred(sample_trees)
  
  
}




# for (f in 1:length(fossil_names)){
#   
#   cols        <- rep("black", length(MCC_trees_trimmed[[f]]$tip.label))
#   index       <- which(MCC_trees_trimmed[[f]]$tip.label == fossil_names[f])
#   cols[index] <- "red"
#   
#   pdf(paste0("figures/Manuscript_figs/MCC_one_fossil_att_TimeHete/", fossil_names[f], ".pdf"), width=8, height=12)
#   
#   plot(MCC_trees_trimmed[[f]], cex = 0.7, tip.color  = cols )
#   nodelabels(round(MCC_trees_trimmed[[f]]$node.label, digits = 2), frame = "none", col = "blue", cex = 0.60, adj = c(1, 1),)
#   title(fossil_names[f])
#   
#   dev.off()
#   
# }
