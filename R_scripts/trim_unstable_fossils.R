library(ape)
library(phangorn)
library(rlist)

# read trace
tree.trace.hetero = read.table("output/final_runs/time_heterogeneous/trees_combined.trees", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)



# fossils to trimm
# fossils_to_trimm <- c("CCN6065_second_seed_Capsicum_pliocenicum",
#                       "K587B_Solanum_foveolatum",
#                       "UF6500_Nephrosemen_reticulatus",
#                       "V_40891_Solanispermum_reniforme",
#                       "V_40898_Solanum_arnense")


# adding Solanum miocenicum
fossils_to_trimm <- c("CCN6065_second_seed_Capsicum_pliocenicum",
                      "K587B_Solanum_foveolatum",
                      "UF6500_Nephrosemen_reticulatus",
                      "V_40891_Solanispermum_reniforme",
                      "V_40898_Solanum_arnense",
                      "K528_Solanum_miocenicus")




for (r in 1:nrow(tree.trace.hetero)){
  
  tree = ape::read.tree(text = tree.trace.hetero$timetree[r])
  new_tree = drop.tip(tree, fossils_to_trimm)
  
  #substitute
  tree.trace.hetero$timetree[r] = as.character(write.tree(new_tree))
  
}


# # write combined traces to file
# write.table(tree.trace.hetero, file=("output/final_runs/time_heterogeneous/unstable_fossils_trimmed/fossil_trimmed_trace.tree"), row.names=FALSE, sep="\t", quote=FALSE)

#write.table(tree.trace.hetero, file=("output/final_runs/time_heterogeneous/unstable_fossils_trimmed/fossil_trimmed_trace_1.tree"), row.names=FALSE, sep="\t", quote=FALSE)



