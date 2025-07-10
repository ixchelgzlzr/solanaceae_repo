#########################################################
# Time heterogeneous tree with not so uncertain fossils
# Making a figure
#########################################################


# path to trimmed time heterogeneous MCC
PATH = "output/final_runs/time_heterogeneous/unstable_fossils_trimmed/MCC_fossil_trimmed_tree.tree"

#PATH = "output/final_runs/time_heterogeneous/unstable_fossils_trimmed_6fos/MCC_fossil_trimmed_tree_6f.tree"

  
# read the tree
tree = readTrees(PATH, tree_name = "timetree")


# Let's change the names of the tips to the names in Rocio's paper

# read in names
tb = read.csv("data/conversion_table.csv")


# get tip labels for this tree
current_labels = tree[[1]][[1]]@phylo$tip.label 

# select relevant names in tb
these_names = tb %>% filter(old_name %in% current_labels)

tree[[1]][[1]]@phylo <- rename_taxa(tree[[1]][[1]]@phylo, 
                                                  these_names, 
                                                  key = "old_name", 
                                                  value = "correct_name")

tree[[1]][[1]]@treetext <- write.tree(tree[[1]][[1]]@phylo)



# Get the fossils for coloring
tips <- tree[[1]][[1]]@phylo$tip.label


fossils <- c("Physalis_hunickenii", "Physalis_infinemundi", "Solanoides_dorofeevii", "Eophysaloides_inflata",
             "Eophysaloides_inflata", "Solanum_miocenicum", "Hyoscyosperma_daturoides", "Hyoscyamus_undulatus", 
             "Thanatosperma_minutum", "Lycianthoides_calycina")


tips <- data.frame(tips)

tips$fossil <- NA

for(i in 1:nrow(tips)){
  
  if (tips$tips[i] %in% fossils){
    tips$fossil[i] <- "#CD4071FF"
  } else {
    tips$fossil[i] <- "black"
  }
  
}


tree.plot <- plotFBDTree(tree = tree, 
                        timeline = T, 
                        geo_units = "epochs",
                        tip_labels = F,
                        tip_labels_italics = T,
                        tip_labels_remove_underscore = T,
                        tip_labels_size = 2, 
                        #tip_labels_color = as.factor(tips$fossil),
                        tip_age_bars = F,
                        age_bars_color = alpha("#26828EFF", 0.5),
                        age_bars_width = 1,
                        node_age_bars = T, 
                        node_pp = T,
                        node_pp_color = c("lemonchiffon", "#440154FF"),
                        #node_pp_shape =19,
                        node_pp_size = 2.5,
                        #age_bars_colored_by = "posterior",
                        label_sampled_ancs = T,
                        line_width = 0.5) + 
  # use ggplot2 to move the legend and make 
  # the legend background transparent
  theme(legend.position=c(.05, .6),
        legend.background = element_rect(fill="transparent")) 



tree.plot + geom_tiplab(color = factor(tips$fossil), size = 2)



# # exporting figure
# pdf("figures/fossil_trimmed_MCC_tree_2.pdf", width=8, height=14)
# 
# tree.plot %<+% tips + geom_tiplab(color = factor(tips$fossil), size = 2, fontface =3)
# 
# dev.off()
# 
# 
# svg("figures/fossil_trimmed_MCC_tree_2.tiff", width=8, height=14, units ="in", res = 300)
# 
# tree.plot + geom_tiplab(color = factor(tips$fossil), size = 2, fontface =3)
# 
# dev.off()

