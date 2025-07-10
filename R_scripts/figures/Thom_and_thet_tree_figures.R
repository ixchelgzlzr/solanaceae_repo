################################
# Using revgadgets for figures #
################################

# Side by side Thom and Thet trees


#install.packages("devtools")
#devtools::install_github("cmt2/RevGadgets")

library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)
library(ape)


# specify files paths
homo_time_trace    <- "output/final_runs/time_homogeneous/params_combined.log"
homo_time_trees    <- "output/final_runs/time_homogeneous/trees_combined.trees"
hetero_times_trace <- "output/final_runs/time_heterogeneous/params_combined.log"
hetero_times_trees <- "output/final_runs/time_heterogeneous/trees_combined.trees"

#-----------------------
# Making a tree figure
#-----------------------

# time homogenerous all
#homoTime_MCC_path <- "output/final_runs/time_homogeneous/MCC_tree_thomo.tre"
#homoTime_MAP_path <- "output/final_runs/time_homogeneous/MAP_tree_thomo.tre"

# time homogenerous extant
extant_homoTime_MCC_path <- "output/final_runs/time_homogeneous/MCC_extant_tree_thomo.tre"
#extant_homoTime_MAP_path <- "output/final_runs/time_homogeneous/MAP_extant_tree_thomo.tre"

# time heterogeneous all
#heteroTime_MCC_path <- "output/final_runs/time_heterogeneous/MCC_tree_thetero.tre"
#heteroTime_MAP_path <- "output/final_runs/time_heterogeneous/MAP_tree_thetero.tre"

# time heterogeneous extant
extant_heteroTime_MCC_path <- "output/final_runs/time_heterogeneous/MCC_extant_tree_thetero.tre"
#extant_heteroTime_MAP_path <- "output/final_runs/time_heterogeneous/MAP_extant_tree_thetero.tre"



# read in the trees
  # time homo
#homoTime_MCC <- readTrees(paths = homoTime_MCC_path)
#homoTime_MAP <- readTrees(paths = homoTime_MAP_path)
ext_homoTime_MCC <- readTrees(paths = extant_homoTime_MCC_path)
#ext_homoTime_MAP <- readTrees(paths = extant_homoTime_MAP_path)


  # time hetero
#heteroTime_MCC <- readTrees(paths = heteroTime_MCC_path)
#heteroTime_MAP <- readTrees(paths = heteroTime_MAP_path)
ext_heteroTime_MCC <- readTrees(paths = extant_heteroTime_MCC_path)
#ext_heteroTime_MAP <- readTrees(paths = extant_heteroTime_MAP_path)


#----
# Let's change the names of the tips to the names in Rocio's paper
# read in names
tb = read.csv("data/conversion_table.csv")


# get tip labels for this tree
current_labels   = ext_heteroTime_MCC[[1]][[1]]@phylo$tip.label 
current_labels_2 = ext_homoTime_MCC[[1]][[1]]@phylo$tip.label

# select relevant names in tb
these_names   = tb %>% filter(old_name %in% current_labels)
these_names_2 = tb %>% filter(old_name %in% current_labels_2)

ext_heteroTime_MCC[[1]][[1]]@phylo <- rename_taxa(ext_heteroTime_MCC[[1]][[1]]@phylo, 
                                                  these_names, 
                                                  key = "old_name", 
                                                  value = "correct_name")

ext_heteroTime_MCC[[1]][[1]]@treetext <- write.tree(ext_heteroTime_MCC[[1]][[1]]@phylo)


ext_homoTime_MCC[[1]][[1]]@phylo <- rename_taxa(ext_homoTime_MCC[[1]][[1]]@phylo, 
                                                these_names_2, 
                                                key = "old_name", 
                                                value = "correct_name")

ext_homoTime_MCC[[1]][[1]]@treetext <- write.tree(ext_homoTime_MCC[[1]][[1]]@phylo)

#----


# plot the FBD tree for time homogeneous MCC
homo_plot <- plotFBDTree(tree = ext_homoTime_MCC, 
            timeline = T, 
            geo_units = "epochs",
            tip_labels = T,
            tip_labels_italics = F,
            tip_labels_remove_underscore = T,
            tip_labels_size = 1.6, 
            line_width = 0.8,
            tip_age_bars = T,
            node_age_bars = T, 
            age_bars_colored_by = "posterior",
            label_sampled_ancs = TRUE) +
  # use ggplot2 to move the legend and make 
  # the legend background transparent
  theme(legend.position=c(.05, .6),
        legend.background = element_rect(fill="transparent")) 
  #ggtitle("Time homogeneous model - MCC tree - ONLY extant")

homo_plot





# plot the FBD tree for time heterogeneous MCC
het_plot <- plotFBDTree(tree = ext_heteroTime_MCC, 
                         timeline = T, 
                         geo_units = "epochs",
                         tip_labels = T,
                         tip_labels_italics = F,
                         tip_labels_remove_underscore = T,
                         tip_labels_size = 1.6, 
                         line_width = 0.8,
                         tip_age_bars = T,
                         node_age_bars = T, 
                         age_bars_colored_by = "posterior",
                         label_sampled_ancs = TRUE) +
  # use ggplot2 to move the legend and make 
  # the legend background transparent
  theme(legend.position=c(.05, .6),
        legend.background = element_rect(fill="transparent")) 
#ggtitle("Time homogeneous model - MCC tree - ONLY extant")

het_plot




# sice by side
grid.arrange(homo_plot,het_plot, ncol = 2)

# # export figure
# tiff("figures/Manuscript_figs/extant_time_homo_VS_het.tiff",
#      width=17,
#      height=23,
#      res=300,
#      units = "cm")
# 
# ggarrange(homo_plot,het_plot,
#           ncol = 2,
#           common.legend = T,
#           legend = "bottom",
#           #labels = "AUTO"
#           labels = list("A) Time-homogeneous FBD", "B) Time-heterogeneous FBD"),
#           font.label = list(size= 11)
#              )
# 
# dev.off()



#------------------------------
# cophyloplot to see differences in topology
#------------------------------

# comparePhylo(ext_homoTime_MCC[[1]][[1]]@phylo, ext_heteroTime_MCC[[1]][[1]]@phylo, 
#              plot = T, 
#              commons = F,
#              force.rooted = T,
#              cex = 0.5, 
#              type = "fan",
#              #use.edge.length = T,
#              location = "bottom")


#creation of the association matrix:


x <- phytools::cophylo(ext_homoTime_MCC[[1]][[1]]@phylo, 
        ext_heteroTime_MCC[[1]][[1]]@phylo, rotate = F)

plot.cophylo(x, fsize=c(0.3,0.3), 
             link.type = "curved",
             link.lty = 'solid', 
             link.col = "#31688EFF",
             link.lwd=2)


# # export figure
# tiff("figures/Manuscript_figs/co_phylo_plot.tiff",
#      width=15,
#      height=20,
#      res=300,
#      units = "cm")
# 
# plot.cophylo(x, fsize=c(0.3,0.3), 
#              link.type = "curved",
#              link.lty = 'solid', 
#              link.col = "#31688EFF",
#              link.lwd=1.5)
# 
# dev.off()




