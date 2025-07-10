####################
# FIG. 1 in paper
####################


# Read the time heterogeneous only extant MCC tree
   #PATH
extant_heteroTime_MCC_path <- "output/final_runs/time_heterogeneous/MCC_extant_tree_thetero.tre"

  # read
ext_heteroTime_MCC <- readTrees(paths = extant_heteroTime_MCC_path)

#----------------------------------------
# update the names to the last version
#----------------------------------------

# read in names
tb = read.csv("data/conversion_table.csv")


# get tip labels for this tree
current_labels = ext_heteroTime_MCC[[1]][[1]]@phylo$tip.label 

# select relevant names in tb
these_names = tb %>% filter(old_name %in% current_labels)

ext_heteroTime_MCC[[1]][[1]]@phylo <- rename_taxa(ext_heteroTime_MCC[[1]][[1]]@phylo, 
                                                  these_names, 
                                                  key = "old_name", 
                                                  value = "correct_name")

# round the posterior to two digits
ext_heteroTime_MCC[[1]][[1]]@data$posterior <- signif(ext_heteroTime_MCC[[1]][[1]]@data$posterior, digits = 2)

# update tree text
ext_heteroTime_MCC[[1]][[1]]@treetext = write.tree(ext_heteroTime_MCC[[1]][[1]]@phylo)


#----------------------------------------
# Make a plot
#----------------------------------------

# plot the FBD tree for time heterogeneous MCC

het_plot <- plotFBDTree(tree = ext_heteroTime_MCC, 
                        timeline = T, 
                        geo_units = "epochs",
                        tip_labels = T,
                        tip_labels_italics = T,
                        tip_labels_remove_underscore = T,
                        tip_labels_size = 1.3, 
                        line_width = 0.5,
                        tip_age_bars = F,
                        node_age_bars = F,
                        node_pp = T,
                        node_pp_color = c("lemonchiffon", "#440154FF"),
                        node_pp_size = 2,
                        node_labels = "posterior",
                        node_labels_size = 1.4,
                        node_labels_offset = 0.7,
                        node_labels_color = "black",
                        label_sampled_ancs = TRUE,
                        branch_color = "grey25") +
    theme(
      legend.position=c(.05, .6),
      legend.background = element_rect(fill="transparent"),
      #plot.margin = margin(2, 0, 1, 2),        # Top, right, bottom, left
      #axis.text.x = element_blank(),    # Smaller time scale text
      axis.title.x = element_text(size = 8),
      #axis.ticks.length = unit(0.15, "cm")     # Shrinks tick marks
    )


het_plot


# # save to file
# png("figures/Manuscript_figs/Fig1_ext_MCC_tHet.png",
#      width=15,
#      height = 20,
#      res = 300,
#      units = "cm")
# 
# het_plot
# 
# dev.off()
# 
# 
# 
# pdf("figures/Manuscript_figs/Fig1_ext_MCC_tHet.pdf",
#     width=5.9,
#     height =8.2)
# 
# het_plot
# 
# dev.off()

