#################################
# Getting ages from time het ext
##################################

library(phytools)
library(dispRity)
library(RevGadgets)
library(dplyr)
library(gridExtra)
library(ggplot2)


# tree path
tree_path = "output/final_runs/time_heterogeneous/MCC_extant_tree_thetero.tre"


# Read the trees
tree = RevGadgets::readTrees(tree_path)

ntips = length(tree[[1]][[1]]@phylo$tip.label)

# Extract the 95 HPD age info
hpd <- tree[[1]][[1]]@data

hpd <- hpd[ , c("age_0.95_HPD", "node")]
hpd <- hpd %>% arrange(as.numeric(node))
hpd <- hpd[(ntips+1):nrow(hpd), ]


# get the ages
ages = dispRity::tree.age(tree[[1]][[1]]@phylo)
ages = ages[(ntips+1):nrow(ages), ]

# combine the tables
table <- full_join(ages, hpd, by= c("elements"="node"))


# the table is in list format, so I need to make it a data
# frame to export:
node <- vector(length = length(table$ages))
age <- vector(length = length(table$ages))
minHPD <- vector(length = length(table$ages))
maxHPD <- vector(length = length(table$ages))


for (i in 1:length(table$ages)){
  node[i]    <- table$elements[[i]]
  age[i]      <- table$ages[[i]]
  minHPD[i]  <- table$age_0.95_HPD[[i]][[1]]
  maxHPD[i]  <- table$age_0.95_HPD[[i]][[2]]
}

ages_table <- data.frame(node, age, minHPD, maxHPD)
rownames(ages_table) <- NULL

# save file
write.csv(ages_table, "R_output/ages_table.csv")



#plot phylogeny with nodes
tree_fig <- RevGadgets::plotTree(tree, 
                                 tip_labels_italics = T,
                                 tip_labels_remove_underscore = T,
                                 tip_labels_size = 2,
                                 node_labels = "node", 
                                 node_labels_size = 2,
                                 node_labels_offset = 1, 
                                 node_labels_color = "#440154FF",
                                 branch_color  = "gray50",
                                 line_width = 0.75, 
                                 timeline =T, 
                                 geo_units = "epochs")


my_table <- tableGrob(ages_table, theme = ttheme_default(base_size = 6, padding = unit(c(2.5, 2.5), "mm")) , rows = NULL)


#grid.arrange(arrangeGrob(tree_fig, my_table, ncol=2), respect = T)

# pdf(file="figures/tree_for_ages.pdf",  width=6, height=12)
# 
# tree_fig
# 
# dev.off()

