########################################################################
# FOSSIL POSITIONS ONE FOSSIL AT THE TIME
#
# NOTE TO SELF: this can be done including all fossils simultaneously, 
# but for consistency we will report these resulst in the manuscript. 
# Including fossils can be run with the script "fossil_position.R"
# Inferences do not change!
#########################################################################



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

tree.trace  <- read.tree("output/final_runs/time_heterogeneous/trees_combined.trees")
sample_size <- 13900
trees <- tree.trace[1:sample_size]

#------------------------------------------
# Extract fossils and its sister clades
#------------------------------------------

# make empty container
trimmed_trace  = vector(mode = "list", length = length(fossil_names))

# for each fossil
for (i in 1:length(fossil_names)){
  
  # get this fossil name
  this_fossil = fossil_names[[i]]
  
  # Get the vector of fossils to drop from sample
  names_to_drop <- setdiff(fossil_names, this_fossil)
  
  # make empty list
  #trimmed_trace[[i]] <- vector(mode = "list", length = length(trees))
  
  # drop tips from the trace 
  for (t in 1:length(trees)){
    
    trimmed_trace[[i]][[t]] <- drop.tip(trees[[t]], names_to_drop)
    
  }
}
    
  
# empty list
fossil_and_sister = vector(mode = "list", length = length(fossil_names))
  
# for each fossil
for (f in 1:length(trimmed_trace)){
  
  # get this fossil name
  this_fossil = fossil_names[[f]]
  
  # for each tree of the trace
  for (t in 1:length(trimmed_trace[[f]])){
    
    # get the tip number for the fossil
    tip_number = which(trimmed_trace[[f]][[t]]$tip.label == this_fossil)
    
    # get the node number for that tip by using the edges
    this_node_row    = which(trimmed_trace[[f]][[t]]$edge[ , 2] == tip_number)
    this_node_number = trimmed_trace[[f]][[t]]$edge[this_node_row, 1] 
    
    # extract the clade originating from that node
    this_clade = extract.clade(trimmed_trace[[f]][[t]], node = this_node_number)
     
    # save fossil + sister in list
    fossil_and_sister[[f]][[t]] = this_clade
    
  } }




#------------------------------------------
# Summarize
#------------------------------------------


# get the tip names per each tree

# empty contained
tips_fossil_and_sister = vector(mode = "list", length = length(fossil_and_sister))

# get tips

# for each fossil
for (f in 1:length(fossil_and_sister)){
  
  # for each tree
  for(t in 1:length(trees)){
    
    # get the tip labels
    tips_fossil_and_sister[[f]][[t]] = fossil_and_sister[[f]][[t]]$tip.label
    
  }
}



# get the unique clades to which each fossil is sister to
summary_sis_groups = vector("list", length = length(tips_fossil_and_sister))


# for each fossil
for (i in 1:length(tips_fossil_and_sister)){
  
  # unlist the tips names
  these_clades  = lapply(tips_fossil_and_sister[[i]], unlist)
  these_clades  = lapply(these_clades, sort)
  # extract a df with the sister clades for these groups
  df_groups     = as.data.frame(as.character(these_clades))
  # summarize the results
  tb_groups = as.data.frame(table(df_groups))
  tb_groups = arrange(tb_groups, desc(Freq))
  tb_groups$Freq = tb_groups$Freq/sample_size
  # save in a list
  summary_sis_groups[[i]] = tb_groups
  
}


# # Make plots
plots = vector("list", length = length(fossil_names))


# max number of sister groups to consider
n <- 50

# for each fossil
for (i in 1:length(fossil_names)){
  
  # get this df
  this_df           = summary_sis_groups[[i]]
  colnames(this_df) = c("clades", "PP")
  
  # if the dataset is large, trimm it 
  
  if (nrow(this_df) > n){
    
    this_df = summary_sis_groups[[i]][1:n, ]
    colnames(this_df) = c("clades", "PP")
    
    # and plot it  
    plot <- ggplot(this_df, aes(x=factor(clades, levels = clades), y = PP)) +
      geom_col(fill = "#440154FF") +
      theme_bw() +
      theme(axis.text.x = element_blank())+
      xlab(paste0(fossil_names[i], " sister groups")  ) +
      scale_y_continuous(limits = c(0, 0.65))
    
    plots[[i]] <- plot
    
  } else {
    
    # make a plot
    plot <- ggplot(this_df, aes(x=factor(clades, levels = clades), y = PP)) +
      geom_col(fill = "#440154FF") +
      theme_bw() +
      theme(axis.text.x = element_blank())+
      xlab(paste0(fossil_names[i], " sister groups") ) +
      scale_y_continuous(limits = c(0, 0.65))
    
    plots[[i]] <- plot
  }
}


# # save figures
# for (f in 1:length(plots)){
#   
#   pdf(paste0("figures/Manuscript_figs/sister_groups_PP_one_fossil_aat/", fossil_names[f], ".pdf"), width=4, height=2.5)
#   
#   plot(plots[[f]])
#   
#   dev.off()
#   
# }




