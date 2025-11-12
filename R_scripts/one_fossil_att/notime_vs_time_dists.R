########################################################################
# FOSSIL POSITIONS ONE FOSSIL AT THE TIME
# For the no time trees
#
# NOTE TO SELF: this can be done including all fossils simultaneously, 
# but for consistency we will report these resulst in the manuscript. 
# Including fossils can be run with the script "fossil_position.R"
# Inferences do not change!
#########################################################################
library(ape)
library(dplyr)
library(ggplot2)
library(phangorn)

source("R_scripts/one_fossil_att/fossils_positions_one_fossil_aat.R")


# # all the fossils
# fossil_names <- c("2003_0073_Nicotianosemen_minuta", 
#                   "CCN6065_second_seed_Capsicum_pliocenicum",
#                   "f_Eophysaloides_inflata",
#                   "f_Lycianthoides_calycina",
#                   "f_Physalis_hunickenii",
#                   "f_Physalis_infinemundi",
#                   "H4895_Hyoscyamus_undulatus",
#                   "K453A_Hyoscyamus_datureoides",
#                   "K528_Solanum_miocenicus",
#                   "K530_Solanoides_dorofeevii",
#                   "K587B_Solanum_foveolatum",
#                   "UF6500_Nephrosemen_reticulatus",
#                   "V_40891_Solanispermum_reniforme",
#                   "V_40898_Solanum_arnense")

notime.tree.trace  <- read.tree("output_no_time/_moleQ_GTR_G_morphQ_Mk_G_MCMC_uniform_1/tree.trees")

# discard the first half of the trace for burnin
notime.trees <- tree.trace[3000:length(notime.tree.trace)]

notime.sample_size = length(notime.trees)

#------------------------------------------
# Make and MCC tree
#------------------------------------------
# 
# MCC <- maxCladeCred(tree.trace)
# write.tree(MCC, "output_no_time/_moleQ_GTR_G_morphQ_Mk_G_MCMC_uniform_1/R_MCC_tree.tre")
# 
# 
# 
# plot.phylo(MCC, cex = 0.6, label.offset = 0.01)
# 
# rooted_MCC <- root(MCC, node = getMRCA(MCC, c("Schizanthus_grahamii", "Duckeodendron_cestroides")) )
# 
# plot(rooted_MCC, cex = 0.6, label.offset = 0.01)
# 
# 
# pdf("output_no_time/_moleQ_GTR_G_morphQ_Mk_G_MCMC_uniform_1/R_MCC.pdf", width = 7, height = 15)
# plot(rooted_MCC, cex = 0.6, label.offset = 0.01)
# dev.off()

#------------------------------------------
# Extract fossils and its sister clades
#------------------------------------------

# make empty container
notime.trimmed_trace  = vector(mode = "list", length = length(fossil_names))

# for each fossil
for (i in 1:length(fossil_names)){
  
  # get this fossil name
  this_fossil = fossil_names[[i]]
  
  # Get the vector of fossils to drop from sample
  names_to_drop <- setdiff(fossil_names, this_fossil)
  
  # drop tips from the trace 
  for (t in 1:length(notime.trees)){
    
    notime.trimmed_trace[[i]][[t]] <- drop.tip(notime.trees[[t]], names_to_drop)
    
  }
}


# empty list
notime_fossil_and_sister = vector(mode = "list", length = length(fossil_names))

# for each fossil
for (f in 1:length(notime.trimmed_trace)){
  
  # get this fossil name
  this_fossil = fossil_names[[f]]
  
  # for each tree of the trace
  for (t in 1:length(notime.trimmed_trace[[f]])){
    
    # get the tip number for the fossil
    tip_number = which(notime.trimmed_trace[[f]][[t]]$tip.label == this_fossil)
    
    # get the node number for that tip by using the edges
    this_node_row    = which(notime.trimmed_trace[[f]][[t]]$edge[ , 2] == tip_number)
    this_node_number = notime.trimmed_trace[[f]][[t]]$edge[this_node_row, 1] 
    
    # extract the clade originating from that node
    this_clade = extract.clade(notime.trimmed_trace[[f]][[t]], node = this_node_number)
    
    # save fossil + sister in list
    notime_fossil_and_sister[[f]][[t]] = this_clade
    
  } }




#------------------------------------------
# Summarize
#------------------------------------------


# get the tip names per each tree

# empty contained
notime_tips_fossil_and_sister = vector(mode = "list", length = length(notime_fossil_and_sister))

# get tips

# for each fossil
for (f in 1:length(notime_fossil_and_sister)){
  
  # for each tree
  for(t in 1:length(notime.trees)){
    
    # get the tip labels
    notime_tips_fossil_and_sister[[f]][[t]] = notime_fossil_and_sister[[f]][[t]]$tip.label
    
  }
}



# get the unique clades to which each fossil is sister to
notime_summary_sis_groups = vector("list", length = length(notime_tips_fossil_and_sister))


# for each fossil
for (i in 1:length(notime_tips_fossil_and_sister)){
  
  # unlist the tips names
  these_clades  = lapply(notime_tips_fossil_and_sister[[i]], unlist)
  these_clades  = lapply(these_clades, sort)
  # extract a df with the sister clades for these groups
  df_groups     = as.data.frame(as.character(these_clades))
  # summarize the results
  tb_groups = as.data.frame(table(df_groups))
  tb_groups = arrange(tb_groups, desc(Freq))
  tb_groups$Freq = tb_groups$Freq/notime.sample_size
  # save in a list
  notime_summary_sis_groups[[i]] = tb_groups
  
}


# # Make plots
notime_plots = vector("list", length = length(fossil_names))


# max number of sister groups to consider
n <- 50

# for each fossil
for (i in 1:length(fossil_names)){
  
  # get this df
  this_notime_df    = notime_summary_sis_groups[[i]]
  this_df           = summary_sis_groups[[i]]
  colnames(this_notime_df) = c("clades", "PP")
  colnames(this_df) = c("clades", "PP")
  
  
  # check and trim if necessary
  if (nrow(this_notime_df) > n | nrow(this_df) > n) {
    if (nrow(this_notime_df) > n) {
      this_notime_df <- this_notime_df[1:n, ]
    }
    if (nrow(this_df) > n) {
      this_df <- this_df[1:n, ]
    }
  }
  
  # let's rename the
  this_df$clades <- c(1:nrow(this_df))
  this_notime_df$clades <- c(1:nrow(this_notime_df))


  
  # Combine data first
  combined_df <- rbind(
    cbind(this_df, type = "time"),
    cbind(this_notime_df, type = "no time")
  )
  
  # Plot
  plot <- ggplot(combined_df, aes(x = clades, y = PP, fill = type)) +
    geom_col(alpha = 0.7, position = "identity") +
    scale_fill_manual(
      name = "Analysis type",
      values = c("time" = "#440154FF", "no time" = "#26828EFF")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      legend.position = "right"
    ) +
    xlab(paste0(fossil_names[i], " sister groups"))
  
  
  
  notime_plots[[i]] <- plot
  
  
  # # let's rename the 
  # this_df$clades <- c(1:nrow(this_df))
  # this_notime_df$clades <- c(1:nrow(this_notime_df))
  # 
  # 
  #   # then make a plot
  # ggplot() +
  #   geom_col(data = this_df, aes(x=factor(clades, levels = clades), y = PP), fill = alpha("#440154FF", 0.5)) +
  #   geom_col(data = this_notime_df, aes(x=factor(clades, levels = clades), y = PP), fill = alpha("#26828EFF", 0.5)) +
  #   theme_bw() +
  #   theme(axis.text.x = element_blank(), legend.position = "top") +
  #   xlab(paste0(fossil_names[i], " sister groups"))   +
  #   scale_fill_manual(
  #     name = "Analysis type",
  #     values = c("time" = "#440154FF", "no_time" = "#26828EFF")
  #   ) 
    
    notime_plots[[i]] <- plot
    
    
}


library(ggpubr)

ggarrange(plotlist= unlist (notime_plots), 
          ncol = 4, nrow = 3, common.legend = T)


#export figure
pdf("figures/Manuscript_figs/submission_figs/Sup_mat_figs/Fig_S12.pdf", width = 12, height = 12)

ggarrange(plotlist= unlist (notime_plots), 
          ncol = 3, nrow = 4, common.legend = T)

dev.off()




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




