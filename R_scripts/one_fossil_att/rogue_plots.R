install.packages("Rogue")

library(ggpubr)

tree.trace  <- read.tree("output/final_runs/time_heterogeneous/trees_combined.trees")
sample_size <- 13900
trees       <- tree.trace[1:sample_size]


fossil_names <- c(  "f_Physalis_hunickenii",
                    "f_Physalis_infinemundi",
                    "H4895_Hyoscyamus_undulatus",
                    "2003_0073_Nicotianosemen_minuta",
                  "CCN6065_second_seed_Capsicum_pliocenicum",
                  "f_Eophysaloides_inflata",
                  "f_Lycianthoides_calycina",
                  "K453A_Hyoscyamus_datureoides",
                  "K528_Solanum_miocenicus",
                  "K530_Solanoides_dorofeevii",
                  "K587B_Solanum_foveolatum",
                  "UF6500_Nephrosemen_reticulatus",
                  "V_40891_Solanispermum_reniforme",
                  "V_40898_Solanum_arnense")


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




correct_names <- read.csv(file = "data/conversion_table.csv", header = T)

correct_names <- filter(correct_names, correct_names$old_name %in% fossil_names)

correct_names <- correct_names[match(fossil_names, correct_names$old_name), ]


# make rogue plots


for (i in 4:length(fossil_names)){
  
  pdf(paste0("figures/Manuscript_figs/rogue_plots/", fossil_names[i], ".pdf"), width = 7, height = 9 )
  
  plotted <- RoguePlot(trimmed_trace[[i]], fossil_names[i], legend = "none", legend.inset = 0.02, main = paste0(correct_names$correct_name[i]), cex= 0.5)
  
  PlotTools::SpectrumLegend(
    "bottomleft",
    palette = colorRampPalette(c(par("fg"), "#009E73"), space = "Lab")(100),
    legend = plotted$legendLabels)
  
  
  dev.off()
  
}




