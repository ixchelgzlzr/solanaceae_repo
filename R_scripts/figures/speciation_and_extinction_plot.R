###############################
# Plots of Speciation and extinction under time heterogeneous model
##############################

# libraries
library(RevGadgets)
library(deeptime)


# specify the input file
file <- "output/final_runs/time_heterogeneous/params_combined.log"

# read the trace and discard burnin
trace_quant <- readTrace(path = file, burnin = 0.1)

# read time intervals
time_intervals = read.csv("data/epoch_timescale.csv")
time_intervals = time_intervals[1:8, ]

# extract the lambdas and mus
lambdas = paste0("lambda.", c(1:8), ".")
mus     = paste0("mu.", c(1:8), ".")
divs    = paste0("diversification.", c(1:8), ".")
psis    = paste0("psi.", c(1:8), ".")

# Summarize the lambdas and mus independently
summary_lambdas = summarizeTrace(trace = trace_quant, vars =  lambdas)
summary_mus     = summarizeTrace(trace = trace_quant, vars =  mus)
summary_divs    = summarizeTrace(trace = trace_quant, vars =  divs)
summary_psis    = summarizeTrace(trace = trace_quant, vars =  psis)



# make empty vectors
mean_l          = vector(length=length(time_intervals))
quantile_2.5_l  = vector(length=length(time_intervals))
quantile_97.5_l = vector(length=length(time_intervals))
mean_m          = vector(length=length(time_intervals))
quantile_2.5_m  = vector(length=length(time_intervals))
quantile_97.5_m = vector(length=length(time_intervals))
mean_d          = vector(length=length(time_intervals))
quantile_2.5_d  = vector(length=length(time_intervals))
quantile_97.5_d = vector(length=length(time_intervals))
mean_p          = vector(length=length(time_intervals))
quantile_2.5_p  = vector(length=length(time_intervals))
quantile_97.5_p = vector(length=length(time_intervals))



# get the summaries in a for loop
for (i in 1:nrow(time_intervals)){
  
  # get lambdas summaries
  mean_l[i]          = summary_lambdas[[i]][[1]]["mean"]
  quantile_2.5_l[i]  = summary_lambdas[[i]][[1]]["quantile_2.5"]
  quantile_97.5_l[i] = summary_lambdas[[i]][[1]]["quantile_97.5"]
  
  # get mu summaries
  mean_m[i]          = summary_mus[[i]][[1]]["mean"]
  quantile_2.5_m[i]  = summary_mus[[i]][[1]]["quantile_2.5"]
  quantile_97.5_m[i] = summary_mus[[i]][[1]]["quantile_97.5"]
  
  # get divs summaries
  mean_d[i]          = summary_divs[[i]][[1]]["mean"]
  quantile_2.5_d[i]  = summary_divs[[i]][[1]]["quantile_2.5"]
  quantile_97.5_d[i] = summary_divs[[i]][[1]]["quantile_97.5"]
  
  # get psis summaries
  mean_p[i]          = summary_psis[[i]][[1]]["mean"]
  quantile_2.5_p[i]  = summary_psis[[i]][[1]]["quantile_2.5"]
  quantile_97.5_p[i] = summary_psis[[i]][[1]]["quantile_97.5"]
  
}

# make data frames
df_lambdas = data.frame(time_intervals$start, time_intervals$end, mean_l, quantile_2.5_l, quantile_97.5_l)
df_mus     = data.frame(time_intervals$start, time_intervals$end, mean_m, quantile_2.5_m, quantile_97.5_m)
df_divs    = data.frame(time_intervals$start, time_intervals$end, mean_d, quantile_2.5_d, quantile_97.5_d)
df_psis    = data.frame(time_intervals$start, time_intervals$end, mean_p, quantile_2.5_p, quantile_97.5_p)


# plot lambdas
spc_plot <- ggplot(df_lambdas) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = mean_l), data = df_lambdas[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_2.5_l), linetype = "dashed", data = df_lambdas[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_97.5_l), linetype = "dashed", data = df_lambdas[1:7, ]) +
 geom_stepribbon(aes(x=time_intervals.end,
                 ymin=quantile_2.5_l, 
                 ymax=quantile_97.5_l), alpha=0.5,
                 fill = "#31688EFF") + 
  ggplot2::scale_x_reverse() +
  theme_bw() +
  theme(axis.title.x=element_blank()) +
  ylab("Speciation")  +
    coord_geo(dat = list("periods", "epochs"), height = unit(1, "line")) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=9))

spc_plot

# plot mus
extplot <- ggplot(df_mus) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = mean_m), data = df_mus[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_2.5_m), linetype = "dashed", data = df_mus[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_97.5_m), linetype = "dashed", data = df_mus[1:7, ]) +
  geom_stepribbon(aes(x=time_intervals.end,
                      ymin=quantile_2.5_m, 
                      ymax=quantile_97.5_m), alpha=0.5,
                  fill = "#440154FF") + 
  ggplot2::scale_x_reverse() +
  theme_bw() +
  xlab("Geological time (Mya)") +
  ylab("Extinction") + coord_geo(dat = list("periods", "epochs"), height = unit(1, "line")) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=9))

extplot


# plot divs
divplot <- ggplot(df_divs) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = mean_d), data = df_divs[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_2.5_d), linetype = "dashed", data = df_divs[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_97.5_d), linetype = "dashed", data = df_divs[1:7, ]) +
  geom_stepribbon(aes(x=time_intervals.end,
                      ymin=quantile_2.5_d, 
                      ymax=quantile_97.5_d), alpha=0.5,
                  fill = "#35B779FF") + 
  ggplot2::scale_x_reverse() +
  theme_bw() +
  xlab("Geological time (Mya)") +
  ylab("Diversification") + coord_geo(dat = list("periods", "epochs"), height = unit(1, "line")) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=9))

divplot



# plot psis
psiplot <- ggplot(df_psis) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = mean_p), data = df_psis[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_2.5_p), linetype = "dashed", data = df_psis[1:7, ]) +
  geom_segment(aes(x=time_intervals.end, xend = time_intervals.start, y = quantile_97.5_p), linetype = "dashed", data = df_psis[1:7, ]) +
  geom_stepribbon(aes(x=time_intervals.end,
                      ymin=quantile_2.5_p, 
                      ymax=quantile_97.5_p), alpha=0.5,
                  fill = "orchid3") + 
  ggplot2::scale_x_reverse() +
  theme_bw() +
  xlab("Geological time (Mya)") +
  ylab("Fossilization") + 
  coord_geo(dat = list("periods", "epochs"), height = unit(1, "line")) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=9))

psiplot

  

# make a single plot
ggarrange2(nrow = 4, divplot, spc_plot, extplot, psiplot, labels = c("A", "B", "C", "D"))


# # export plot
# tiff("figures/Manuscript_figs/FBD_rates.tiff", 
#      width = 16,
#      height = 18,
#      res = 300,
#      units = "cm")
# 
# ggarrange(nrow = 4, divplot, spc_plot, extplot, psiplot, 
#           labels = c("A", "B", "C", "D"),
#           font.label = list(size = 11, color = "black", face = "bold", family = NULL))
# 
# 
# dev.off()


