#############################
# Read and merge the traces #
#############################

library(ape)

# burnin per trace, gotta define manually
#burnin <- c(500, 500) #homogeneous
burnin <- c(500, 1600) #heterogeneous

# set up the directory for the type of analysis
  # this assumes that each directory has subdirectories with the same 
  # type of analyses, and with files named the same way, 
  # e.g. trees.tree in each folder
#dir = "output/final_runs/time_homogeneous"
dir      = "output/final_runs/time_heterogeneous"
sub.dirs = list.dirs(dir, recursive = F)

# Get names of traces to read
trace.files             = vector("list", length(sub.dirs))
tree.trace.files        = vector("list", length(sub.dirs))
extant.tree.trace.files = vector("list", length(sub.dirs))

# create file names
for (i in 1:length(sub.dirs)){
  trace.files[[i]]             = paste0( sub.dirs[i], "/params.log" )
  tree.trace.files[[i]]        = paste0( sub.dirs[i], "/tree.trees")
  extant.tree.trace.files[[i]] = paste0( sub.dirs[i], "/extant_tree.trees")
}


  # empty vectors
traces       = vector("list", length(trace.files))
rows_to_drop = vector("list", length(trace.files))

##################
# MERGE THE LOGS #
##################

# for each sample
for (i in 1:length(trace.files)){
  
  # read in the trace
  traces[[i]] = read.table(trace.files[[i]], header = T)
  
  # decide which samples to drop
  rows_to_drop[[i]] = which(apply(traces[[i]], 1, function(x) any(is.nan(x)) ))
  rows_to_drop[[i]] = sort(c(1:burnin[i], rows_to_drop[[i]]))
  
  # drop NA and burnin
  traces[[i]] = traces[[i]][-rows_to_drop[[i]],]
  
}  


# combine traces in one
combined_traces = do.call(rbind, traces)

# re-number iterations
combined_traces$Iteration = 1:nrow(combined_traces)

# write combined traces to file
write.table(combined_traces, file=paste0(dir,"/params_combined.log"), row.names=FALSE, sep="\t", quote=FALSE)


########################
# MERGE THE FULL TREES #
########################

# empty containers
tree.traces    = vector("list", length(tree.trace.files))
trees_to_drop  = vector("list", length(tree.trace.files))

# for each sample
for (i in 1:length(tree.trace.files)){
  
  # read in the tree traces
  tree.traces[[i]] = read.table(tree.trace.files[[i]], header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

  
  # decide which samples to drop
  trees_to_drop[[i]] = which(apply(tree.traces[[i]], 1, function(x) any(is.nan(x)) ))
  trees_to_drop[[i]] = sort(c(1:burnin[i], trees_to_drop[[i]]))
  
  # drop NA and burnin
  tree.traces[[i]] = tree.traces[[i]][-trees_to_drop[[i]],]
  
}  


# combine traces in one
combined_tree_traces = do.call(rbind, tree.traces)

# re-number iterations
combined_tree_traces$Iteration = 1:nrow(combined_tree_traces)

# write combined traces to file
write.table(combined_tree_traces, file=paste0(dir,"/trees_combined.trees"), row.names=FALSE, sep="\t", quote=FALSE)



##########################
# MERGE THE EXTANT TREES #
##########################


# empty containers
extant.tree.traces    = vector("list", length(extant.tree.trace.files))
extant_trees_to_drop  = vector("list", length(extant.tree.trace.files))

# for each sample
for (i in 1:length(extant.tree.trace.files)){
  
  # read in the tree traces
  extant.tree.traces[[i]] = read.table(extant.tree.trace.files[[i]], header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  
  
  # decide which samples to drop
  extant_trees_to_drop[[i]] = which(apply(extant.tree.traces[[i]], 1, function(x) any(is.nan(x)) ))
  extant_trees_to_drop[[i]] = sort(c(1:burnin[i], extant_trees_to_drop[[i]]))
  
  # drop NA and burnin
  extant.tree.traces[[i]] = extant.tree.traces[[i]][-extant_trees_to_drop[[i]],]
  
}  


# combine traces in one
combined_extant_tree_traces = do.call(rbind, extant.tree.traces)

# re-number iterations
combined_extant_tree_traces$Iteration = 1:nrow(combined_extant_tree_traces)

# write combined traces to file
write.table(combined_extant_tree_traces, file=paste0(dir,"/extant_trees_combined.trees"), row.names=FALSE, sep="\t", quote=FALSE)

