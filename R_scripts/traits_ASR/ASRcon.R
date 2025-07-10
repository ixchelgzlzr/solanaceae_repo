# Set working directory

library(phytools)
library(ape)
library(geiger)
library(picante)
library(ggtree)
library(ggplot2)


contdata<-read.csv("continuous_extantonly.csv", row.names=1, header=T)
contdata[contdata == "?"]<-NA


# Load tree file must be in newick format
read.newick ("MCC_extant_tree.nwk") -> Tree
class(Tree)


## 1.Seed length
trait<-as.numeric(contdata$seed_L) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop

# Plot tree
plot(Tree, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines

# Reconstruct node data - not working
node_df <- data.frame(node = as.numeric(names(fit$ace)),
                      state = fit$ace)

# Plot directly - not working
ggplot(Tree) +
  geom_tree(aes(x = x, y = y, color = edge_state), size = 1) +
  geom_tiplab(aes(label = label), size = 2.5) +
  scale_color_viridis_c(option = "plasma", name = "Seed length") +
  theme_tree2() +
  theme(legend.position = "right")



## 2.Seed ratio
trait<-as.numeric(contdata$seed_ratio) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop

# Plot tree
plot(Tree, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines



## 3. Length fruiting calyx
trait<-as.numeric(contdata$frt_calyx_L) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop

# Plot tree
plot(Tree, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines



## 4. Fruiting calyx ratio
trait<-as.numeric(contdata$frt_calyx_ratio) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop

# Plot tree
plot(Tree, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines



## 5. Length calyx lobes
trait<-as.numeric(contdata$calyx_lobe_L) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop

# Plot tree
plot(Tree, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines



## 6. Calyx lobes: total calyx length ratio
trait<-as.numeric(contdata$calyx_lobe_ratio) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop
Tree_cut<-drop.tip(Tree,to.drop) #drop tips that had NA

difference<-name.check(Tree_cut,trait)
difference

# Plot tree
plot(Tree_cut, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree_cut,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree_cut,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree_cut)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines



## 7. Fruit pedicel
trait<-as.numeric(contdata$fruit_pedicel) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop

# Plot tree
plot(Tree, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines



## 7. Fruit ratio
trait<-as.numeric(contdata$fruit_ratio) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs

to.drop<-setdiff(Tree$tip.label,names(trait))
to.drop

# Plot tree
plot(Tree, cex = 0.6)

attach(contdata)

## ASRcon
fit<-fastAnc(Tree,trait,vars=TRUE,CI=TRUE)
fit #point estimates, variance, and CIs for nodes

mapped<-contMap(Tree,trait,plot=FALSE) #these reconstructions are often displayed as heatmaps

plot(mapped,
     legend = 0.7 * max(nodeHeights(Tree)),  # where to place the legend
     fsize = c(0.4, 0.9),                    # font sizes (tip, internal)
     lwd = 2,                                # thicker branches
     outline = FALSE)                        # remove dark outlines
