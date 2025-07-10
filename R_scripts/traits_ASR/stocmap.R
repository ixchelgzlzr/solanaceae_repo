
# Set working directory

# Run this every time
# 	 install.packagees("phytools")
# 	 when running the first time
require(phytools)
library(ape)
library(geiger)
library(phangorn)
library(picante)
library(phylometrics)

# Load tree file must be in newick format
read.tree("MCC_extant_tree.nwk") -> Tree
# Plot tree
plot(Tree, cex = 0.5)


# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

Tree$tip.label # its tip labels

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ] # sort data along tips order of first tree
Data_SORT<-as.data.frame(Data_SORT) # make data frame
row.names(Data_SORT)<-Tree$tip.label # give it row names


## 1. Compression
names(Data_SORT)[1]<-"compression"
head(Data_SORT)

# full dataset
vec<-Data_SORT$compression
names(vec)<-row.names(Data_SORT)

to.drop<-setdiff(Tree$tip.label,names(vec))
to.drop

Tree_cut<-drop.tip(Tree,to.drop) #drop tips that had NA
difference<-name.check(Tree_cut,vec)
difference

###### Stochastic Mapping
cols <- c("Black" , "Red")
names(cols)<-c("flatenned","not flattened")

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = vec, model ="ARD",nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("compression.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 2. Hilum position
names(Data_SORT)[2]<-"hilum_position"
head(Data_SORT)

# full dataset
vec<-Data_SORT$hilum_position
names(vec)<-row.names(Data_SORT)

to.drop<-setdiff(Tree$tip.label,names(vec))
to.drop

Tree_cut<-drop.tip(Tree,to.drop) #drop tips that had NA
difference<-name.check(Tree_cut,vec)
difference

###### Stochastic Mapping
cols <- c("Black" , "Red")
names(cols)<-c("lateral or sub-lateral","terminal")

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = vec, model ="ARD",nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("hilum position.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 3. Embryo shape
# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[3]<-"embryo"
head(Data_SORT)

# full dataset
vec<-Data_SORT$embryo
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:3]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.33 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
    if(vec_NO3[i]==2){
      mat1[i,3]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("embryo shape.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 4. Hilar-chalazal cavity
# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[4]<-"cavity"
head(Data_SORT)

# full dataset
vec<-Data_SORT$cavity
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Hilar-chalazal cavity.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 5. Exotestal cell walls
# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[5]<-"testal_walls"
head(Data_SORT)

# full dataset
vec<-Data_SORT$testal_walls
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("straight","sinuate-cerebelloid")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Exotestal cell walls.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 6. Seed wings
# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[6]<-"wings"
head(Data_SORT)

# full dataset
vec<-Data_SORT$wings
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Seed wings.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 7. Seed elaiosomes
# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[7]<-"elaiosome"
head(Data_SORT)

# full dataset
vec<-Data_SORT$elaiosome
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Seed elaiosomes.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 8. Fruiting calyx base invagination

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[8]<-"base_invaginated"
head(Data_SORT)

# full dataset
vec<-Data_SORT$base_invaginated
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx base invagination.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 9. Fruiting calyx angled

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[9]<-"angled"
head(Data_SORT)

# full dataset
vec<-Data_SORT$angled
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx angled.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 10. Fruiting calyx lobe sinus

# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[10]<-"lobe_sinus"
head(Data_SORT)

# full dataset
vec<-Data_SORT$lobe_sinus
names(vec)<-row.names(Data_SORT)
vec

class(vec)
class(Tree)

cols <- c("Black" , "Red", "Blue")
names(cols)<-c("rounded","angular", "absent")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:3]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.33 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
    if(vec_NO3[i]==2){
      mat1[i,3]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx lobe sinus.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 11. Fruiting calyx widest veins

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[11]<-"widest_veins"
head(Data_SORT)

# full dataset
vec<-Data_SORT$widest_veins
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx widest veins.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 12. Fruiting calyx secondary veins

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[12]<-"secondary_veins"
head(Data_SORT)

# full dataset
vec<-Data_SORT$secondary_veins
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx secondary veins.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 13. Fruiting calyx secondary forking veins

# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[13]<-"secondary_veins_fork"
head(Data_SORT)

# full dataset
vec<-Data_SORT$secondary_veins_fork
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx secondary forking veins.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 15. Fruiting calyx inflated

# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[15]<-"calyx_inflated"
head(Data_SORT)

# full dataset
vec<-Data_SORT$calyx_inflated
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx inflated.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 16. Fruiting calyx venation pattern

# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[16]<-"venation_pattern"
head(Data_SORT)

# full dataset
vec<-Data_SORT$venation_pattern
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("parallelodromous","other type")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx venation pattern.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device



## 17. Fruiting calyx teeth

# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[17]<-"calyx_teeth"
head(Data_SORT)

# full dataset
vec<-Data_SORT$calyx_teeth
names(vec)<-row.names(Data_SORT)

class(vec)
class(Tree)

cols <- c("Black" , "Red")
names(cols)<-c("absent","present")

vec_NO3<-vec # save to new vector
mat1<-to.matrix(vec_NO3, letters[1:2]) # make a matrix with 2 columns to store probability of each state
for(i in 1:length(vec_NO3)){
  if(is.na(vec_NO3[i])){
    mat1[i,]<-0.5 # if there is NA, put 0.5 in each column
  }else{
    if(vec_NO3[i]==0){
      mat1[i,1]<-1 # if in state 0, put a 1 in first column
    }
    if(vec_NO3[i]==1){
      mat1[i,2]<-1 # if in state 2 put a 1 in second column
    }
  }
}

class(mat1)

FCSimmapTree_ARD1000 <- make.simmap(tree = Tree, x = mat1, model ="ARD", nsim=1000, Q = "empirical")
sum.SIMMAP_ARD<-describe.simmap(FCSimmapTree_ARD1000)

counts<-countSimmap(FCSimmapTree_ARD1000,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,1]) #for 0 to 1 transitions, for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,1]) # gives 95% credibility interval

pdf("Fruiting calyx teeth.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sum.SIMMAP_ARD$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=7, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sum.SIMMAP_ARD$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device


## 14. Fruit type

# Load trait data
Data<-read.csv("morphology.csv", row.names='species')
head(Data)

# Sort taxa orders
Data_SORT <- Data[Tree$tip.label, ]
Data_SORT<-as.data.frame(Data_SORT)
row.names(Data_SORT)<-Tree$tip.label
names(Data_SORT)[14]<-"fruit_type"
head(Data_SORT)

# full dataset
vec<-Data_SORT$fruit_type
names(vec)<-row.names(Data_SORT)

#stochasting mapping
simmap<-make.simmap(tree = Tree, x = vec, model = "ARD", nsim = 1000, Q = "empirical") # this takes a while 
sumobj<-describe.simmap(simmap) # summarize shifts
sumobj

#to summarize median and 95% credibility interval of number of transitions between states
counts<-countSimmap(simmap,message=FALSE) 
head(counts) #you can see that each column is one transition type 
median(counts[,9]) #for example, generally preferred over means for non-normal distributions
diversitree:::hdr->hdr #to pull hdr function from diversitree 
hdr(counts[,9]) # gives 95% credibility interval

cols <- c("Black" , "Red", "Blue", "Green", "Purple", "Orange")
names(cols)<-c("berry","capsule", "diclesium or mericarp", "drupe", "pyxidium", "non-capsular dehiscent fruits")

pdf("Fruit type.pdf", width = 10, height = 10, useDingbats = F) # make pdf for it
plotTree(Tree, fsize=0.5, ftype="i", lwd = 2, plot = T, offset = 0.5)  # plot tree (this is the 2nd map from the posterior sample)
nodelabels(pie=sumobj$ace,piecol=cols,cex=0.30) # add pies at nodes 
add.simmap.legend(x=0.1*max(nodeHeights(Tree)),y=25, colors=cols,prompt=F, y=5, cex = 2) # add legend
tiplabels(pie = sumobj$tips, piecol = cols, cex = 0.1, adj = 0.5) # add tip pies with probabilities from posterior
dev.off() # turn off pdf device

