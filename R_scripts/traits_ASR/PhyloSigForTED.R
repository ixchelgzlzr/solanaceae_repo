library(phytools)
library(TreeTools)
library(ape)
library(expm)
library(caper)

setwd("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/")

##### Phylogenetic Signal for Continuous Traits ################

contdata<-read.csv("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/continuous_extantonly.csv", row.names=1, header=T)
contdata[contdata == "?"]<-NA

trees<-read.tree("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/sample_100_extant_tHet_trees.tree")
tree<-trees[[1]] #this is just for grabbing tip labels later
trait<-as.numeric(contdata$fruit_pedicel) #this is where we select the trait
names(trait)<-row.names(contdata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs
to.drop<-setdiff(tree$tip.label,names(trait))
trees<-lapply(trees,drop.tip,tip=to.drop) #drop tips that had NA

#goal is to run tests for phy signal on all characters across all trees
#and create a summary table with min across trees, max across trees, average, and number of signif
#I think I'll do character by character since there are only 8 and compile separately

Kvals<-c()
Pvals<-c()

for(i in 1:length(trees)) {
	tree<-trees[[i]]
	test<-phylosig(tree,trait,test=T)
	Kvals[i]<-test$K
	Pvals[i]<-test$P
}

results<-cbind(Kvals,Pvals)
write.csv(results,"phylosig_fruit_pedicel.csv") #this needs to match trait name above

min(Kvals)
max(Kvals)
mean(Kvals)
min(Pvals)
max(Pvals)
mean(Pvals) #I pasted this values into a table.

##### Phylogenetic Signal for Binary Traits ################

# will compute FPD for the subset of traits that are binary

catedata<-read.csv("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/morphology.csv", row.names=1, header=T)
catedata[catedata == "?"]<-NA # turn ? into NA
trees<-read.tree("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/sample_100_extant_tHet_trees.tree")

trait<-catedata$calyx_teeth #this is where we select the trait
names(trait)<-row.names(catedata)
trait<-trait[!is.na(trait)] # remove NAs
data<-as.data.frame(trait) # caper wants a df
data[,2]<-row.names(data) # and it wants taxon names in a column called 'binomial' by default
colnames(data)<-c("trait","binomial")
to.drop<-setdiff(tree$tip.label,names(trait))
trees<-lapply(trees,drop.tip,tip=to.drop) #drop tips that had NA

DEstimate<-c()
Pval1s<-c()
Pval0s<-c()

for(i in 1:length(trees)) {
	tree<-trees[[i]]
	input<-comparative.data(tree,data=data, names.col="binomial") # special S3 object
	FPD<-phylo.d(input, binvar=trait, permut = 1000, rnd.bias=NULL)
	DEstimate[i]<-FPD$DEstimate
	Pval1s[i]<-FPD$Pval1
	Pval0s[i]<-FPD$Pval0
}

results<-cbind(DEstimate,Pval1s,Pval0s)
write.csv(results,"phylosig_calyx_teeth.csv") #this needs to match trait name above

min(DEstimate)
max(DEstimate)
mean(DEstimate)
min(Pval1s)
max(Pval1s)
mean(Pval1s) 
min(Pval0s)
max(Pval0s)
mean(Pval0s) #I pasted this values into a table.


##### Phylogenetic Signal for Categorical Traits ################

source('~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/code.R', chdir = TRUE)
catedata<-read.csv("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/morphology.csv", row.names=1, header=T)
catedata[catedata == "?"]<-NA

# for Borges' delta, the tip data need to be in the same order as the tree
# I think that may not be the same across all trees so I will sort every time

#SOME PARAMETERS... 
lambda0 <- 0.1   #rate parameter of the proposal 
se      <- 0.5   #standard deviation of the proposal
sim     <- 1000 #number of iterations
thin    <- 10    #we kept only each 10th iterate 
burn    <- 100   #100 iterates are burned-in

# clunky but I don't feel like prettying -- will just swap out character names as I go
# compression,hilum_position,embryo,cavity,testal_walls,wings,elaiosomes,base_invaginated,angled,
# lobe_sinus,widest_veins,secondary_veins,secondary_veins_fork,fruit_type,calyx_inflated,venation_pattern,calyx_teeth

# the code below generates a list of deltaA values and p-values across 100 trees for one character
# and it takes forever

trees<-read.tree("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/sample_100_extant_tHet_trees.tree")
tree<-trees[[1]] #this is just for grabbing tip labels later
trait<-catedata$testal_walls #this is where we select the trait
names(trait)<-row.names(catedata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs
to.drop<-setdiff(tree$tip.label,names(trait))
trees<-lapply(trees,drop.tip,tip=to.drop) #drop tips that had NA

deltaA<-c()
for(i in 1:100) {
	tree<-trees[[i]]
	match(tree$tip.label, names(trait))->Sort
	trait<-trait[order(match(names(trait),tree$tip.label))] #resorting every time to match tree
	deltaA[i]<-delta(trait,tree,lambda0,se,sim,thin,burn)
}

pvals<-c()
for(i in 1:100) {
	tree<-trees[[i]]
	random_delta<-rep(NA,100)
	for(j in 1:100){
		deltaA[j]->origdelta
		rtrait<-sample(trait) #we don't need to sort since we are shuffling anyway
		random_delta[j]<-delta(rtrait,tree,lambda0,se,sim,thin,burn)
		}
	pvals[i]<-sum(random_delta>origdelta)/length(random_delta)
}

deltaresults<-cbind(deltaA,pvals)
write.csv(deltaresults,"delta_cavity.csv")

min(deltaA)
max(deltaA)
mean(deltaA)
min(pvals)
max(pvals)
mean(pvals)

# results with the above are *weird* -- tons of variation across trees, check by replicating just one tree

#SOME PARAMETERS... 
lambda0 <- 0.1   #rate parameter of the proposal 
se      <- 0.5   #standard deviation of the proposal
sim     <- 10000 #number of iterations
thin    <- 10    #we kept only each 10th iterate 
burn    <- 100   #100 iterates are burned-in

trees<-read.tree("~/Library/CloudStorage/OneDrive-UCB-O365/Documents/Papers/TED/PhylSignal/sample_100_extant_tHet_trees.tree")
tree<-trees[[4]]
trait<-catedata$cavity #this is where we select the trait
names(trait)<-row.names(catedata) #attach names to data
trait<-trait[!is.na(trait)] #remove NAs
to.drop<-setdiff(tree$tip.label,names(trait))
tree<-drop.tip(tree, to.drop)
deltaA<-delta(trait,tree,lambda0,se,sim,thin,burn)
deltaA
random_delta<-rep(NA,100)
	for(j in 1:100){
		rtrait<-sample(trait) #we don't need to sort since we are shuffling anyway
		random_delta[j]<-delta(rtrait,tree,lambda0,se,sim,thin,burn)
		}
pval<-sum(random_delta>deltaA)/length(random_delta)
pval

# conclusion: the wide variation across trees is real BUT ALSO there is a lot of variation 
# across repeated searches of the same tree, enough that sometimes it would seem like delta
# is signficant and other times not.  This instability must mean that the proposal or mcmc settings 
# are bad. Would need a good deal of troubleshooting. More simulutions was not enough to fix. Boo.









