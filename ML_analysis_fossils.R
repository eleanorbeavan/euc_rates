# a script for:
    # randomly pruning the branches in a phylogeny to >5 substitutions per branch (200 random repeats)
    # estimating rates using chronos
    # merging the trait data with rate data for each species
    # running PGLS on all 200 trees and extracting stats results

# input: a tree, trait dataframe, translation dataframe (phylogeny species to accepted names)
# output: 5 dataframes (one for each model) with p-values and R^2 for all 200 trees

library(ape)
library(phytools)
library(phylobase)
library(parallel)
library(caper)
library(geiger)
library(nlme)


# read in data
tree = read.tree("~/cp_tree.phy")
seqdata = read.csv("~/Dropbox/euc_sr/CSV_files/alignment_names.csv")
lhtdata = read.csv("~/Dropbox/euc_sr/CSV_files/merged_LHS.csv")
# keep species with complete height data only
height = lhtdata[complete.cases(lhtdata$mean_height),]

# change the tip.labels in the tree so they match accepted_names in the traits dataframe
tree$tip.label=seqdata$accepted_name[match(tree$tip.label, seqdata$alignment_name)]
tree$tip.label = as.character(tree$tip.label)

# three functions to 
# 1. find shortest branch length
# 2. remove shortest branch 
# 3. drop all branches <5/1767

# function to find the shortest branch length in the tree
shortest.tip = function(tree) {
    n.tips = length(tree$tip.label)
    descendant.nodes = tree$edge[,2]
    terminal.edges = tree$edge.length[descendant.nodes<=n.tips]
    return(min(terminal.edges))
}

## write function to loop through the following function 200 times to produce file with 200 pruned trees

# function to remove shortest bl
prune.shortest.branch = function(tree) {
    # take the edge length of those <= the number of tip labels (terminal nodes)
    # vector with terminal branch lengths
    n.tips = length(tree$tip.label)
    descendant.nodes = tree$edge[,2]
    terminal.edges = tree$edge.length[descendant.nodes<=n.tips]
    # round terminal.edges to 5 decimal places
    terminal.edges = round(terminal.edges, digits = 5)
    # keep random shortest tips
    shortest.tips = which(terminal.edges == min(terminal.edges))
    # the loop is needed because if there is only 1 value for shortest tip, the function will return any number from 1:x
    if(length(shortest.tips)>1){ 
        shortest.tip = sample(shortest.tips, 1)
    }else{
        shortest.tip = shortest.tips
    }
    
    new.tree = drop.tip(tree, shortest.tip) 
    return(new.tree)
}

# use a loop to remove all branches <5/1767
x = (5/1767)
prune.tree = function(tree, x) {
    while (shortest.tip(tree) < x) {
        tree = prune.shortest.branch(tree)
    }
    return(tree)
}

# loop to perform "prune.tree" 200 times and create list of 200 pruned trees
pruned.trees = list()
times = 200
for (i in 1:times) {
    newtree = prune.tree(tree, x)
    pruned.trees = c(pruned.trees, list(newtree))
}

#------------------------------------------------------------
# FOSSILS
#------------------------------------------------------------
# create dataframe with node ages
# root age = 52-85
# corymbia/angophora crown = 52+ 
# eucalyptus crown = 37+
tree1 = pruned.trees[[1]]
tree1$tip.label = as.character(tree1$tip.label)
fossils = makeChronosCalib(tree1, node="root", age.min = 52, age.max = 85, interactive = T)

## for loop for fossil data
chronograms = list()
for (i in 1:200) {
    c = chronos(pruned.trees[[i]], lambda = 1, calibration = fossils)
    print(c)
    chronograms = c(chronograms, list(c))
}

# extract rates from 200 trees and create dataframe
extract.rates = function(c) {
    # create dataframe with all nodes in tree
    f = data.frame(c$edge)
    # create a new column with the rate for each branch
    f$rate = attr(c, "rate")
    # keep only terminal branches
    e = f[f$X2<=length(c$tip.label),]
    # add species names
    e$accepted_name = c$tip.label
    
    return(e)
}

rates = mclapply(chronograms, extract.rates, mc.cores = 6)

# merge rate data with LHT data
# separate dataframe for complete cases of each life history trait

height.rate = list()
for (i in 1:200) {
    merged = merge(rates[[i]], height, by = "accepted_name")
    height.rate = c(height.rate, list(merged))
}
write.csv(height.rate[[1]], file = "cp_height_rate.csv")

# perform PGLS on 200 dataframes for each of:
    # rate ~ height
    # rate ~ genome size
    # rate ~ seed mass
    # rate ~ SLA
    # rate ~ all traits

# names in tree must be same as names in dataframe
# check which species are not in dataframe and drop tips from tree

###### HEIGHT ########

tips.to.drop = list()
for (i in 1:100) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(height$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

results.mod.1 = data.frame()
for (i in 1:100) {
    tree = pruned.trees[[i]]
    
    # drop tips from the tree that aren't in the dataset
    tree = drop.tip(tree, tips.to.drop[[i]])
    
    # PGLS using caper
    compdata1 = comparative.data(data = height.rate[[i]], phy = tree, names.col = "accepted_name", vcv = TRUE)
    model1 = pgls(log(rate)~log(mean_height), compdata1, delta = "ML", lambda = 1.0, kappa = 1.0)
    mod1 = summary(model1)
    
    results.mod.1 = rbind(results.mod.1, c(i, mod1$coefficients[2,4], mod1$r.squared))
}     

colnames(results.mod.1) = c("tree", "p_value", "R_squared")
write.csv(results.mod.1, file = "results.height.cp.fossils.csv")


## REPEAT FOR EACH OF THE 4 TRAITS AND FOR ALL TRAITS COMBINED