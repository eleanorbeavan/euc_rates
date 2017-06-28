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

options(stringsAsFactors = FALSE)

# read in data
tree = read.tree("~/Dropbox/euc_sr/old_stuff/chloroplast_data/chloroplast_tree_one_outgroup.phy")
seqdata = read.csv("~/Dropbox/euc_sr/CSV_files/alignment_names.csv")
lhtdata = read.csv("~/Dropbox/euc_sr/CSV_files/merged_LHS.csv")
height = lhtdata[complete.cases(lhtdata$mean_height),]
SLA = lhtdata[complete.cases(lhtdata$mean_SLA),]
seedmass = lhtdata[complete.cases(lhtdata$sm),]
genomesize = lhtdata[complete.cases(lhtdata$mean_gs),]
alldata = lhtdata[complete.cases(lhtdata),]

# change the tip.labels in the tree so they match accepted_names in the traits dataframe
tree$tip.label=seqdata$accepted_name[match(tree$tip.label, seqdata$alignment_name)]

is.rooted(tree)

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

"multiPhylo" = class(pruned.trees)

#------------------------------------------------------------

# apply chronos to file with 200 pruned trees to get rates for all 200 trees

# create dataframe with node ages
makeChronosCalib(tree, node="root", age.min = 1, age.max = 1, interactive = F)

# create chronogram for all 200 pruned trees
make.chronogram = function(tree) {
    c = chronos(tree, lambda = 1, calibration = makeChronosCalib(tree))
    return(c)
}

chronograms = mclapply(pruned.trees, make.chronogram, mc.cores = 10)

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

rates = mclapply(chronograms, extract.rates, mc.cores = 10)

########################
# creates function to 
    # 1. make 200 chronograms
    # 2. extract rates and create 200 dataframes
get.rates = function(tree, c) {
    chronograms = mclapply(pruned.trees, make.chronogram, mc.cores = 10)
    rates = mclapply(chronograms, extract.rates, mc.cores = 10)
    return(rates)
}

rates = mclapply(pruned.trees, get.rates, mc.cores = 10)
########################

# merge rate data with LHT data
# separate dataframe for complete cases of each life history trait

height.rate = list()
for (i in 1:200) {
    merged = merge(rates[[i]], height, by = "accepted_name")
    height.rate = c(height.rate, list(merged))
}

gs.rate = list()
for (i in 1:200) {
    merged = merge(rates[[i]], genomesize, by = "accepted_name")
    gs.rate = c(gs.rate, list(merged))
}

sm.rate = list()
for (i in 1:200) {
    merged = merge(rates[[i]], seedmass, by = "accepted_name")
    sm.rate = c(sm.rate, list(merged))
}

sla.rate = list()
for (i in 1:200) {
    merged = merge(rates[[i]], SLA, by = "accepted_name")
    sla.rate = c(sla.rate, list(merged))
}

LHD.rate = list()
for (i in 1:200) {
    merged = merge(rates[[i]], alldata, by = "accepted_name")
    LHD.rate = c(LHD.rate, list(merged))
}

# perform PGLS on 200 dataframes for 5 models
    # rate ~ height
    # rate ~ genome size
    # rate ~ seed mass
    # rate ~ SLA
    # rate ~ all traits

# names in tree must be same as names in dataframe
# check which species are not in dataframe and drop tips from tree
###### HEIGHT ########

tips.to.drop = list()
for (i in 1:200) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(height$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

class(chronograms) <- "multiPhylo"

mod.1 = list()
results.mod.1 = data.frame()

for (i in 1:200) {
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
write.csv(results.mod.1, file = "results.mod.1.csv")

###### GENOME SIZE ########

tips.to.drop = list()
for (i in 1:200) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(genomesize$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

mod.2 = list()
results.mod.2 = data.frame()

for (i in 1:200) {
    tree = pruned.trees[[i]]
    
    # drop tips from the tree that aren't in the dataset
    tree = drop.tip(tree, tips.to.drop[[i]])
    
    # PGLS using caper
    compdata2 = comparative.data(data = gs.rate[[i]], phy = tree, names.col = "accepted_name", vcv = TRUE)
    model2 = pgls(log(rate)~log(mean_gs), compdata2, delta = "ML", lambda = 1.0, kappa = 1.0)
    mod2 = summary(model2)
    
    results.mod.2 = rbind(results.mod.2, c(i, mod2$coefficients[2,4], mod2$r.squared))
}     

colnames(results.mod.2) = c("tree", "p_value", "R_squared")
write.csv(results.mod.2, file = "results.mod.2.csv")

###### SEED MASS ########

tips.to.drop = list()
for (i in 1:200) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(seedmass$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

mod.3 = list()
results.mod.3 = data.frame()

for (i in 1:200) {
    tree = pruned.trees[[i]]
    
    # drop tips from the tree that aren't in the dataset
    tree = drop.tip(tree, tips.to.drop[[i]])
    
    # PGLS using caper
    compdata3 = comparative.data(data = sm.rate[[i]], phy = tree, names.col = "accepted_name", vcv = TRUE)
    model3 = pgls(log(rate)~log(sm), compdata3, delta = "ML", lambda = 1.0, kappa = 1.0)
    mod3 = summary(model3)
    
    results.mod.3 = rbind(results.mod.3, c(i, mod3$coefficients[2,4], mod3$r.squared))
}     

colnames(results.mod.3) = c("tree", "p_value", "R_squared")
write.csv(results.mod.3, file = "results.mod.3.csv")

###### SLA ########

tips.to.drop = list()
for (i in 1:200) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(SLA$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

mod.4 = list()
results.mod.4 = data.frame()

for (i in 1:200) {
    tree = pruned.trees[[i]]
    
    # drop tips from the tree that aren't in the dataset
    tree = drop.tip(tree, tips.to.drop[[i]])
    
    # PGLS using caper
    compdata4 = comparative.data(data = sla.rate[[i]], phy = tree, names.col = "accepted_name", vcv = TRUE)
    model4 = pgls(log(rate)~log(mean_SLA), compdata4, delta = "ML", lambda = 1.0, kappa = 1.0)
    mod4 = summary(model4)
    
    results.mod.4 = rbind(results.mod.4, c(i, mod4$coefficients[2,4], mod4$r.squared))
}     

colnames(results.mod.4) = c("tree", "p_value", "R_squared")
write.csv(results.mod.4, file = "results.mod.4.csv")

###### ALL DATA ########

tips.to.drop = list()
for (i in 1:200) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(alldata$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

mod.5 = list()
results.mod.5 = data.frame()

for (i in 1:200) {
    tree = pruned.trees[[i]]
    
    # drop tips from the tree that aren't in the dataset
    tree = drop.tip(tree, tips.to.drop[[i]])
    
    # PGLS using caper
    compdata5 = comparative.data(data = LHD.rate[[i]], phy = tree, names.col = "accepted_name", vcv = TRUE)
    model5 = pgls(log(rate)~log(mean_height)*log(sm)*log(mean_SLA)*log(mean_gs), compdata5, delta = "ML", lambda = 1.0, kappa = 1.0)
    mod5 = summary(model5)
    
    results.mod.5 = rbind(results.mod.5, c(i, mod5$coefficients[6,4], mod5$coefficients[7,4], mod5$coefficients[8,4], mod5$coefficients[9,4], mod5$coefficients[10,4], mod5$coefficients[11,4], mod5$r.squared)) ## will need to add more coeffs
}     

colnames(results.mod.5) = c("tree", "p_value_h*sm", "p_value_h*sla", "p_value_sla*sm", "p_value_h*gs", "p_value_gs*sm", "p_value_sla*gs", "R_squared")
write.csv(results.mod.5, file = "results.mod.5.csv")
