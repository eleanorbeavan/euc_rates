# a script for:
    # randomly pruning the branches in a phylogeny to >5 substitutions per branch (200 random repeats)
    # estimating rates using chronos
    # merging the trait data with rate data for each species
    # running PGLS on all 200 trees and extracting stats results

# input: a tree, trait dataframe, translation dataframe (phylogeny species to accepted names)
# output: 5 dataframes (one for each model) with p-values and R^2 for all 200 trees

library(ape)
library(phytools)
library(caper)


# read in data
tree = read.tree("your path to the two-gene chloroplast or nuclear tree")
seqdata = read.csv("your path the the alignment names file", StringAsFactors = FALSE)
lhtdata = read.csv("your path to the life-history trait data", StringAsFactors = FALSE)
# keep species with complete height data only
height = lhtdata[complete.cases(lhtdata$mean_height),]

# change the tip.labels in the tree so they match accepted_names in the traits dataframe
tree$tip.label=seqdata$accepted_name[match(tree$tip.label, seqdata$alignment_name)]
tree$tip.label = as.character(tree$tip.label)

# use prune.tree() function. This is found in the 'prune_short_branches_fn' file

# loop to perform "prune.tree" 200 times and create list of 200 pruned trees
pruned.trees = list()
times = 200
for (i in 1:times) {
    newtree = prune.tree(tree, x)
    pruned.trees = c(pruned.trees, list(newtree))
}

#------------------------------------------------------------
# NO FOSSILS
#------------------------------------------------------------
# apply chronos to file with 200 pruned trees to get rates for all 200 trees

# create chronogram for all 200 pruned trees
l = value of lambda found in 'cross_validation_lambda' file
make.chronogram = function(tree) {
    c = chronos(tree, lambda = l, calibration = makeChronosCalib(tree, node="root", age.min =55, age.max = 55, interactive = F))
    return(c)
}

chronograms = mclapply(pruned.trees, make.chronogram, mc.cores = 6)

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
write.csv(rates, file = "name of rate file")

# merge rate data with LHT data
# separate dataframe for complete cases of each life history trait

height.rate = list()
for (i in 1:200) {
    merged = merge(rates[[i]], height, by = "accepted_name")
    height.rate = c(height.rate, list(merged))
}
write.csv(height.rate[[1]], file = "name of rate & trait file")

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
for (i in 1:200) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(height$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

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
write.csv(results.mod.1, file = "name of results file")


## REPEAT FOR EACH OF THE 4 TRAITS AND FOR ALL TRAITS COMBINED

