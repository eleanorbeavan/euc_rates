# a script for:
    # randomly pruning the branches in a phylogeny to >5 substitutions per branch 
    # estimating rates using chronos
    # merging the trait data with rate data for each species
    # running PGLS and extracting stats results

# input: a tree, trait dataframe, translation dataframe (phylogeny species to accepted names)
# output: 5 dataframes (one for each model) with p-values and R^2 

library(ape)
library(phytools)
library(caper)

# read in data
tree = read.tree("your path to whole chloroplast tree")
lhtdata = read.csv("your path to life-history trait data")
# only keep species with height data
height = lhtdata[complete.cases(lhtdata$mean_height),]

# ONLY USE when you are interested in average trait data
lhtdata = read.csv("your path to average trait data")

# use prune.tree() function. This is found in the 'prune_short_branches_fn' file

# prune branches with <5 substitutions
x = 5/(length of alignment)
pruned.tree = prune.tree(tree, x)

#------------------------------------------------------------
# NO FOSSILS
#------------------------------------------------------------
# create chronogram
# lambda set to 10000
fossils = makeChronosCalib(tree, node="root", age.min = 1, age.max = 1, interactive = F)

chronogram = chronos(pruned.tree, lambda = 10000, calibration = fossils)

# extract rates 
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

rates = extract.rates(chronogram)

# merge rate data with LHT data
# separate dataframe for complete cases of each life history trait

height.rate = merge(rates, height, by = "accepted_name")

# perform PGLS for 5 models:
    # rate ~ height
    # rate ~ genome size
    # rate ~ seed mass
    # rate ~ SLA
    # rate ~ all traits

# names in tree must be same as names in dataframe
# check which species are not in dataframe and drop tips from tree

###### HEIGHT ########

treetips = chronogram$tip.label
datatips = as.character(height$accepted_name)
tips = setdiff(treetips, datatips)
tree1 = drop.tip(tree, tips)

cdat = comparative.data(data = height.rate, phy= tree1, names.col = "accepted_name", vcv = TRUE)
mod1 = pgls(log(rate)~log(mean_height), cdat, delta = "ML", lambda = 1.0, kappa = 1.0)
model1 = summary(mod1)
p_val = model1$coefficients[2,4]
r_squared = model1$r.squared

write.csv(model1, file = "wholecp_height.csv")

# repeat for the other 3 traits plus all traits combined
