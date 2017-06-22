# a script for:
    # randomly pruning the branches in a phylogeny to >5 substitutions per branch (200 random repeats)
    # estimating rates using chronos
    # merging the trait data with rate data for each species
    # running PGLS on all 200 trees and extracting stats results

# input: a tree, trait dataframe, translation dataframe (phylogeny species to accepted names)
# output: 5 dataframes (one for each model) with p-values and R^2 for all 200 trees

library(ape)
install.packages("phytools")
library(phytools)
install.packages("phylobase")
library(phylobase)

tree = read.tree("~/Dropbox/euc_sr/old_stuff/chloroplast_tree_one_outgroup.phy")

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
# output?? file with 200 different pruned trees?
# specify where the output needs to go (create dataframe)

# function to remove shortest bl
prune.shortest.branch = function(tree) {
    # take the edge length of those <= the number of tip labels (terminal nodes)
    # vector with terminal branch lengths
    n.tips = length(tree$tip.label)
    descendant.nodes = tree$edge[,2]
    terminal.edges = tree$edge.length[descendant.nodes<=n.tips]
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
        new.tree = prune.shortest.branch(tree)
    }
    return(new.tree)
}

# loop to perform "prune.tree" 200 times and create list of 200 pruned trees
output = list()
times = 200
for (i in 1:times) {
    newtree = prune.tree(tree)
    output = c(output, list(newtree))
}

# apply chronos to file with 200 pruned trees to get rates for all 200 trees
# create chronogram for all 200 pruned trees
make.chronogram = function(tree) {
    c = chronos(tree, lambda = 1, calibration = makeChronosCalib(tree))
    print(c)
    return(c)
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

# creates function to 
    # 1. make 200 chronograms
    # 2. extract rates and create 200 dataframes
get.rates = function(tree, c) {
    chronograms = mclapply(pruned.trees, make.chronogram, mc.cores = 10)
    rates = mclapply(chronograms, extract.rates, mc.cores = 10)
    return(rates)
}

rates = mclapply(pruned.trees[1:3], get.rates, mc.cores = 10)

# merge rate data with LHT data
merge.data = function(rate, trait) {
    merged = merge(rate_data, lhtdata, by = "accepted_name")
    return(merged)
}

LHD = mclapply(rates, merge.data, mc.cores = 10)

# perform PGLS on 200 dataframes for 5 models
    # rate ~ genome size
    # rate ~ height
    # rate ~ genome size
    # rate ~ SLA
    # rate ~ all traits

# names in tree must be same as names in dataframe
# check which species are not in dataframe and drop tips from tree

# perform PGLS on 200 dataframes for 5 models
    # rate ~ genome size
    # rate ~ height
    # rate ~ genome size
    # rate ~ SLA
    # rate ~ all traits

# names in tree must be same as names in dataframe
# check which species are not in dataframe and drop tips from tree
tips.to.drop = list()
for (i in 1:200) {
    treetips = chronograms[[i]]$tip.label
    treetips = na.omit(treetips)
    datatips = as.character(traits$accepted_name)
    tips = setdiff(treetips, datatips)
    tips = append(tips, NA)
    tips.to.drop = c(tips.to.drop, list(tips))
}

mod.1 = list()
results.mod.1 = data.frame(tree_number = "", p_value = "", r_squared = "", StringAsFactors = FALSE)

for (i in 1:200) {
    tree = chronograms[[i]]
    
    # drop tips from the tree that aren't in the dataset
    tree = drop.tip(tree, tips.to.drop[[i]])
    
    # PGLS using caper
    compdata1 = comparative.data(data = traits, phy = tree, names.col = "accepted_name", vcv = TRUE)
    model1 = pgls(log(rate)~log(mean_gs), compdata1, delta = "ML", lambda = 1.0, kappa = 1.0)
    mod1 = summary(model1)
    mod.1 = c(mod.1, list(model1)
              
    results.mod.1 = rbind(results.mod.1, c("i", mod1$coefficients[??], mod1$r.squared) ###### CHECK THIS ######             
 }     
              
   
