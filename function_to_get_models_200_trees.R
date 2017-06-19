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
times = 2
for (i in 1:times) {
    newtree = prune.tree(tree)
    output = c(output, list(newtree))
}

# apply chronos to file with 200 pruned trees to get rates for all 200 trees

# create dataframe with node ages
makeChronosCalib(tree, node="root", age.min = 1, age.max = 1, interactive = F)

# create chronogram for all 200 pruned trees
make.chronogram = function(tree) {
    c = chronos(tree, lambda = 1, calibration = makeChronosCalib(tree))
    return(c)
}

chronograms = lapply(pruned.trees, make.chronogram)

