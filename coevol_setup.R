# A script for preparing alignment, phylogeny and trait data to use in coevol
# species present in trait data and alignments must be in phylogeny
# there must be < 100 species present or coevol takes too long to run

library(ape)
library(phytools)
library(phylobase)

# read in data
lhtdata = read.csv("your path to trait data")
seqdata = read.csv("your path to alignment names data")
tree = read.tree("your path to the tree file")
# prior to importing, save alignemnt file as CSV
alignment = read.csv("your path to csv alignment file")

# change all NAs in trait file to -1 (required for Coevol)
lhtdata[is.na(lhtdata)] <- (-1)

# change names in tree to match LHT data
tree$tip.label=seqdata$coevol[match(tree$tip.label, seqdata$alignment_name)]

# remove NAs
cptree = drop.tip(cptree, c(3, 41, 47, 115, 116, 118, 119, 162, 164, 172, 190, 192, 197, 198, 216, 299, 376, 401, 438, 479, 481, 502, 504, 520, 521, 547, 553, 557, 597, 645, 646, 650))

#-------------------------------------------------------------------------------
## FUNCTIONS TO PRUNE SHORTEST BRANCH
# function to find the shortest branch length in the tree
shortest.tip = function(tree) {
    n.tips = length(tree$tip.label)
    descendant.nodes = tree$edge[,2]
    terminal.edges = tree$edge.length[descendant.nodes<=n.tips]
    return(min(terminal.edges))
}

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

# use a loop to remove all branches <5 substitutions
prune.tree = function(tree, x) {
    while (shortest.tip(tree) < x) {
        tree = prune.shortest.branch(tree)
    }
    return(tree)
}
#-------------------------------------------------------------------------------

# TREE
# drop tips with less than 5 substitutions
# x = 5/(length of alignment - different for each dataset)
x = 5/"alignment length"
pruned.tree = prune.tree(tree,x)
# drop species without complete LHT data
# determine which species have complete LHT data
alldata = lhtdata[complete.cases(lhtdata),]
tokeep = as.character(alldata$coevol_name)
# append the outgroup species
tokeep = append(tokeep, "E_papuana")
## use set diff to find species in tree without complete LHT data
tipstodrop = setdiff(pruned.tree$tip.label, tokeep)
tree = drop.tip(pruned.tree, tipstodrop)
# continue to prune species with the shortest branches until 100 species are left
pruned.tree = prune.tree(tree, (14.5/1767))
write.tree(pruned.tree, file = "tree file name")

# LIFE HISTORY TRAITS
speciestokeep = tree$tip.label
pruned.data = lhtdata[lhtdata$coevol_name %in% speciestokeep, ]
write.csv(pruned.data, file = "trait file name")

# ALIGNMENT
pruned.alignment = alignment[alignment$V1 %in% speciestokeep, ]
write.csv(pruned.alignment, file = "alignment file name")
# alignment file needs to be changed back into phylip format in a text editor

# repeat with all data sets
