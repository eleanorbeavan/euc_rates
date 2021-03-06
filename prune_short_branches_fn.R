# three functions to 
# 1. find shortest branch length
# 2. remove shortest branch 
# 3. drop all branches <5 substitutions


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
x = (5/"length of alignment(different for each dataset")
prune.tree = function(tree, x) {
    while (shortest.tip(tree) < x) {
        tree = prune.shortest.branch(tree)
    }
    return(tree)
    
new.tree = prune.tree(tree, x)
