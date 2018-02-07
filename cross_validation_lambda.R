## A script for choosing a value of lambda, the rate smoothing parameter
library(ape)

l = c(0, 0.001, 0.01, 1, 10, 1000, 10000)
cv = (length(l))

## 2-gene cp
tree = read.tree("path to 2-gene cp tree")
tree = prune.tree(tree, (5/1767))

cvs = data.frame()
for (i in 1:length(l)) {
    cv = chronopl(tree, lambda = l[i], CV=TRUE, age.min = 1, age.max = 1)
                 
    cvs = rbind(cvs, c(i, sum(attr(cv, "D2")))) 
}

min(cvs$X2522394.63841202) ## lambda = 1

## 2-gene nuc
tree = read.tree("path to 2-gene nuc tree")
tree = prune.tree(tree, (5/1391))

cvs2 = data.frame()
for (i in 1:length(l)) {
    cv = chronopl(tree, lambda = l[i], CV=TRUE, age.min = 1, age.max = 1)
    
    cvs2 = rbind(cvs2, c(i, sum(attr(cv, "D2")))) 
}

min(cvs2$X2522394.63841202) ## all give the same result

## whole cp
tree = read.tree("path to whole cp tree")

cvs3 = data.frame()
for (i in 1:length(l)) {
    cv = chronopl(tree, lambda = l[i], CV=TRUE, age.min = 1, age.max = 1)
    
    cvs3 = rbind(cvs3, c(i, sum(attr(cv, "D2")))) 
}

min(cvs3$X80717211.1889796) ## lambda = 10000

