library(coda)
library(ggplot2)
library(reshape2)
# first remove the '#' from the first row of the trace file by hand...

# then...

d = read.table("path to your file from coevol '.trace'", header = TRUE)
# add a column with generation number
d$generation = 1:nrow(d)
dm = melt(d, id.vars = 'generation')

# choose only the correlations we are interesed in: sub rate with all four traits (2,3,4,5)
i.like = c('sigma_1_2', 'sigma_1_3', 'sigma_1_4', 'sigma_1_5')
# subset based on these four correlations
dm.i.like = dm[which(dm$variable %in% i.like),]

# plot the trace file and visually specify an appropriate burn-in.
# this should be when distribution is stationary
ggplot(dm.i.like, aes(x=generation, y=value)) + geom_line() + facet_wrap(~variable, scales='free_y', ncol=2)

## repeat this for all analyses to choose an appropriate burn-in
