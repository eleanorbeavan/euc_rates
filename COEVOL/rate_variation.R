## a script for extracting substitution rate estimates from coevol output 
  # and plotting them with 95% confidence intervals

library(ggplot2)

## when running 'readcoevol' specify +mean as an option to obtain mean substitution rate estimates in the '.postmeanbranchsynrate.tab' file

# read in data
rates = read.table(file = "path to your coevol output file .postmeanbranchsynrate.tab")
# keep only the terminal branch rates (because this is for comparison to ML analysis, which only uses terminal branches)
rates = rates[rates$V1==rates$V2,]
# column 4 is mean sub rate per branch
# divide rate value by 70 million to get rate/site/year (divide by 1 million if fossil calibrations were used)
rates$V7 = (rates$V4/70000000)
# divide upper and lower confidence intervals by 70 million (or 1 million if fossils were used)
rates$V8 = (rates$V5/70000000)
rates$V9 = (rates$V6/70000000)

# plot the rate values with 95% confidence intervals
ggplot(data = rates, aes(x=V1, y=V7)) + geom_point() + geom_errorbar(aes(ymin = V8, ymax = V9)) + xlab("Species") + ylab("Rate") + theme(axis.text.x = element_blank()) 

## repeat this for all analyses
