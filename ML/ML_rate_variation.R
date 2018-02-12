# a script for plotting variance in rate data
    # 3 datasets - all with fossils
# average species rates across 200 dataframes
    # plot averages

library(ggplot2)

##### NO FOSSILS #######
rates = read.csv("your path to the rate file (created in ML_analysis.R)")
# take mean of all rate values
rates$mean = rowMeans(subset(rates, select = c(rate,rate.1,rate.2,rate.3,rate.4,rate.5,rate.6,rate.7,rate.8,rate.9,rate.10,rate.11,rate.12,rate.13,rate.14,rate.15,rate.16,rate.17,rate.18,rate.19,rate.20,rate.21,rate.22,rate.23,rate.24,rate.25,rate.26,rate.27,rate.28,rate.29,rate.30,rate.31,rate.32,rate.33,rate.34,rate.35,rate.36,rate.37,rate.38,rate.39,rate.40,rate.41,rate.42,rate.43,rate.44,rate.45,rate.46,rate.47,rate.48,rate.49,rate.50,rate.51,rate.52,rate.53,rate.54,rate.55,rate.56,rate.57,rate.58,rate.59,rate.60,rate.61,rate.62,rate.63,rate.64,rate.65,rate.66,rate.67,rate.68,rate.69,rate.70,rate.71,rate.72,rate.73,rate.74,rate.75,rate.76,rate.77,rate.78,rate.79,rate.80,rate.81,rate.82,rate.83,rate.84,rate.85,rate.86,rate.87,rate.88,rate.89,rate.90,rate.91,rate.92,rate.93,rate.94,rate.95,rate.96,rate.97,rate.98,rate.99)))
rates = rates[,c(5,402)]
# determine the substitution rate per site per year by dividing by 70 million
rates$sr_persite = (rates$mean/70000000)


##### PLOT #######
qplot(rates$sr_persite, geom = "histogram", xlab = "chloroplast/nuclear rate") + scale_x_log10()

# repeat with each dataset (2-gene cp, 2-gene nuc and whole cp)
# repeat with fossil data, dividing by 1 million insstead of 70 million.
