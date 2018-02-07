# a script for plotting variance in rate data
    # 3 datasets - all with fossils
# average species rates across 200 dataframes
    # plot averages

library(ggplot2)

##### NO FOSSILS #######
## cp data
cp_rates = read.csv("~/Dropbox/euc_sr/CSV_files/Rates/cp_rates.csv")
# take mean of all rate values
cp_rates$mean = rowMeans(subset(cp_rates, select = c(rate,rate.1,rate.2,rate.3,rate.4,rate.5,rate.6,rate.7,rate.8,rate.9,rate.10,rate.11,rate.12,rate.13,rate.14,rate.15,rate.16,rate.17,rate.18,rate.19,rate.20,rate.21,rate.22,rate.23,rate.24,rate.25,rate.26,rate.27,rate.28,rate.29,rate.30,rate.31,rate.32,rate.33,rate.34,rate.35,rate.36,rate.37,rate.38,rate.39,rate.40,rate.41,rate.42,rate.43,rate.44,rate.45,rate.46,rate.47,rate.48,rate.49,rate.50,rate.51,rate.52,rate.53,rate.54,rate.55,rate.56,rate.57,rate.58,rate.59,rate.60,rate.61,rate.62,rate.63,rate.64,rate.65,rate.66,rate.67,rate.68,rate.69,rate.70,rate.71,rate.72,rate.73,rate.74,rate.75,rate.76,rate.77,rate.78,rate.79,rate.80,rate.81,rate.82,rate.83,rate.84,rate.85,rate.86,rate.87,rate.88,rate.89,rate.90,rate.91,rate.92,rate.93,rate.94,rate.95,rate.96,rate.97,rate.98,rate.99)))
cp_rates = cp_rates[,c(5,402)]
# determine the substitution rate per site per year by dividing by 55 million
cp_rates$sr_persite = (cp_rates$mean/55000000)

## nuc data
nuc_rates = read.csv("~/Dropbox/euc_sr/CSV_files/Rates/nuc_rates.csv")
# take mean of all rate values
nuc_rates$mean = rowMeans(subset(nuc_rates, select = c(rate,rate.1,rate.2,rate.3,rate.4,rate.5,rate.6,rate.7,rate.8,rate.9,rate.10,rate.11,rate.12,rate.13,rate.14,rate.15,rate.16,rate.17,rate.18,rate.19,rate.20,rate.21,rate.22,rate.23,rate.24,rate.25,rate.26,rate.27,rate.28,rate.29,rate.30,rate.31,rate.32,rate.33,rate.34,rate.35,rate.36,rate.37,rate.38,rate.39,rate.40,rate.41,rate.42,rate.43,rate.44,rate.45,rate.46,rate.47,rate.48,rate.49,rate.50,rate.51,rate.52,rate.53,rate.54,rate.55,rate.56,rate.57,rate.58,rate.59,rate.60,rate.61,rate.62,rate.63,rate.64,rate.65,rate.66,rate.67,rate.68,rate.69,rate.70,rate.71,rate.72,rate.73,rate.74,rate.75,rate.76,rate.77,rate.78,rate.79,rate.80,rate.81,rate.82,rate.83,rate.84,rate.85,rate.86,rate.87,rate.88,rate.89,rate.90,rate.91,rate.92,rate.93,rate.94,rate.95,rate.96,rate.97,rate.98,rate.99)))
nuc_mean = nuc_rates[,c(5,402)]
# determine the mutation rate per site per year by dividing by 55 million
nuc_mean$sr_persite = (nuc_mean$mean/55000000)

### whole cp dataset
wholecp_rates = read.csv("~/Dropbox/euc_sr/CSV_files/Rates/rates_wholecp.csv")
wholecp_rates = wholecp_rates[,c(4,5)]
# determine the mutation rate per site per year by dividing by 55 million
wholecp_rates$sr_persite = (wholecp_rates$rate/55000000)

##### LOG PLOTS #######
cplog = qplot(cp_rates$sr_persite, geom = "histogram", xlab = "chloroplast rate") + scale_x_log10()
nuclog = qplot(nuc_mean$sr_persite, geom = "histogram", xlab = "nuclear rate") + scale_x_log10()
wholecplog = qplot(wholecp_rates$sr_persite, geom = "histogram", xlab = "whole cp rate") + scale_x_log10() + scale_x_continuous(breaks = 7e-11)

# repeat with fossil data, dividing by 1 million.