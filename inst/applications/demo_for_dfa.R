library(reshape2)
library(ggplot2)
library(rstan)
library(loo)
library(knitr)
# Load ddata once
load("users/eric.ward/documents/ewidata/data/ewidata")

# Demo for a single EWI dataset - 
# extract data quantifying latitude shifts
lat_data = ewidata[grep(".LAT",ewidata$code),]

# Other data subsets to apply this to:
#blki_data = ewidata[grep("BLKI.",ewidata$code),]
#chi_data = ewidata[grep("CHI.",ewidata$code),]
#pav_data = ewidata[grep("PAV.",ewidata$code),]
#coho_data = ewidata[grep(".CO",ewidata$code),]
#pink_data = ewidata[grep(".PI",ewidata$code),]
#icy_data = ewidata[grep("ICY.",ewidata$code),]
#ich_data = ewidata[grep("ICH.",ewidata$code),]
#juv_data = ewidata[grep("JUVHERRING.",ewidata$code),]

# convert long to wide format
melted = melt(lat_data[,c("code","year","value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y)=="code")])

# generate a data frame of 1:5 trends
dfa_summary = find_dfa_trends(y=Y, iter=1000)

# print the table of LOOIC values and models compared,
kable(dfa_summary$summary)

# plot the estimated rotated trends for the best model
rotated = rotate_trends(dfa_summary$best_model)
plot_trends(rotated)

# plot the estimated loadings for each 
# I have a couple ways we could do this:
plot_loadings(rotated, threshold=0.2, facet=FALSE, names=names)

plot_loadings(rotated, threshold=0.2, facet=TRUE, names=names)
