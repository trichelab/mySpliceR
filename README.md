# install.packages("devtools")
library(devtools)

# Bioconductor needs to be activated when automatically installing dependencies.
options(repos = BiocManager::repositories())
devtools::install_github("trichelab/mySpliceR")

library(mySpliceR)
