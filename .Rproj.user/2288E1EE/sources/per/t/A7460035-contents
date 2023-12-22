## ---------------------------
##
## Script name: table1.R
##
## Purpose of script: To create Table 1 in paper from simulation results
##
## Author: Jeppe Madsen
##
## Date Created: 2023-11-27
##

# Package -----------------------------------------------------------------
library(xtable)

# Read data ---------------------------------------------------------------
d <- readRDS("intermediate_results//simresults.rds")

# Create table ------------------------------------------------------------
beta <- -log(1.5)
Bias <- apply(d[,c(4,7,10,1)], 2, mean) - beta
SD <- apply(d[,c(4,7,10,1)], 2, sd)
RelativeBias <- Bias / SD
AvgSE <- apply(d[,c(2,5,8,11)], 2, mean)
CI_coverage <- apply(d[,c(6,9,12,3)], 2, mean)

# Export table ------------------------------------------------------------
print(xtable(cbind(Bias, SD, RelativeBias, AvgSE, CI_coverage)), file = "results/table1.tex")



