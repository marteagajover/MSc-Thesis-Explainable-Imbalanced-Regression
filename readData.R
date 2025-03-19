## Instructions

# The working directory must include three archives and one folder:

# - readData.R:    This R script, which calls ReadDataAux.R and DIBSRegress.R functions
#                  to generate the partitions.
# + ReadDataAux.R: R script to generate the partitions in a new folder called "OUTPUT".
#                  It also creates a txt file which contains the percentage of rares.
# + DIBSRegress.R: R script with the SMOGN method.
# + Datasets:      It MUST contain each dataset in KEEL format (*.dat).

################################################


## Installing packages

#install.packages("UBL")
#install.packages("performanceEstimation")
#install.packages("gstat")
#install.packages("intervals")
#install.packages("caret")

## Loading packages

library(UBL)
library(performanceEstimation)
library(gstat)
library(intervals)
require(caret)

### I M P O R T A N T  :  S U B J E C T  T O   C H A N G E (seed) ###

set.seed(12345678)


## Setting working directory

### I M P O R T A N T  :  S U B J E C T  T O   C H A N G E (working directory) ###

setwd("C:/Users/Maria/Desktop/TFM/Experimentos iniciales/Imbalanceo")


# Auxiliar functions (reading datasets)

source("ReadDataAux.R") # generating partitions
source("DIBSRegress.R") # SMOGN method

# Datasets

### I M P O R T A N T  :  S U B J E C T  T O   C H A N G E (datasets directory) ###

# DSlocation = "/Datasets" # directory with 2 datasets
DSlocation = "/zWC-formatooriginal-tra" # directory with 2 datasets

DSnames = list.files(paste(".", DSlocation, sep = ""))

# Reading datasets

### I M P O R T A N T  :  S U B J E C T  T O   C H A N G E (number of datasets) ###

DS = readData(DSnames[[1]], DSlocation) # 2 stands for 2nd dataset

# cv = cvStrat(DS, fold = 10)


############### porcentajes de raros:
percentages = getPercentages(DSnames = DSnames, DSlocation)
fileConn <- file("imbalance_percentage.txt", 'a')
header <- sprintf("%16s", percentages)
write(header, fileConn, append = TRUE)
close(fileConn)


# Stratified, normal and with SMOGN partitions

# For one dataset
cvPartitions(DSnames[1], times = 2, folds = 10, strat = TRUE) # 1 stands for 1st dataset

# For all datasets 
lapply(DSnames, cvPartitions)
