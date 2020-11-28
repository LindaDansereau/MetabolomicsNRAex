# load libraries
library(here)
library(tidyverse)
library(readr)
library(janitor) # fixes column names
library(RColorBrewer)
#library(patchwork)

# function for Standard error calculations
se <- function(x){sd(x)/sqrt(length(x))} 

# concentration calculations
conc_c <- function(y,m,b=0){(y-b)/m}

source("FilteredMetGraph.R") 
# Filters group data based on metabolite list and outputs a bar graph 
# Inputs are the dataframe and metabolite list, eg: highMet

source("ReorderingData.R")
# Reorders data with new levels
# Currently new column is "Name"

source("FilteredMetGraph_Spike.R") 
# Filters individual data based on metabolite list and/or spike level and outputs a bar graph 
# Inputs are the dataframe, spike level, metabolite list, eg: keepMet