# _"All_modules.R"_ script notes 

This script has been made to calculate all primary biodegradation values from the 4 modules present in the HC-BioSIM model (water, seawater, soil and sediment environmental compartments)

It contains a source() function, hence make sure the script is within the **working directory** 

This script runs using the fingerprints coming directly from the chemotyper results (.tsv file) and pulls the "Training sets" (.csv files) of all modules from the working directory. Thus, please make sure these datasets and the script are stored within the same **working directory** 

Prior running the script, please make sure the next libraries are installed:

library("rpart")  
library("rpart.plot")  
library("Cubist")  
library("caret")   
library("stringi")  
library("stringr")   
library("scales")  
library("RColorBrewer")  
library("dplyr")   
library("ggplot2")   
library("tidyverse")  
library("plyr")  
library("RCurl")



