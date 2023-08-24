# HC-BioSIM-Biodegradation-QSAR-Public-

This model predicts primary biodegradation (DT50) using a supervised model tree machine-learning model. 

Further details of this model can be found in the Open access paper: 

_"Predicting Primary Biodegradation of Petroleum Hydrocarbons in Aquatic Systems: Integrating System and Molecular Structure Parameters using a Novel Machine-Learning Framework"_

**Available in: https://setac.onlinelibrary.wiley.com/doi/10.1002/etc.5328**

The scripts require the installation of R and RStudio with the next libraries:

        library(rpart)
        library(rpart.plot)
        library(Cubist)
        library(caret)
        library(stringi)
        library(tidyverse)
        library(scales)
        library(dplyr)
        library(plyr)
        library(RcolorBrewer) 

The fingerprints (Toxprints) required for the model must come from the Chemotyper available here: 
https://chemotyper.org/

Reference: 
Yang C, Tarkhov A, Marusczyk J, et al. "New Publicly Available Chemical Query Language, CSRML, 
To Support Chemotype Representations for Application to Data Mining and Modeling." J. Chem. Inf. Model.. 2015;55(3):510-528.
