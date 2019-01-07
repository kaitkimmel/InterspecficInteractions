## Diversity-Interaction Modeling: BioCON Tissue N ##
# Author: Kaitlin Kimmel
# Date: January, 7 2018

#######################################################################
## Purpose: To determine the extent to which interspecific interactions and environmental
# context modify predicted community weighted traits for tissue N

######################################################################
## Methods: I will be using Diversity- Interaction modeling 
# (refs: Kirwan et al 2009 Ecology & Connolly et al 2013 Journal of Ecology)
# Some of the below code is borrowed from Forest Isbell 

# DI models predict an ecosystem function (in our case, community tissue nitrogen)
# by modeling identity effects, diversity effects, and environmental context effects
# The functional form is: 
# y = sum(BiPi) + sum(dij*(PiPj)^theta) + sum(gammaiPi)*T + sum(gammaij*PiPj)xT
# Bi is the average tissue N in monoculture for species i
# Pi is the proportion of species 1
# dij is the potential for sp i & j to interact such that dijPiPj is the contribution
# of the interaction to the community tissue N
# theta allows the interaction to take on different functional forms, bounded between 0 & 1
# gamma terms x T are how environmental context modifies the identity effect and the diveristy effect

# Pi can either be the planted proportions (as the experimental design is testing)
# or the proportions of the species 1 year prior. It cannot be the proportions of that year
# because that information is captured in the y variable already
# We will test both of these cases. 

# We will also be looking at the model using all 16 species and just the 4 functional groups 
# in the experiment. 

# We will start just using data from the last year [emailed Dan for 2018 data]

###########################################################################
# load libraries
library(here)
###########################################################################
# load data

# Data 
# dataset with % n and % c from aboveground biomass
TisN <- read.csv(here("data", "TisN_clean.csv"),row.names = 1)
names(TisN) <- c("Year", "Date", "Plot", "Ring", "CO2Trt", "NTrt", "SR", "FGR",
                 "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Comments",
                 "Carbon", "Nitrogen", "CNRatio")
Biomass <- read.delim(here("data","Biomass_BioCON.txt"),na.strings=c("","NA"))
names(Biomass) <- c("Sample", "Date", "Plot", "Ring", "CO2Trt", "NTrt", "SR", "FGR",
                    "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Species", "Biomass")
PercCover<- read.delim(here("data","PercCover_BioCON.txt"), na.strings=c("NA",""))
names(PercCover) <- c("Sample", "Season", "Year", "Plot", "Ring", "CO2Trt", "NTrt", "SR", "FGR",
                      "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Species", "PercCov")
CDRSPDat <- read.csv(here::here("data","CDRSPDat.csv"), na.strings=c("NA",""))


##########################################################################
# For the species identity model: I will need to create a dataframe that 
# includes the plot, ring, CO2 Treatment, N treatment, the proportions 
# of the species in the plot, and the measured tissue N
########################################################################

# Subset TisN for last year of experiment, plot, ring, SR, FR, %N, CO2, & N treatments
SIMod.df <- TisN[TisN$ExpYear == max(TisN$ExpYear), c(3:8,16)]
