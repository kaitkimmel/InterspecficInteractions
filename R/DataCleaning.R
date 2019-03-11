# Data cleaning: 
library(dplyr)
library(here)

## 1994-2016 data 

TisN <- read.csv(here("data", "PercN.csv"), na.strings=c("","NA"))
names(TisN) <- c("Year", "Date", "Plot", "Ring", "CO2", "N", "SR", "FGR",
                 "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Comments",
                 "Carbon", "Nitrogen", "CNRatio")
CDRsp <- read.csv(here("data", "CDRSPDat.csv"))
## Data cleaning - Species names ##
# Need to change species names to Genus species that are in GenusSpecies format
TisN$monospecies<- gsub(pattern = "AchilleaMillefolium", replacement = "Achillea millefolium", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "BoutelouaGracilis", replacement = "Bouteloua gracilis", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "AsclepiasTuberosa", replacement = "Asclepias tuberosa", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "SchizachyriumScoparium", replacement = "Schizachyrium scoparium", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "AmorphaCanescens", replacement = "Amorpha canescens", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "BromusInermis", replacement = "Bromus inermis", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "AgropyronRepens", replacement = "Agropyron repens", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "LespedezaCapitata", replacement = "Lespedeza capitata", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "PetalostemumVillosum", replacement = "Petalostemum villosum", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "PoaPratensis", replacement = "Poa pratensis", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "SolidagoRigida", replacement = "Solidago rigida", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "KoeleriaCristata", replacement = "Koeleria cristata", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "LupinusPerennis", replacement = "Lupinus perennis", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "AndropogonGerardi", replacement = "Andropogon gerardi", TisN$monospecies)
TisN$monospecies<- gsub(pattern = "SorghastrumNutans", replacement = "Sorghastrum nutans", TisN$monospecies)

## Data cleaning - getting columns and rows correct ##
TisN <- TisN[-which(is.na(TisN$Nitrogen)),] # get rid of rows without tissue N values
TisN <- subset(TisN, Comments == "Total Above") # Aboveground tissue N only
TisN$ExpYear <- TisN$Year-1997 #Create experiment year column
TisN$l.year <- log(TisN$ExpYear) #Log year column
TisN$YearSq <- TisN$ExpYear^2 # Year squared column

# Get rid of drought plots
TisN <- TisN[-which(TisN$WaterTrt == 'H2Oneg'),]
# Get rid of warmed plots; some years Ht v HT
TisN$TempTrt <- gsub(pattern = "Htelv", replacement = "HTelv", TisN$TempTrt)
TisN$TempTrt <- gsub(pattern = "Htamb", replacement = "HTamb", TisN$TempTrt)
TisN <- TisN[-which(TisN$TempTrt == 'HTelv'),]

# save cleaned dataset
write.csv(TisN, here("data", "TisN_clean.csv"))

# monoplots
monoplots <- TisN[which(TisN$SR ==1),c(3,10)]
monoplots <- unique(monoplots)
colnames(CDRsp)[1] <- "monospecies"
CDRsp$monospecies <- as.character(CDRsp$monospecies)
CDRsp$monospecies[5] <- "Achillea millefolium"
monoplots <- merge(monoplots, CDRsp[,c(1,2)])
write.csv(monoplots, here("data", "monoplots.csv"))

## 2017 data

# load dataset
TisN17 <- read.csv(here("data", "CN_dat2017.csv"))

# Get rid of drought and elevated temp plots 

TisN17 <- TisN17[TisN17$Water.Treatment != "H2Oneg",]
TisN17 <- TisN17[TisN17$Temp.Treatment != "HTelv",]

# Get rid of rows with no % N values
TisN17 <- TisN17[-which(is.na(TisN17$Aboveground.Nitrogen...)),]

# Subset out only columns needed
# Plot, Ring, CO2 Treat, N Treat, SR, %N
TisN17 <- TisN17[,c(4:7,12,18)]
names(TisN17) <- c("Ring", "Plot",  "CO2", "N", "SR", "Nitrogen")

# Save cleaned data

write.csv(TisN17, here("data", "TisN17_clean.csv"))


# Combine 2017 with rest
TisN17$Year = as.numeric(2017)
plotsub <- unique(TisN17$Plot) # want to look at the same set of plots through time
TisN <- TisN[TisN$Plot %in% plotsub,]
TisN_sub <- TisN[,c(1,3:7,16)]

totdat <- bind_rows(TisN17, TisN_sub)
totdat$ExpYear <- totdat$Year-1997

write.csv(totdat,here("data", "total_clean.csv"), row.names = TRUE)
