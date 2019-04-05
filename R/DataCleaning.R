# Data cleaning: 
library(dplyr)
library(here)


## Tissue Nitrogen and Carbon 1994-2016 data 

TissueNut <- read.csv(here("data", "PercN.csv"), na.strings=c("","NA"))
names(TissueNut) <- c("Year", "Date", "Plot", "Ring", "CO2", "N", "SR", "FGR",
                 "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Comments",
                 "Carbon", "Nitrogen", "CNRatio", "Notes")
CDRsp <- read.csv(here("data", "CDRSPDat.csv"))
## Data cleaning - Species names ##
# Need to change species names to Genus species that are in GenusSpecies format
TissueNut$monospecies<- gsub(pattern = "AchilleaMillefolium", replacement = "Achillea millefolium", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "BoutelouaGracilis", replacement = "Bouteloua gracilis", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "AsclepiasTuberosa", replacement = "Asclepias tuberosa", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "SchizachyriumScoparium", replacement = "Schizachyrium scoparium", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "AmorphaCanescens", replacement = "Amorpha canescens", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "BromusInermis", replacement = "Bromus inermis", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "AgropyronRepens", replacement = "Agropyron repens", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "LespedezaCapitata", replacement = "Lespedeza capitata", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "PetalostemumVillosum", replacement = "Petalostemum villosum", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "PoaPratensis", replacement = "Poa pratensis", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "SolidagoRigida", replacement = "Solidago rigida", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "KoeleriaCristata", replacement = "Koeleria cristata", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "LupinusPerennis", replacement = "Lupinus perennis", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "AndropogonGerardi", replacement = "Andropogon gerardi", TissueNut$monospecies)
TissueNut$monospecies<- gsub(pattern = "SorghastrumNutans", replacement = "Sorghastrum nutans", TissueNut$monospecies)

## Data cleaning - getting columns and rows correct ##
TissueNut <- TissueNut[-which(is.na(TissueNut$Nitrogen)),] # get rid of rows without tissue N values
TissueNut <- subset(TissueNut, Comments == "Total Above") # Aboveground tissue nutrients only
TissueNut$ExpYear <- TissueNut$Year-1997 #Create experiment year column
TissueNut$l.year <- log(TissueNut$ExpYear) #Log year column
TissueNut$YearSq <- TissueNut$ExpYear^2 # Year squared column

# Get rid of drought plots
watersub <- unique(TissueNut[which(TissueNut$WaterTrt == 'H2Oneg'), "Plot"])
# Get rid of warmed plots; some years Ht v HT
TissueNut$TempTrt <- gsub(pattern = "Htelv", replacement = "HTelv", TissueNut$TempTrt)
TissueNut$TempTrt <- gsub(pattern = "Htamb", replacement = "HTamb", TissueNut$TempTrt)
tempsub <- unique(TissueNut[which(TissueNut$TempTrt == 'HTelv'),"Plot"])
trtsub <- unique(c(watersub, tempsub))

TissueNut <- TissueNut[-which(TissueNut$Plot %in% trtsub),]

# save cleaned dataset
write.csv(TissueNut, here("data", "TissueNut_clean.csv"))

# monoplots
monoplots <- TissueNut[which(TissueNut$SR ==1),c(3,10)]
monoplots <- unique(monoplots)
colnames(CDRsp)[1] <- "monospecies"
CDRsp$monospecies <- as.character(CDRsp$monospecies)
CDRsp$monospecies[5] <- "Achillea millefolium"
monoplots <- merge(monoplots, CDRsp[,c(1,2)])
write.csv(monoplots, here("data", "monoplots.csv"))

## 2017 Tissue Nutrient data

# load dataset
TissueNut17 <- read.csv(here("data", "CN_dat2017.csv"))

# Get rid of drought and elevated temp plots 

TissueNut17 <- TissueNut17[-which(TissueNut17$Plot %in%trtsub),]



# Get rid of rows with no % N values
TissueNut17 <- TissueNut17[-which(is.na(TissueNut17$Aboveground.Nitrogen...)),]

# Subset out only columns needed
# Plot, Ring, CO2 Treat, N Treat, SR, %N
TissueNut17 <- TissueNut17[,c(4:7,12,18:20)]
names(TissueNut17) <- c("Ring", "Plot",  "CO2", "N", "SR", "Nitrogen", "Carbon", "CNRatio")

# Save cleaned data

write.csv(TissueNut17, here("data", "TissueNut17_clean.csv"))


# Combine 2017 with rest
TissueNut17$Year = as.numeric(2017)
TissueNut_sub <- TissueNut[,c(1,3:7,15:17)]

TissueNut_dat <- bind_rows(TissueNut17, TissueNut_sub)
TissueNut_dat$ExpYear <- TissueNut_dat$Year-1997

write.csv(TissueNut_dat,here("data", "TissueNut_data.csv"), row.names = TRUE)


### Root Nutrients

rootNut <- read.csv (here("data", "RootNut.csv"), na.strings=c("","NA"))
names(rootNut) <- c("Year", "Date", "Plot", "Ring", "CO2", "N", "SR", "FGR",
                      "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "RootSampleDepth",
                      "Root.Carbon", "Root.Nitrogen", "Root.CNRatio")
# Need to change species names to Genus species that are in GenusSpecies format
rootNut$monospecies<- gsub(pattern = "AchilleaMillefolium", replacement = "Achillea millefolium", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "BoutelouaGracilis", replacement = "Bouteloua gracilis", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "AsclepiasTuberosa", replacement = "Asclepias tuberosa", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "SchizachyriumScoparium", replacement = "Schizachyrium scoparium", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "AmorphaCanescens", replacement = "Amorpha canescens", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "BromusInermis", replacement = "Bromus inermis", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "AgropyronRepens", replacement = "Agropyron repens", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "LespedezaCapitata", replacement = "Lespedeza capitata", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "PetalostemumVillosum", replacement = "Petalostemum villosum", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "PoaPratensis", replacement = "Poa pratensis", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "SolidagoRigida", replacement = "Solidago rigida", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "KoeleriaCristata", replacement = "Koeleria cristata", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "LupinusPerennis", replacement = "Lupinus perennis", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "AndropogonGerardi", replacement = "Andropogon gerardi", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "SorghastrumNutans", replacement = "Sorghastrum nutans", rootNut$monospecies)
rootNut$monospecies<- gsub(pattern = "AnemoneCylindrica", replacement = "Anemone cylindrica", rootNut$monospecies)
rootNut$RootSampleDepth<- gsub(pattern = "0-20cm ", replacement = "0-20cm", rootNut$RootSampleDepth)
rootNut <- rootNut[-which(rootNut$Plot %in%trtsub),]
rootNut <- rootNut[rootNut$RootSampleDepth == "0-20cm",]
rootNut <- rootNut[,c(1,3:7,10,15:17)]

write.csv(rootNut,here("data", "rootNutClean.csv"))

### Amax, SLA, gs
AmaxEtc <- read.csv(here("data","Amax_SLA_gs.csv"), na.strings = c("", "NA"))
names(AmaxEtc) <- c("Year", "Plot", "Ring", "CO2", "N", "SR", "Species",
                    "FG", "Rep", "Date", "Time", "AtGrowth", "PAR",
                    "Trans", "Amax.area", "tleaf", "gs", "A.gs", "SLA", "Amax.mass",
                    "system.used")
AmaxEtc <- AmaxEtc[AmaxEtc$AtGrowth == "atgrwth",]
AmaxEtc<- AmaxEtc[is.na(AmaxEtc$Rep) | AmaxEtc$Rep == 1,]
AmaxEtc <- AmaxEtc[-which(AmaxEtc$Plot %in% trtsub),]
AmaxEtc <- AmaxEtc[,c(1:7,14,15,17:20)]
AmaxEtc$SLA<- gsub(pattern = "#DIV/0!", replacement = "NA", AmaxEtc$SLA)
AmaxEtc$SLA <- as.numeric(AmaxEtc$SLA)
write.csv(AmaxEtc,here("data", "AmaxEtcClean.csv"))

### Light, I*

light <- read.csv(here("data", "Light.csv"), na.strings = c("", "NA"))
names(light) <- c("SampleNum", "Year", "Date", "Time", "Plot", "Ring", "CO2", "N", "SR",
                  "FGR", "Exp", "monospecies", "monoFunGroup","WaterTrt", "TempTrt", "Clip.Strip", "propTrans",
                  "Notes")
light <- light[-which(light$Plot %in% trtsub),]
light$Date <- as.Date(light$Date,"%m/%d/%Y")
light$Month <- as.numeric(format(light$Date, format = "%m"))
light <- light[light$Month == 8,] # August data
light <- aggregate(light$propTrans, by = list(Year = light$Year, Plot = light$Plot, CO2 = light$CO2, N = light$N),
                  FUN = mean)
names(light)[5]<-"propLightTrans"
write.csv(light, here("data", "lightClean.csv"))

## Soil Nitrate , R*
soilN <- read.csv(here("data", "SoilNitrate.csv"))#, na.strings = c("", "NA"))
soilN <- soilN[,-c(15:17)]
soilN <- soilN[-which(is.na(soilN[,c(1:14)])),]
names(soilN) <- c("SampleNum", "Date", "Plot", "Ring", "CO2", "N", "SR",
                  "FGR", "Exp", "monospecies", "monoFunGroup","Depth", 
                  "Nitrate", "Ammonium")
soilN <- soilN[-which(soilN$Plot %in% trtsub),]
soilN$Date <- as.Date(soilN$Date,"%m/%d/%y")
soilN$Month <- as.numeric(format(soilN$Date, format = "%m"))
soilN$Year <- as.numeric(format(soilN$Date, format = "%Y"))
soilN <- soilN[soilN$Month == 8,]
soilN <- soilN[soilN$Depth == "0-20",]
soilN <- soilN[,c(3:7,10,13,16)]
write.csv(soilN, here("data", "soilNClean.csv"))

## Seed weight

seed <- read.delim(here("data", "txt files", "extraseedwt.txt")) # using data from big bio to keep consistent
seed <- seed[,c(2,7)]
#names(seed) <- c("Year", "Plot", "Species", "Ring", "CO2", "N", "SR", "FGR", "Exp", 
#                 "monospecies", "MonoFunGroup", "Stalk", "FlowersOnStalk", "Flower", 
#                 "SeedCount", "SeedWt", "AvgSeedWt", "FlowerWt")
#seed <- seed[-which(seed$Plot %in% trtsub),]
#seed <- aggregate(seed$SeedWt, by = list(Year = seed$Year, Plot = seed$Plot, X5Lspecid = seed$Species, CO2 = seed$CO2, N = seed$N),
#                  FUN = mean)
names(seed) <- c("monospecies", "SeedWt")
#seed <- merge(seed,CDRsp[,c(1,8)])
#seed <- seed[,-1]
seed <- merge(monoplots, seed)
write.csv(seed, here("data", "seedWtClean.csv"))
seed$Year <- 2008

## Root Biomass
rootBiomass <- read.csv(here("data", "RootBiomass.csv"), na.strings = c("","NA"))
names(rootBiomass) <- c("Sample", "Date", "Plot", "Ring", "CO2", "N", "SR", "FGR", 
                        "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt",
                        "Depth", "SubSamp", "Mass")
rootBiomass$SubSamp<- gsub(pattern = "Crowns", replacement = "Crown", rootBiomass$SubSamp)
rootBiomass <- rootBiomass[rootBiomass$SubSamp != "Crown",]
rootBiomass <- rootBiomass[-which(rootBiomass$Plot %in% trtsub),]
rootBiomass$Date <- as.Date(rootBiomass$Date,"%m/%d/%y")
rootBiomass$Month <- as.numeric(format(rootBiomass$Date, format = "%m"))
rootBiomass$Year <- as.numeric(format(rootBiomass$Date, format = "%Y"))
rootBiomass <- rootBiomass[rootBiomass$Month == "8",] # Get august data only
rootBiomass$Ring  <- NULL
# Need to change species names to Genus species that are in GenusSpecies format
rootBiomass$monospecies<- gsub(pattern = "AchilleaMillefolium", replacement = "Achillea millefolium", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "BoutelouaGracilis", replacement = "Bouteloua gracilis", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "AsclepiasTuberosa", replacement = "Asclepias tuberosa", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "SchizachyriumScoparium", replacement = "Schizachyrium scoparium", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "AmorphaCanescens", replacement = "Amorpha canescens", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "BromusInermis", replacement = "Bromus inermis", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "AgropyronRepens", replacement = "Agropyron repens", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "LespedezaCapitata", replacement = "Lespedeza capitata", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "PetalostemumVillosum", replacement = "Petalostemum villosum", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "PoaPratensis", replacement = "Poa pratensis", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "SolidagoRigida", replacement = "Solidago rigida", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "KoeleriaCristata", replacement = "Koeleria cristata", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "LupinusPerennis", replacement = "Lupinus perennis", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "AndropogonGerardi", replacement = "Andropogon gerardi", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "SorghastrumNutans", replacement = "Sorghastrum nutans", rootBiomass$monospecies)
rootBiomass$monospecies<- gsub(pattern = "AnemoneCylindrica", replacement = "Anemone cylindrica", rootBiomass$monospecies)
rootBiomass$Depth<- gsub(pattern = "0-20 ", replacement = "0-20", rootBiomass$Depth)

rootBiomass <- rootBiomass[rootBiomass$SubSamp != "Unsorted", ]

# Calculate total biomass over all depth
TotalRootBiomass <- aggregate(rootBiomass$Mass, 
                              by = list(Plot = rootBiomass$Plot,
                                        Year = rootBiomass$Year), FUN = sum)
names(TotalRootBiomass)[3] <- "TotalRootBiomass"
# Calfulate Fine root biomass in shallow
FineRootBiomass <- rootBiomass[rootBiomass$SubSamp == "Fine" & rootBiomass$Depth == "0-20", ]
FineRootBiomass <- aggregate(FineRootBiomass$Mass, 
                              by = list(Plot = FineRootBiomass$Plot,
                                        Year = FineRootBiomass$Year), FUN = sum)
names(FineRootBiomass)[3] <- "FineRootBiomass"
# Calculate total shallow biomass
ShallowRootBiomass <- rootBiomass[rootBiomass$Depth == "0-20",]
ShallowRootBiomass <- aggregate(ShallowRootBiomass$Mass, 
                             by = list(Plot = ShallowRootBiomass$Plot,
                                       Year = ShallowRootBiomass$Year), FUN = sum)
names(ShallowRootBiomass)[3] <- "ShallowRootBiomass"
RootBiomass <- plyr::join_all(list(FineRootBiomass, TotalRootBiomass, ShallowRootBiomass))
RootBiomass$FineRootAllo <- RootBiomass$FineRootBiomass/RootBiomass$ShallowRootBiomass
RootBiomass$PropShallow <- RootBiomass$ShallowRootBiomass/ RootBiomass$TotalRootBiomass
# Root biomass is only measured in every plot for years 2000, 2001, 2004, 2005, 2006
# Setting proportion shallow to NA for other years
yearsub <- c(2000,2001,2004,2005,2006)
RootBiomass$PropShallow[-which(RootBiomass$Year %in% yearsub)] <- NA
write.csv(RootBiomass, here("data", "RootsClean.csv"))

## Biomass data

ab.mat <- read.delim(here("data","txt files","Biomass_BioCON.txt"),na.strings=c("","NA"))
names(ab.mat) <- c("Sample", "Date", "Plot", "Ring", "CO2", "N", "SR", "FGR",
                    "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Species", "Biomass")

# Fix species names
ab.mat$Species<- gsub(pattern = "Achillea millefolium ", replacement = "Achillea millefolium", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Bouteloua gracilis ", replacement = "Bouteloua gracilis", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Asclepias tuberosa ", replacement = "Asclepias tuberosa", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Schizachyrium scoparium ", replacement = "Schizachyrium scoparium", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Amorpha canescens ", replacement = "Amorpha canescens", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Bromus inermis ", replacement = "Bromus inermis", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Agropyron repens ", replacement = "Agropyron repens", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Lespedeza capitata ", replacement = "Lespedeza capitata", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Petalostemum villosum ", replacement = "Petalostemum villosum", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Poa pratensis ", replacement = "Poa pratensis", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "poa pratensis", replacement = "Poa pratensis", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Solidago rigida ", replacement = "Solidago rigida", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Koeleria cristata ", replacement = "Koeleria cristata", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Lupinus perennis ", replacement = "Lupinus perennis", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Andropogon gerardi ", replacement = "Andropogon gerardi", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Sorghastrum nutans ", replacement = "Sorghastrum nutans", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "Anemone cylindrica ", replacement = "Anemone cylindrica", ab.mat$Species)
ab.mat$Species<- gsub(pattern = "bromus inermis", replacement = "Bromus inermis", ab.mat$Species)
ab.mat$monospecies<- gsub(pattern = "AchilleaMillefolium", replacement = "Achillea millefolium", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "BoutelouaGracilis", replacement = "Bouteloua gracilis", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "AsclepiasTuberosa", replacement = "Asclepias tuberosa", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "SchizachyriumScoparium", replacement = "Schizachyrium scoparium", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "Amorpha  canescens", replacement = "Amorpha canescens", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "AmorphaCanescens", replacement = "Amorpha canescens", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "BromusInermis", replacement = "Bromus inermis", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "AgropyronRepens", replacement = "Agropyron repens", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "LespedezaCapitata", replacement = "Lespedeza capitata", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "PetalostemumVillosum", replacement = "Petalostemum villosum", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "PoaPratensis", replacement = "Poa pratensis", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "Solidago sigida", replacement = "Solidago rigida", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "SolidagoRigida", replacement = "Solidago rigida", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "KoeleriaCristata", replacement = "Koeleria cristata", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "LupinusPerennis", replacement = "Lupinus perennis", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "AndropogonGerardi", replacement = "Andropogon gerardi", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "SorghastrumNutans", replacement = "Sorghastrum nutans", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "AnemoneCylindrica", replacement = "Anemone cylindrica", ab.mat$monospecies)
ab.mat$monospecies<- gsub(pattern = "bromusInermis", replacement = "Bromus inermis", ab.mat$monospecies)
# Change Cenriched to Cenrich
ab.mat$CO2 <- gsub(pattern = "Cenriched", replacement = "Cenrich", ab.mat$CO2)
ab.mat$CO2 <- gsub(pattern = "Cenrich ", replacement = "Cenrich", ab.mat$CO2)
# Get only August data
ab.mat$Date <- as.Date(ab.mat$Date,"%m/%d/%Y")
ab.mat$Year <- as.numeric(format(ab.mat$Dat, format = "%Y"))
ab.mat$Month <- as.numeric(format(ab.mat$Date, format = "%m"))
ab.mat <- ab.mat[-which(ab.mat$Month == 6),]
ab.mat$ExpYear <- ab.mat$Year - 1997

ab.mat <- ab.mat[-which(ab.mat$Species == 'Unsorted Biomass'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Miscellaneous Litter'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Miscellaneous litter'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Real Weeds'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Real Weeds '),]
ab.mat <- ab.mat[-which(ab.mat$Species == '16 Species Weeds'),]
ab.mat <- ab.mat[-which(ab.mat$Species == '16 Species Weeds '),]

# get same plots used in totdat - no elevated temp or drought
ab.mat <- ab.mat[-which(ab.mat$Plot %in% trtsub),]

# Only considering species originally planted 
###### Check if this is ok by looking at perc. biomass of these species
sp <- unique(monoplots$monospecies)
# Biomass of planted species 
planted <- ab.mat[which(ab.mat$Species %in% sp),]
planted_bio <- aggregate(planted$Biomass, by = list(Plot = planted$Plot, Year = planted$Year), FUN = sum)
# Get total Biomass for each plot and year
tot_bio <- aggregate(ab.mat$Biomass, by = list(Plot = ab.mat$Plot, Year = ab.mat$Year), FUN = sum)
colnames(tot_bio)[3] <- "TotalBiomass"
x <- merge(tot_bio, planted_bio)
x$diff <- (x$x- x$TotalBiomass)/x$TotalBiomass
length(which(x$diff > -0.1))/nrow(tot_bio) ## 99% of plots have 90% of biomass from 16 species in experiment

## Only consider 16 species for rest of analysis ##
ab.mat <- ab.mat[which(ab.mat$Species %in% sp),]
planted_bio <- aggregate(planted$Biomass, by = list(Plot = planted$Plot, Year = planted$Year), FUN = sum)
colnames(planted_bio)[3] <- "TotalBiomass"
ab.mat <- merge(ab.mat, planted_bio)
ab.mat <- ab.mat[-20617,]
biomass.mat <- tidyr::spread(data = ab.mat, key = Species, fill = 0, value = Biomass)
write.csv(biomass.mat, here::here("data", "biomass_mat.csv"), row.names = TRUE)

ab.mat$propbio <- ab.mat$Biomass/ab.mat$TotalBiomass
ab.mat.sub <- ab.mat[,c(1,2,5:9,15,18,20)]
# Duplicate identifiers for rows (20383, 20384)
#ab.mat.sub <- ab.mat.sub[-20617,]
ab.mat.wide <- tidyr::spread(data = ab.mat.sub, key = Species, fill = 0, value = propbio)


write.csv(ab.mat.wide, here::here("data", "abund_mat.csv"), row.names = TRUE)


# Combine all datasets for monocultures

mono.dat <- plyr::join_all(list(monoplots, rootNut, RootBiomass , TissueNut_dat, soilN, seed), type = "left")
hold <- plyr::join_all(list(monoplots, AmaxEtc))
hold1 <- plyr::join_all(list(monoplots, light))
mono.dat <- plyr::join(hold,mono.dat, by = c("Plot", "Year"), type = "full")
mono.dat <- plyr::join(hold1,mono.dat, by = c("Plot", "Year"), type = "full")
mono.dat <- plyr::join(mono.dat, planted_bio, type = "left")  
mono.dat <- mono.dat[,c(1:6,8,9,28,7,11:27,29:31)]
mono.dat$CO2 <- as.factor(mono.dat$CO2)
mono.dat$N <- as.factor(mono.dat$N)
mono.dat$SLA <- as.numeric(mono.dat$SLA)
mono.dat$SeedWt <- as.numeric(mono.dat$SeedWt)
mono.dat$ExpYear <- mono.dat$Year-1997

write.csv(mono.dat, here("data", "monocultureTraitsBiomass.csv"))
