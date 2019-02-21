## Do species interactions modify trait predictions?

## Kaitlin Kimmel
## November 28,2018

## Goal: To determine how much species interactions and environmental conditions impact
# trait predictions. 


#libraries

library(plyr)
library(dplyr)
library(tidyr)
library(nlme)
library(purrr)
library(reshape)
library(here)
library(ggplot2)

# Data 
TisN <- read.csv(here("data", "TisN_clean.csv"), row.names = 1)
Biomass <- read.delim(here("data","Biomass_BioCON.txt"),na.strings=c("","NA"))
names(Biomass) <- c("Sample", "Date", "Plot", "Ring", "CO2", "N", "SR", "FGR",
                 "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Species", "Biomass")
PercCover<- read.delim(here("data","PercCover_BioCON.txt"), na.strings=c("NA",""))
names(PercCover) <- c("Sample", "Season", "Year", "Plot", "Ring", "CO2", "N", "SR", "FGR",
                    "Exp", "monospecies", "monoFunGroup", "WaterTrt", "TempTrt", "Species", "PercCov")
CDRSPDat <- read.csv(here::here("data","CDRSPDat.csv"), na.strings=c("NA",""))




## How does tissue N depend on time, treatment, and species richness? ##
mod <- lme(Nitrogen ~ CO2*N*SR*ExpYear, random = ~1|Ring/Plot, 
           correlation = corCAR1(form = ~ 1 | Ring/Plot), data = TisN,  method = "ML")

mod1 <- lme(Nitrogen ~ CO2*N*SR*ExpYear, random = ~1|Ring/Plot, 
           correlation = corCompSymm(form = ~ 1 | Ring/Plot), data = TisN,  method = "ML")

mod2 <- lme(Nitrogen ~ CO2*N*SR*ExpYear, random = ~1|Ring/Plot,data = TisN,  method = "ML")

anova(mod, mod1, mod2)

mod3 <- lme(Nitrogen ~ CO2 + N + SR + ExpYear+ CO2:N + CO2:SR + CO2:ExpYear +
              N:SR + N:ExpYear + SR:ExpYear + CO2:N:SR + CO2:N:ExpYear +
              CO2:SR:ExpYear + N:SR:ExpYear, random = ~1|Ring/Plot, 
              correlation = corCAR1(form = ~ 1 | Ring/Plot), data = TisN, method = "ML")

mod4 <- lme(Nitrogen ~ CO2 + N + SR + l.year + CO2:N + CO2:SR + CO2:l.year +
              N:SR + N:l.year + SR:l.year + CO2:N:SR + CO2:N:l.year +
              CO2:SR:l.year + N:SR:l.year, random = ~1|Ring/Plot, 
              correlation = corCAR1(form = ~ 1 | Ring/Plot), data = TisN, method = "ML")
anova(mod3, mod4)

mod5 <- lme(Nitrogen ~ CO2 + N + SR + ExpYear+ CO2:N + CO2:SR + CO2:ExpYear +
              N:SR + N:ExpYear + SR:ExpYear, random = ~1|Ring/Plot, 
              correlation = corCAR1(form = ~ 1 | Ring/Plot), data = TisN, method = "ML")

mod6 <- lme(Nitrogen ~ CO2 + N + SR + l.year+ CO2:N + CO2:SR + CO2:l.year +
              N:SR + N:l.year + SR:l.year, random = ~1|Ring/Plot, 
              correlation = corCAR1(form = ~ 1 | Ring/Plot), data = TisN, method = "ML")

mod7 <- lme(Nitrogen ~ CO2 + N + SR + poly(ExpYear,2)+ CO2:N + CO2:SR + CO2:poly(ExpYear,2) +
              N:SR + N:poly(ExpYear,2) + SR:poly(ExpYear,2), random = ~1|Ring/Plot, 
              correlation = corCAR1(form = ~ 1 | Ring/Plot), data = TisN, method = "ML")
mod8 <- lme(Nitrogen ~ CO2 + N + SR + ExpYear + YearSq + CO2:N + 
              CO2:SR + CO2:YearSq + CO2:ExpYear + N:SR + 
              N:YearSq + N:ExpYear + SR:YearSq + SR:ExpYear, random = ~1|Ring/Plot, 
            correlation = corCAR1(form = ~ 1 | Ring/Plot), data = TisN, method = "ML")

anova(mod5, mod6, mod7, mod8) # polynomial fit the best, but would conclude mostly the same with any model
summary(mod7) 
summary(mod5) # Significant negative year:sr interaction, positive CO2:Year interaction
# Significant negative effect of Year and species richness, poistive effect of N addition


# Getting average trait values based on monocultures# 

MonoSp <- TisN[-which(is.na(TisN$monospecies)),]

ggplot(data = MonoSp, aes(x = ExpYear, y = Nitrogen)) +
  geom_line(aes(color = monospecies)) + 
  facet_grid(CO2~N)+
  ylab("Nitrogen")+
  theme_linedraw()

ggplot(data = MonoSp, aes(x = monospecies, y = Nitrogen)) +
  geom_boxplot(aes(color = CO2:N)) + 
  ylab("Nitrogen") +
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = TisN, aes(x=SR, y = Nitrogen)) +
  geom_point(aes(color = N:CO2)) +
  geom_smooth(aes(color = N:CO2), method = "lm")
## DOES VAR CHANGE WITH MEAN? 

# Need to look at monoculture values in order to predict the community weighted mean
# of polyculture. But there are many ways we can group them based on environment & year.

# First we will look at the average values over each treatment and year
# Grouped by species, year, and treatment
MonoSpMeans_ty <- aggregate(MonoSp$Nitrogen, list(ExpYear = MonoSp$ExpYear, N = MonoSp$N, CO2 = MonoSp$CO2, 
                                     sp = MonoSp$monospecies), mean, na.action = na.omit)

ggplot(data = MonoSpMeans_ty, aes(x = ExpYear, y = x)) +
  geom_line(aes(color = sp))+
  facet_grid(N~CO2)+
  ylab("Nitrogen")+
  theme_linedraw()

# Next we will only consider ambient plots - but still look at those changing through time
# Monoculture ambient means for each year
MonoSpMeans_amby <- MonoSpMeans_ty[which(MonoSpMeans_ty$N == "Namb" & MonoSpMeans_ty$CO2 == "Camb"),]

ggplot(data = MonoSpMeans_amby, aes(x = ExpYear, y = x)) +
  geom_line(aes(color = sp))+
  ylab("Nitrogen")+
  theme_linedraw()

# Now we will just consider ambient plots in the first year
#Ambient year 1 means
MonoSpMeans_amb1 <- MonoSpMeans_ty[which(MonoSpMeans_ty$N == 'Namb' & 
                                           MonoSpMeans_ty$CO2 == 'Camb' &
                                           MonoSpMeans_ty$ExpYear == 1),]
asc <- MonoSpMeans_amb1[order(MonoSpMeans_amb1$x),]
asc$sp <- factor(asc$sp, levels = asc$sp[order(asc$x)])
ggplot(data = asc, aes(x = sp, y = x)) +
  geom_point()+
  ylab("Nitrogen")+
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Now we can get an average for each treatment, but compiled over all years.9
# Grouped by species & treatment, regardless of year
MonoSpMeans_t <- aggregate(MonoSp$Nitrogen, list(N = MonoSp$N, CO2 = MonoSp$CO2, 
                                                 sp = MonoSp$monospecies), mean, na.action = na.omit)

# Averages for ambient plots only
MonoSpMeans_amb <- MonoSpMeans_t[which(MonoSpMeans_t$N == 'Namb' & MonoSpMeans_t$CO2 == 'Camb'),]

asc <- MonoSpMeans_amb[order(MonoSpMeans_amb$x),]
asc$sp <- factor(asc$sp, levels = asc$sp[order(asc$x)])
ggplot(data = asc, aes(x = sp, y = x)) +
  geom_point()+
  ylab("Nitrogen")+
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Grouped by species & year, regardless of treatment
MonoSpMeans_y <- aggregate(MonoSp$Nitrogen, list(ExpYear = MonoSp$ExpYear, sp = MonoSp$monospecies), 
                           mean, na.action = na.omit)

ggplot(data = MonoSpMeans_y, aes(x = ExpYear, y = x)) +
  geom_line(aes(color = sp))+
  ylab("Nitrogen")+
  theme_linedraw()

# Grouped by species, regarless of treatment or year
MonoSpMeans <- aggregate(MonoSp$Nitrogen, list(sp = MonoSp$monospecies), mean, na.action = na.omit)

asc <- MonoSpMeans[order(MonoSpMeans$x),]
asc$sp <- factor(asc$sp, levels = asc$sp[order(asc$x)])
ggplot(data = asc, aes(x = sp, y = x)) +
  geom_point()+
  ylab("Nitrogen")+
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




# Proportional abundance of polyculture species
ab.mat <- Biomass[which(Biomass$SR >0),]
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

# Change "Cenriched " to "Cenrich"
ab.mat$CO2<- gsub(pattern = "Cenriched ", replacement = "Cenrich", ab.mat$CO2)
# Get rid of drought plots
ab.mat <- ab.mat[-which(ab.mat$WaterTrt == 'H2Oneg'),]
# Get rid of warmed plots
ab.mat <- ab.mat[-which(ab.mat$TempTrt == 'HTelv'),]

# Get only August data
ab.mat$Date <- as.Date(ab.mat$Date,"%m/%d/%Y")
ab.mat$Year <- as.numeric(format(ab.mat$Dat, format = "%Y"))
ab.mat$Month <- as.numeric(format(ab.mat$Date, format = "%m"))
ab.mat <- ab.mat[-which(ab.mat$Month == 6),]
ab.mat$ExpYear <- ab.mat$Year - 1997

# Get rid of Misc. Lit and Unsorted Biomass
ab.mat <- ab.mat[-which(ab.mat$Species == 'Unsorted Biomass'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Miscellaneous Litter'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Miscellaneous litter'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Real Weeds'),]
ab.mat <- ab.mat[-which(ab.mat$Species == 'Real Weeds '),]
ab.mat <- ab.mat[-which(ab.mat$Species == '16 Species Weeds'),]
ab.mat <- ab.mat[-which(ab.mat$Species == '16 Species Weeds '),]

# Only considering species originally planted 
###### Check if this is ok by looking at perc. biomass of these species
sp <- unique(TisN$monospecies)
# Biomass of planted species 
planted <- ab.mat[which(ab.mat$Species %in% sp),]
planted_bio <- aggregate(planted$Biomass, by = list(Plot = planted$Plot, Year = planted$Year), FUN = sum)
# Get total Biomass for each plot and year
tot_bio <- aggregate(ab.mat$Biomass, by = list(Plot = ab.mat$Plot, Year = ab.mat$Year), FUN = sum)
colnames(tot_bio)[3] <- "TotalBiomass"
x <- join(tot_bio, planted_bio)
x$diff <- (x$x- x$TotalBiomass)/x$TotalBiomass
length(which(x$diff > -0.1))/nrow(tot_bio) ## 98% of plots have 90% of biomass from 16 species in experiment

## Only consider 16 species for rest of analysis ##
ab.mat <- ab.mat[which(ab.mat$Species %in% sp),]
planted_bio <- aggregate(planted$Biomass, by = list(Plot = planted$Plot, Year = planted$Year), FUN = sum)
colnames(planted_bio)[3] <- "TotalBiomass"
ab.mat <- join(ab.mat, planted_bio)
ab.mat <- ab.mat[,c(3:8,11,14,15,17,18,19)]
ab.mat$propbio <- ab.mat$Biomass/ab.mat$TotalBiomass
ab.mat.sub <- ab.mat[,c(1:6,8,11,13)]
# Duplicate rows (21904, 21905)
ab.mat.sub <- ab.mat.sub[-21905,]
ab.mat.wide <- spread(data = ab.mat.sub, key = Species, fill = 0, value = propbio)

# add functional grouping
ab.mat <- join(ab.mat, CDRSPDat[,c(1:2)])
ab.mat[is.na(ab.mat$FunctionalGroup),]$FunctionalGroup <- "F"
FG.mat <- aggregate(ab.mat$propbio, by=list(Plot = ab.mat$Plot, Ring = ab.mat$Ring, 
                                            CO2 = ab.mat$CO2, N =ab.mat$N, 
                                            SR = ab.mat$SR, FGR = ab.mat$FGR, ExpYear = ab.mat$ExpYear,
                                            FunGroup = ab.mat$FunctionalGroup),
                    FUN= sum)
colnames(FG.mat)[length(FG.mat)] <- "propbio"
FG.wide <- spread(data = FG.mat, key = FunGroup, fill = 0, value = propbio)

# CWM using MonoSpMeans

enviro.mat.full <- ab.mat.wide[,c(1:7)]
ab.mat.full <- ab.mat.wide[,-c(1:7)]

rownames(MonoSpMeans) <- MonoSpMeans[,1]
MonoSpMeans[,1]<- NULL
attempt1 <- as.matrix(ab.mat.full)%*%as.matrix(MonoSpMeans)
data.mat <- cbind(enviro.mat.full, attempt1)
colnames(data.mat)[ncol(data.mat)] <- "N_avg"

# CWM using MonoSpMeans_amb
enviro.mat.full <- ab.mat.wide[,c(1:7)]
ab.mat.full <- ab.mat.wide[,-c(1:7)]

rownames(MonoSpMeans_amb) <- MonoSpMeans_amb[,3]
MonoSpMeans_amb[,c(1:3)]<- NULL
attempt1 <- as.matrix(ab.mat.full)%*%as.matrix(MonoSpMeans_amb)
amb.mat <- cbind(enviro.mat.full, attempt1)
colnames(amb.mat)[ncol(amb.mat)] <- "N_amb"

# CWM using MonoSpMeans_t
# Make 4 matrices (one for each treatment)

# Trait matrices
amb.t <- MonoSpMeans_t[which(MonoSpMeans_t$N == "Namb" & MonoSpMeans_t$CO2 == "Camb"),c(3,4)]
rownames(amb.t) <- amb.t[,1]
amb.t[,1]<- NULL
Nadd.t <- MonoSpMeans_t[which(MonoSpMeans_t$N == "Nenrich" & MonoSpMeans_t$CO2 == "Camb"),c(3,4)]
rownames(Nadd.t) <- Nadd.t[,1]
Nadd.t[,1]<- NULL
Cadd.t <- MonoSpMeans_t[which(MonoSpMeans_t$N == "Namb" & MonoSpMeans_t$CO2 == "Cenrich"),c(3,4)]
rownames(Cadd.t) <- Cadd.t[,1]
Cadd.t[,1]<- NULL
En.t <- MonoSpMeans_t[which(MonoSpMeans_t$N == "Nenrich" & MonoSpMeans_t$CO2 == "Cenrich"),c(3,4)]
rownames(En.t) <- En.t[,1]
En.t[,1]<- NULL

# environmental and abundance matrices
amb.ab <- ab.mat.wide[which(ab.mat.wide$N == "Namb" & ab.mat.wide$CO2 == "Camb"),]
env.amb <- amb.ab[,c(1:7)]
amb.ab <- amb.ab[,-c(1:7)]
Nadd.ab <- ab.mat.wide[which(ab.mat.wide$N == "Nenrich" & ab.mat.wide$CO2 == "Camb"),]
env.Nadd <- Nadd.ab[,c(1:7)]
Nadd.ab <- Nadd.ab[,-c(1:7)]
Cadd.ab <- ab.mat.wide[which(ab.mat.wide$N == "Namb" & ab.mat.wide$CO2 == "Cenrich"),]
env.Cadd <- Cadd.ab[,c(1:7)]
Cadd.ab <- Cadd.ab[,-c(1:7)]
En.ab <- ab.mat.wide[which(ab.mat.wide$N == "Nenrich" & ab.mat.wide$CO2 == "Cenrich"),]
env.En <- En.ab[,c(1:7)]
En.ab <- En.ab[,-c(1:7)]

#CMW calcs. 
attempt1 <- as.matrix(amb.ab)%*%as.matrix(amb.t)
amb.means <- cbind(env.amb, attempt1)
attempt1 <- as.matrix(Nadd.ab)%*%as.matrix(Nadd.t)
Nadd.means <- cbind(env.Nadd, attempt1)
attempt1 <- as.matrix(Cadd.ab)%*%as.matrix(Cadd.t)
Cadd.means <- cbind(env.Cadd, attempt1)
attempt1 <- as.matrix(En.ab)%*%as.matrix(En.t)
En.means <- cbind(env.En, attempt1)

Treat.means <- rbind(amb.means, Nadd.means, Cadd.means, En.means)
colnames(Treat.means)[8] <- "N_trt"


# CWM using MonoSpMeans_y
# Make 18 matrices (one for each year!)
years <- unique(MonoSpMeans_y$ExpYear)
trait.mat <- data.frame()
ab.mat <- data.frame()
env.mat <- data.frame()
cwm.mat <- data.frame()
year.means <- data.frame()

for(i in 1:length(years)){
  trait.mat <- MonoSpMeans_y[which(MonoSpMeans_y$ExpYear == years[i]), c(2,3)]
  ab.mat <- ab.mat.wide[which(ab.mat.wide$ExpYear ==years[i]),]
  env.mat <- ab.mat[,c(1:7)]
  ab.mat <- ab.mat[,-c(1:7)]
  tr.sub <- as.character(trait.mat$sp) %in% as.character(colnames(ab.mat))
  trait.mat <- trait.mat[tr.sub,]
  tr.sub1 <- as.character(colnames(ab.mat)) %in% as.character(trait.mat$sp)
  ab.mat <- ab.mat[,tr.sub1]
  rownames(trait.mat) <- trait.mat[,1]
  trait.mat[,1]<- NULL
  attempt1 <- as.matrix(ab.mat)%*%as.matrix(trait.mat)
  cwm.mat <- cbind(env.mat, attempt1)
  year.means <-rbind(year.means, cwm.mat)
}
colnames(year.means)[8]<- "N_year"

# CWM using MonoSpMeans_y
# Make 18 matrices (one for each year!)
years <- unique(MonoSpMeans_amby$ExpYear)
trait.mat <- data.frame()
ab.mat <- data.frame()
env.mat <- data.frame()
cwm.mat <- data.frame()
ambyear.means <- data.frame()

for(i in 1:length(years)){
  trait.mat <- MonoSpMeans_amby[which(MonoSpMeans_amby$ExpYear == years[i]), c(4,5)]
  ab.mat <- ab.mat.wide[which(ab.mat.wide$ExpYear ==years[i]),]
  env.mat <- ab.mat[,c(1:7)]
  ab.mat <- ab.mat[,-c(1:7)]
  tr.sub <- as.character(trait.mat$sp) %in% as.character(colnames(ab.mat))
  trait.mat <- trait.mat[tr.sub,]
  tr.sub1 <- as.character(colnames(ab.mat)) %in% as.character(trait.mat$sp)
  ab.mat <- ab.mat[,tr.sub1]
  rownames(trait.mat) <- trait.mat[,1]
  trait.mat[,1]<- NULL
  attempt1 <- as.matrix(ab.mat)%*%as.matrix(trait.mat)
  cwm.mat <- cbind(env.mat, attempt1)
  ambyear.means <-rbind(ambyear.means, cwm.mat)
}
colnames(ambyear.means)[8]<- "N_ambyear"


# CWM using MonoSpMeans_ty

ab.mat.wide$Treat <- paste(ab.mat.wide$N, ab.mat.wide$CO2, sep = "_")
MonoSpMeans_ty$Treat <- paste(MonoSpMeans_ty$N, MonoSpMeans_ty$CO2, sep = "_")
years <- unique(MonoSpMeans_y$ExpYear)
treatments <- unique(MonoSpMeans_ty$Treat)
tr.mat <- data.frame()
abund.mat <- data.frame()
enviro.mat <- data.frame()
cwms.mat <- data.frame()
treatyear.means <- data.frame()


for(i in 1:length(years)){
  for(j in 1:length(treatments)){
    tr.mat <- MonoSpMeans_ty[which(MonoSpMeans_ty$ExpYear == years[i] & MonoSpMeans_ty$Treat == treatments[j]), c(4,5)]
    abund.mat <- ab.mat.wide[which(ab.mat.wide$ExpYear ==years[i] & ab.mat.wide$Treat == treatments[j]),]
    enviro.mat <- abund.mat[,c(1:7)]
    abund.mat <- abund.mat[,-c(1:7)]
    tr.sub <- as.character(tr.mat$sp) %in% as.character(colnames(abund.mat))
    tr.mat <- tr.mat[tr.sub,]
    tr.sub1 <- as.character(colnames(abund.mat)) %in% as.character(tr.mat$sp)
    abund.mat <- abund.mat[,tr.sub1]
    rownames(tr.mat) <- tr.mat[,1]
    tr.mat[,1]<- NULL
    attempt1 <- as.matrix(abund.mat)%*%as.matrix(tr.mat)
    cwms.mat <- cbind(enviro.mat, attempt1)
    treatyear.means <-rbind(treatyear.means, cwms.mat)
    }
}

ab.mat.wide$Treat <- NULL
colnames(treatyear.means)[8]<-"N_trtyear"

full.data <- list(TisN, data.mat, amb.mat, Treat.means, year.means, ambyear.means, treatyear.means) %>% 
              reduce(join)
# get rid of years without individual species biomass
full.data <- full.data[-which(is.na(full.data$N_amb)),]
# full.data <- full.data[which(is.na(full.data$monospecies)),]
### Biomass * nitrogen content ??? ###

ggplot(data = full.data[full.data$SR != 1,], aes(x = Nitrogen, y = N_amb)) +
  geom_point(aes(color = as.factor(SR)))+
  geom_smooth(method = "lm", aes(group=SR,linetype=factor(SR), color = factor(SR))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data[full.data$SR != 1,], aes(x = Nitrogen, y = N_ambyear)) +
  geom_point(aes(color = as.factor(SR)))+
  geom_smooth(method = "lm", aes(group=SR,linetype=factor(SR), color = factor(SR))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~ExpYear)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data[full.data$SR != 1,], aes(x = Nitrogen, y = N_avg)) +
  geom_point(aes(color = as.factor(SR)))+
  geom_smooth(method = "lm", aes(group=SR,linetype=factor(SR), color = factor(SR))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data, aes(x = Nitrogen, y = N_trt)) +
  geom_point(aes(color = as.factor(SR)))+
  geom_smooth(method = "lm", aes(group=SR, colour=factor(SR),linetype=factor(SR))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()


ggplot(data = full.data[full.data$SR != 1, ], aes(x = Nitrogen, y = N_trtyear)) +
  geom_point(aes(color = as.factor(SR)))+
  geom_smooth(method = "lm", aes(group=SR, colour=factor(SR),linetype=factor(SR))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data, aes(x = Nitrogen, y = N_year)) +
  geom_point(aes(color = as.factor(SR)))+
  geom_smooth(method = "lm", aes(group=SR, colour=factor(SR),linetype=factor(SR))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data, aes(x = Nitrogen, y = N_ambyear)) +
  geom_smooth(method = lm, aes(color = as.factor(CO2:N)))+
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~ExpYear)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data, aes(x = Nitrogen, y = N_trt)) +
  geom_point(aes(color = as.factor(CO2:N)))+
  geom_smooth(method = lm, aes(color = as.factor(CO2:N)))+
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~ExpYear)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()


ggplot(data = full.data, aes(x = Nitrogen, y = N_year)) +
  geom_point(aes(color = as.factor(CO2:N)))+
  geom_smooth(method = "lm", aes(group=SR, colour=factor(SR),linetype=factor(SR))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~ExpYear)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data, aes(x = Nitrogen, y = N_trtyear)) +
  geom_point(aes(color = as.factor(CO2:N)))+
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~ExpYear)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data, aes(x = Nitrogen, y = N_trt)) +
  geom_point(aes(color = as.factor(CO2:N)))+
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~ExpYear)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = full.data, aes(x = Nitrogen, y = N_year)) +
  geom_point(aes(color = as.factor(CO2:N)))+
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~ExpYear)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()


# Functional group effects
# SR = 4
datsub <- subset(full.data, FGR ==1)

ggplot(data = datsub, aes(x = Nitrogen, y = N_avg)) +
  geom_point(aes(color = as.factor(monoFunGroup)))+
  geom_smooth(method = "lm", aes(colour=factor(monoFunGroup))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = datsub, aes(x = Nitrogen, y = N_amb)) +
  geom_point(aes(color = as.factor(monoFunGroup)))+
  geom_smooth(method = "lm", aes(colour=factor(monoFunGroup))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = datsub, aes(x = Nitrogen, y = N_trt)) +
  geom_point(aes(color = as.factor(monoFunGroup)))+
  geom_smooth(method = "lm", aes(colour=factor(monoFunGroup))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()

ggplot(data = datsub, aes(x = Nitrogen, y = N_trtyear)) +
  geom_point(aes(color = as.factor(monoFunGroup)))+
  geom_smooth(method = "lm", aes(colour=factor(monoFunGroup))) +
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(CO2~N)+
  ylab("Predicted")+
  xlab("Actual") +
  theme_linedraw()


full.data$diffN_avg <- full.data$N_avg - full.data$Nitrogen
full.data$diffN_amb <- full.data$N_amb - full.data$Nitrogen
full.data$diffN_trt <- full.data$N_trt - full.data$Nitrogen
full.data$diffN_year <- full.data$N_year - full.data$Nitrogen
full.data$diffN_trtyear <- full.data$N_trtyear - full.data$Nitrogen

long.data <- melt(full.data, id=c('Plot', 'Ring', 'CO2', 'N', 'SR', 'FGR',
                             'Nitrogen', 'ExpYear'), 
             measure.vars = c('diffN_avg', 'diffN_amb', 'diffN_trt', 'diffN_year',
                              'diffN_trtyear'))
long.data$absval <- abs(long.data$value)

ggplot(data = long.data[long.data$SR !=1,], aes(x = variable, y = absval)) +
  geom_boxplot(aes(color = as.factor(FGR)))+
  facet_grid(CO2~N)+
  ylab("Deviation")+
  xlab("Trait") +
  theme_linedraw() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = long.data, aes(x = variable, y = absval)) +
  geom_boxplot(aes(color = as.factor(SR)))+
  ylab("Deviation")+
  xlab("Trait") +
  theme_linedraw() +   
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = long.data[long.data$SR !=1,], aes(x = variable, y = absval)) +
  geom_boxplot(aes(color = N:CO2))+
  ylab("Deviation")+
  xlab("Trait") +
  theme_linedraw() +   
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


## CWM with PercCover --- Just in case.
## Actually decided against this because Bare ground is taken into consideration, so that the perc cover isn't reflective of aboveground biomass
PercCover$Species<- gsub(pattern = "Poa Pratensis", replacement = "Poa pratensis", PercCover$Species)
PercCover$Species<- gsub(pattern = "Poa pratensis ", replacement = "Poa pratensis", PercCover$Species)
PercCover$Species<- gsub(pattern = "schizachyrium scoparium", replacement = "Schizachyrium scoparium", PercCover$Species)
PercCover$Species<- gsub(pattern = "amorpha canescens", replacement = "Amorpha canescens", PercCover$Species)

PercCover <- PercCover[which(PercCover$Season == "August"),]
PercCover <- PercCover[-which(PercCover$WaterTrt == "H2Oneg"),]
PercCover <- PercCover[-which(PercCover$TempTrt == "HTelv"),]

PercCover<- PercCover[which(PercCover$Species %in% sp),]
PercCover$Prop <- PercCover$PercCov/100
PercCover <- PercCover[-c(28525,28529, 28507, 28527),]

PercCov.long <- spread(data = PercCover, key = Species, fill = 0, value = Prop)



#####################################
## Diversity-Interaction Modeling ##
####################################

# References: Kirwan et al 2009, Ecology & Connolly et al 2013, Ecology

# Need to see where NA values are coming from. 
spmat <- list(ab.mat.wide, full.data[,c(3:9,18,25:29)],data.mat, amb.mat, Treat.means, year.means, treatyear.means) %>% 
  reduce(join)
spmat[,(25:29)] <- abs(spmat[,c(25:29)])

fgmat <- list(FG.wide, full.data[,c(3:9,18,25:29)],data.mat, amb.mat, Treat.means, year.means, treatyear.means) %>% 
                reduce(join)
fgmat[,(14:17)] <- abs(fgmat[,c(14:17)])

# Planted proportion
# If biomass >0 then newval = 1/SR
# for (i in 8:23){
#   for (j in 1:nrow(DivIntMod)){
#     if(DivIntMod[j,i] > 0){
#       DivIntMod[j,i] = 1/DivIntMod[j,5]
#     }
#   }
# }

# DivIntMod$PropCheck <- rowSums(DivIntMod[,c(8:23)])
# nrow(DivIntMod[which(DivIntMod$PropCheck == 1),]) #2818... out of 5124... ooh. no. 


## Code below modified from Connolly et al 2013, Ecology

#################################################################
#### DATA MANIPULATION

#### FUNCTION TO CREATE PAIRWISE INTERACTIONS 
#### X IS DATAFRAME; NC IS NUMBER COLUMNS IN DATAFRAME BEFORE SPECIES PROPORTIONS BEGIN
#### BIOMASS IS THE LAST COLUMN
PairwiseIntTerms <- function(x,nc){
  vk <- NULL
  for (k in 1:nrow(x)){
    vi <- NULL
    for (i in (nc+1):(ncol(x)-2)){
      vj <- NULL
      
      for ( j in (i+1):(ncol(x)-1) ) {
        x[k,i]*x[k,j] -> temp
        names(temp)<- paste(colnames(x)[i], colnames(x)[j], sep=".")	
        vj=c(vj, temp)
      }
      vi <- c(vi, vj)
      vi		
    }
    vk <- rbind(vk, vi)
    row.names(vk) <- c(1:k)	
    print(paste("row", k, "done"))
  }
  x <- data.frame(x, vk)
  return(x)
}


# Because the last column is the function, I will pass a new df in for each 
# of the ways I calculated the CWM
# Start with just one year and work up to every year

try1 <- fgmat[,c(1:11,14)]
try1 <- try1[which(try1$ExpYear == 18),]
#### APPLY THE PAIRWISE INTERACTIONS FUNCTION
try1 <- PairwiseIntTerms(try1, 7)
View(try1)

#### COMPUTE SUM OF PAIRWISE INTERACTIONS OVERALL AND WITHIN 
#### AND BETWEEN FUNCTIONAL GROUP

try1$PPsum <- apply(try1[13:18], MARGIN=1, FUN=sum)
parm1 <- paste("b", 1:9, sep="", paste("*P",1:9,sep=""))
parm2 <- paste("a", 1:3, sep="", paste("*BLB",1:3,sep=""))
parm3 <- paste("d", 1:36, sep="", paste("*PP",1:36,sep=""))	
f2b <- as.formula(paste("Biomass ~ ", paste(parm1, collapse= "+"),paste("+"),
                        paste(parm2, collapse= "+"),paste("+"),
                        paste(parm3, paste("**theta"), collapse= "+")))

