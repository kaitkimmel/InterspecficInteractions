## Diversity-Interaction Modeling: BioCON Tissue N ##
# Author: Kaitlin Kimmel
# Date: February 20 2018

#######################################################################
## Purpose: To determine the extent to which interspecific interactions and environmental
# context modify predicted community weighted traits for tissue N

######################################################################
## Methods: I will be using Diversity- Interaction modeling 
# (refs: Kirwan et al 2009 Ecology & Connolly et al 2013 Journal of Ecology)
# Some of the below code is borrowed from Connolley et al 2013 & Forest Isbell 

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
library(tidyr)
library(nlme)
library(dplyr)
library(nls2)
library(lmtest)
library(MuMIn)
library(AICcmodavg)
library(ggplot2)
###########################################################################
# load data

# Data 
totdat <- read.csv(here("data", "total_clean.csv"),row.names = 1)
# Experimental design of BioCON
ExpDes <- read.csv(here("data", "BioCONExpDes.csv"))
# Monocultures
monoplots <- read.csv(here("data", "monoplots.csv"), row.names = 1)

###########################################################################
# Create a datafame with the functional groups to be used later
fungroup <- data.frame(Species = as.character(colnames(ExpDes)[3:18]))
fungroup$fungroup <- c("F", "C3", "L", "C4", "F", "F", "C4", "C3", "C3", 
                       "L", "L", "L", "C3", "C4", "F", "C4")
df1 <- data.frame(Sp1= character(), Sp2 = character(), Fg1 = character(), Fg2 = character())
for(i in 1:(nrow(fungroup)-1)){
  Sp1 = as.character(rep(fungroup$Species[i],16-i))
  Sp2 = as.character(fungroup$Species[c((i+1):nrow(fungroup))])
  Fg1 = as.character(rep(fungroup$fungroup[i], 16-i))
  Fg2 = as.character(fungroup$fungroup[c((i+1):nrow(fungroup))])
  x <- cbind(Sp1, Sp2, Fg1, Fg2)
  df1 <- rbind(df1, x)
}
df1$colnums <- seq(26,145, by = 1) #PP1 - PP120 columns

for(i in 1:nrow(df1)){
  if(df1$Fg1[i] == df1$Fg2[i]){
    df1$Between[i] = 0
    if(df1$Fg1[i] == "F"){
      df1$Within[i] = 1
    }
    if(df1$Fg1[i] == "L"){
      df1$Within[i] = 2
    }
    if(df1$Fg1[i] == "C3"){
      df1$Within[i] = 3
    }
    if(df1$Fg1[i] == "C4"){
      df1$Within[i] = 4
    }
  }
  else{
    df1$Within[i] = 0
    if(df1$Fg1[i] == "F" & df1$Fg2[i] == "L" | 
       df1$Fg1[i] == "L" & df1$Fg2[i] == "F"){
      df1$Between[i] = 12
    }
    if(df1$Fg1[i] == "F" & df1$Fg2[i] == "C3" | 
       df1$Fg1[i] == "C3" & df1$Fg2[i] == "F"){
      df1$Between[i] = 13
    }
    if(df1$Fg1[i] == "F" & df1$Fg2[i] == "C4" | 
       df1$Fg1[i] == "C4" & df1$Fg2[i] == "F"){
      df1$Between[i] = 14
    }
    if(df1$Fg1[i] == "L" & df1$Fg2[i] == "C3" | 
       df1$Fg1[i] == "C3" & df1$Fg2[i] == "L"){
      df1$Between[i] = 23
    }
    if(df1$Fg1[i] == "L" & df1$Fg2[i] == "C4" | 
       df1$Fg1[i] == "C4" & df1$Fg2[i] == "L"){
      df1$Between[i] = 24
    }
    if(df1$Fg1[i] == "C3" & df1$Fg2[i] == "C4" | 
       df1$Fg1[i] == "C4" & df1$Fg2[i] == "C3"){
      df1$Between[i] = 34
    }
  }
}



##########################################################################
# For the species identity model: I will need to create a dataframe that 
# includes the plot, ring, CO2 Treatment, N treatment, the proportions 
# of the species in the plot, and the measured tissue N
########################################################################
# Subset TisN17 for plot, ring, SR, %N, CO2, & N treatments
#df <- TisN[TisN$ExpYear == max(TisN$ExpYear), c(4:7,16)]
df <- totdat
df <- merge(df, ExpDes, by =c("Plot", "Ring"))
df$Ring <- as.factor(df$Ring)
# Check the SR column = number of planted species
which(df$SR != rowSums(df[,c(9:24)])) # Equals 0 - all good!
# Create mixture column for the unique species combinations
mixtures <- unique.data.frame(df[,c(9:24)])
mixtures$mixture <- seq(1:nrow(mixtures))
df <- inner_join(df, mixtures)

# Get proportion planted
df[c(9:24)] <- df[c(9:24)]/df$SR
# Need to put Tissue N as last column
df <- df[,c(1:5,7,8,25,9:24,6)]
# Rename species to P1:16
colnames(df)[9:24] <- paste("P", 1:16, sep = "")

# Below code is adapted from Connolly et al 2013 and Forest Isbell

#################################################################
# Data manipulation
# Function to create pairwise interactions
# X is the data; NC is the number of columns before species proportions
# Ecosystem function is the last column (here, tissue N)
#################################################################

PairwiseInt <- function(x,nc){
  vk <- NULL
  for (k in 1:nrow(x)){
    vi <- NULL
    for (i in (nc+1):(ncol(x)-2)){
      vj <- NULL
      
      for (j in (i+1):(ncol(x)-1)) {
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

#############################################################
# Apply pairwise interaction function - assuming theta = 1
df <- PairwiseInt(df,8)
# Rename interaction terms to PP1:120
colnames(df)[26:145] <- paste("PP", 1:120, sep="")

## Compute sum of pairwise interactions
df$PPsum <- apply(df[26:145], MARGIN=1, FUN=sum)

#########################################################
## Compute within and between functional group interaction sums
# FG1 = Forbs; FG2 = Legumes; FG3= C3; FG4 = C4
df$PPwfg1 <- apply(df[df1$colnums[which(df1$Within == 1)]],MARGIN=1,FUN=sum)
df$PPwfg2 <- apply(df[df1$colnums[which(df1$Within == 2)]],MARGIN=1,FUN=sum)
df$PPwfg3 <- apply(df[df1$colnums[which(df1$Within == 3)]],MARGIN=1,FUN=sum)
df$PPwfg4 <- apply(df[df1$colnums[which(df1$Within == 4)]],MARGIN=1,FUN=sum)
df$PPbfg12 <- apply(df[df1$colnums[which(df1$Between == 12)]],MARGIN=1,FUN=sum)
df$PPbfg13 <- apply(df[df1$colnums[which(df1$Between == 13)]],MARGIN=1,FUN=sum)
df$PPbfg14 <- apply(df[df1$colnums[which(df1$Between == 14)]],MARGIN=1,FUN=sum)
df$PPbfg23 <- apply(df[df1$colnums[which(df1$Between == 23)]],MARGIN=1,FUN=sum)
df$PPbfg24 <- apply(df[df1$colnums[which(df1$Between == 24)]],MARGIN=1,FUN=sum)
df$PPbfg34 <- apply(df[df1$colnums[which(df1$Between == 34)]],MARGIN=1,FUN=sum)

# Test that PP variables sum to df$PPsum to test arithmetic
test= df$PPsum-df$PPwfg1 -df$PPwfg2 -df$PPwfg3 -df$PPwfg4-df$PPbfg12-df$PPbfg13 -df$PPbfg14 -df$PPbfg23 -df$PPbfg24 -df$PPbfg34 
mean(test)
max(df$PPsum)

df$l.year <- log(df$ExpYear)
df$f.year <- as.factor(df$ExpYear)
df$s.year <- df$ExpYear^2
####################################################################
## Model fitting
## 1. Find best model with time as a main effect
## 2. Find best model with time and environment interactions

#####################################################################

#################################################################
## M1a: IDENTITY + + AVERAGE PAIRWISE INTERACTION + CO2 + N + CO2:N + LINEAR YEAR. 
#################################################################
M1a <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
            P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N + CO2 +
            N:CO2 + ExpYear, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
anova(M1a)						
summary(M1a)
###########################################################
## M1b: IDENTITY + + AVERAGE PAIRWISE INTERACTION + CO2 + N + CO2:N + LOG YEAR. 
###########################################################
M1b <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
anova(M1b)						
summary(M1b)

###########################################################
## M1c: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 + N + CO2:N + FACTOR YEAR. 
###########################################################
M1c <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N + CO2 +
             N:CO2 + f.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
anova(M1c)						
summary(M1c)

###########################################################
## M1d: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N * LINEAR YEAR. 
###########################################################
M1d <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*ExpYear, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1d)

###########################################################
## M1e: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N * LOG YEAR. 
###########################################################
M1e <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*l.year, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1e)

###########################################################
## M1f: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N * LOG YEAR. 
###########################################################
M1f <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*f.year, 
           random = ~1|Ring, data=df, method = "ML", 
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1f)
###########################################################
## M1f: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N + QUAD YEAR. 
###########################################################
M1g <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2 + ExpYear*s.year, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1g)

###########################################################
## M1h: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N + QUAD YEAR. 
###########################################################
M1h <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*ExpYear*s.year, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1h)
#############################################################
# Compare M1a,b,c,d,e,f
#############################################################
anova(M1a, M1b, M1c, M1d, M1e, M1f, M1g, M1h)
# M1f has the lowest AIC value and the greatest log likelihood 

############################################################
## M2a: Interactions with environment and year on identity and diversity
############################################################
### Linear Year for interaction and main
M2a <- lme(Nitrogen ~ P1*N*CO2*ExpYear + P2*N*CO2*ExpYear + P3*N*CO2*ExpYear + 
             P4*N*CO2*ExpYear + P5*N*CO2*ExpYear + P6*N*CO2*ExpYear + P7*N*CO2*ExpYear + 
             P8*N*CO2*ExpYear + P9*N*CO2*ExpYear + P10*N*CO2*ExpYear + 
             P11*N*CO2*ExpYear + P12*N*CO2*ExpYear + P13*N*CO2*ExpYear + P14*N*CO2*ExpYear +
             P15*N*CO2*ExpYear + PPsum*N*CO2*ExpYear + N + CO2 +
             N:CO2 + ExpYear, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Linear year for interaction and log year for main
M2b <- lme(Nitrogen ~ P1*N*CO2*ExpYear + P2*N*CO2*ExpYear + P3*N*CO2*ExpYear + 
             P4*N*CO2*ExpYear + P5*N*CO2*ExpYear + P6*N*CO2*ExpYear + P7*N*CO2*ExpYear + 
             P8*N*CO2*ExpYear + P9*N*CO2*ExpYear + P10*N*CO2*ExpYear + 
             P11*N*CO2*ExpYear + P12*N*CO2*ExpYear + P13*N*CO2*ExpYear + P14*N*CO2*ExpYear +
             P15*N*CO2*ExpYear + PPsum*N*CO2*ExpYear + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Log year for interaction and linear year main
M2c <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum*N*CO2*l.year + N + CO2 +
             N:CO2 + ExpYear, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Log year for interaction and main
M2d <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum*N*CO2*l.year + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Quad year for interaction and main 
M2e <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + P3*N*CO2*ExpYear*s.year + 
             P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year + P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year + P14*N*CO2*ExpYear*s.year +
             P15*N*CO2*ExpYear*s.year + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

## Quad year for interaction log year for main
M2f <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + 
             P3*N*CO2*ExpYear*s.year + P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year +
             P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year +
             P14*N*CO2*ExpYear*s.year + P15*N*CO2*ExpYear*s.year + PPsum*N*CO2*ExpYear*s.year + 
             N + CO2 + N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

anova(M1f, M2a, M2b, M2c, M2d, M2e, M2f)

############################################################
## M3: Interactions with environment and year on identity
############################################################
M3a <- lme(Nitrogen ~ P1*N*CO2*ExpYear + P2*N*CO2*ExpYear + P3*N*CO2*ExpYear + 
             P4*N*CO2*ExpYear + P5*N*CO2*ExpYear + P6*N*CO2*ExpYear + P7*N*CO2*ExpYear + 
             P8*N*CO2*ExpYear + P9*N*CO2*ExpYear + P10*N*CO2*ExpYear + 
             P11*N*CO2*ExpYear + P12*N*CO2*ExpYear + P13*N*CO2*ExpYear + P14*N*CO2*ExpYear +
             P15*N*CO2*ExpYear + PPsum + N + CO2 +
             N:CO2 + ExpYear, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3b <- lme(Nitrogen ~ P1*N*CO2*ExpYear + P2*N*CO2*ExpYear + P3*N*CO2*ExpYear + 
             P4*N*CO2*ExpYear + P5*N*CO2*ExpYear + P6*N*CO2*ExpYear + P7*N*CO2*ExpYear + 
             P8*N*CO2*ExpYear + P9*N*CO2*ExpYear + P10*N*CO2*ExpYear + 
             P11*N*CO2*ExpYear + P12*N*CO2*ExpYear + P13*N*CO2*ExpYear + P14*N*CO2*ExpYear +
             P15*N*CO2*ExpYear + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3c <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum+ N + CO2 +
             N:CO2 + ExpYear, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3d <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum+ N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3e <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + P3*N*CO2*ExpYear*s.year + 
             P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year + P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year + P14*N*CO2*ExpYear*s.year +
             P15*N*CO2*ExpYear*s.year + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3f <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + P3*N*CO2*ExpYear*s.year + 
             P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year + P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year + P14*N*CO2*ExpYear*s.year +
             P15*N*CO2*ExpYear*s.year + PPsum + N + CO2 +
             N:CO2 + ExpYear*s.year, random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))


anova(M1f, M3a, M3b, M3c, M3d, M3e, M3f)

#############################################################################
## M4 : Interactions with environment and year on diversity
#############################################################################
M4a <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*ExpYear, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4b <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*l.year, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*f.year, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4d <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*ExpYear*s.year, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

anova(M1f, M4a, M4b, M4c, M4d)
anova(M4c, M1f)

###################################################################
## Moving forward with M4c
##################################################################
M4c1 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2 + PPsum:ExpYear +
            PPsum:N:ExpYear + PPsum:CO2:ExpYear + PPsum:CO2:N:ExpYear + f.year, 
           random = ~1|Ring, data=df, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c2 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2 + PPsum:l.year +
              PPsum:N:l.year + PPsum:CO2:l.year + PPsum:CO2:N:l.year+ f.year, 
            random = ~1|Ring, data=df, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c3 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
            P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2 + PPsum:f.year +
            PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c4 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:f.year + P2:f.year + P3:f.year + P4:f.year + P5:f.year + 
              P6:f.year + P7:f.year + P8:f.year + P9:f.year + 
              P10:f.year + P11:f.year + P12:f.year + P13:f.year + P14:f.year + 
              P15:f.year + PPsum*N*CO2 + PPsum:f.year +
              PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c5 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:f.year + P2:f.year + P3:f.year + P4:f.year + P5:f.year + 
              P6:f.year + P7:f.year + P8:f.year + P9:f.year + 
              P10:f.year + P11:f.year + P12:f.year + P13:f.year + P14:f.year + 
              P15:f.year + PPsum*N*CO2 + PPsum:l.year +
              PPsum:N:l.year + PPsum:CO2:l.year + f.year, 
            random = ~1|Ring, data=df, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c6 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:f.year + P2:f.year + P3:f.year + P4:f.year + P5:f.year + 
              P6:f.year + P7:f.year + P8:f.year + P9:f.year + 
              P10:f.year + P11:f.year + P12:f.year + P13:f.year + P14:f.year + 
              P15:f.year + PPsum*N*CO2 + PPsum:ExpYear +
              PPsum:N:ExpYear + PPsum:CO2:ExpYear + f.year, 
            random = ~1|Ring, data=df, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c7 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:N*CO2 + P2:N*CO2 + P3:N*CO2 + P4:N*CO2 + P5:N*CO2 + 
              P6:N*CO2 + P7:N*CO2 + P8:N*CO2 + P9:N*CO2 + 
              P10:N*CO2 + P11:N*CO2 + P12:N*CO2 + P13:N*CO2 + P14:N*CO2 + 
              P15:N*CO2 + PPsum*N*CO2 + PPsum:ExpYear +
              PPsum:N:ExpYear + PPsum:CO2:ExpYear + f.year, 
            random = ~1|Ring, data=df, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c8 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:f.year + P2:f.year + P3:f.year + P4:f.year + P5:f.year + 
              P6:f.year + P7:f.year + P8:f.year + P9:f.year + 
              P10:f.year + P11:f.year + P12:f.year + P13:f.year + P14:f.year + 
              P15:f.year + PPsum*N*CO2 + PPsum:ExpYear*s.year +
              PPsum:N:ExpYear*s.year + PPsum:CO2:ExpYear*s.year + ExpYear*s.year, 
            random = ~1|Ring, data=df, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

## compare models
anova(M4c, M4c1, M4c2, M4c3, M4c4, M4c5, M4c6, M4c7, M4c8)


# M4c4 has the lowest AIC value by 100 AIC points
# This model includes Identity X year factor effects and 
# diveristy x year factor x N x CO2 (no 4-way interaction)

## Graphs

x <- predictSE.lme(M4c4, newdata = df)
newdf <- cbind(x$fit, df)
colnames(newdf)[1] <- "fit"
avgpred<- aggregate(newdf[,c(1)], by = list(newdf[,"CO2"], newdf[,"N"], newdf[,"SR"], newdf[,"ExpYear"]), FUN = "mean")
names(avgpred) <- c("CO2", "N", "SR", "Year", "fit")
avgpred1<- aggregate(df[,c(25)], by = list(newdf[,"CO2"], newdf[,"N"], newdf[,"SR"], newdf[,"ExpYear"]), FUN = "mean")
names(avgpred1) <- c("CO2", "N", "SR", "Year", "fit")

cbp1 <- c("#999999","#009E73", "#56B4E9",  "#E69F00")
ggplot(aes(x = Nitrogen, y = x$fit), data = newdf) + 
  geom_point(aes(color = N:CO2)) +
  geom_abline(aes(slope = 1, intercept = 0)) + 
  scale_color_manual( values = cbp1) + 
  labs(x = "Actual", y = "Predicted", color = "Treatment") + 
  theme_classic()

ggplot(aes(x = Year, y = fit), data = avgpred) +
  geom_point(aes(color = factor(SR))) +
  geom_line(aes(color = factor(SR))) + 
  facet_wrap(~N:CO2) + 
  geom_point(aes(color = factor(SR), x = Year, y = fit), shape = 2, data = avgpred1) +
  scale_color_manual( values = cbp1) +
  labs(x = "Experiment Year", y = "Tissue %N", color = "Richness") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
ggplot(aes(x = Year, y = fit), data = avgpred) +
  geom_point(aes(color = N:CO2)) +
  geom_line(aes(color = N:CO2)) + 
  facet_wrap(~SR) + 
  geom_point(aes(color = N:CO2, x = Year, y = fit), shape = 2, data = avgpred1) +
  scale_color_manual( values = cbp1) +
  labs(x = "Experiment Year", y = "Tissue %N", color = "Richness") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


monocul <- merge(monoplots, df)
monocul<- aggregate(monocul[,"Nitrogen"], 
                    by = list(monocul[,"CO2"], monocul[,"N"], 
                              monocul[,"monospecies"], monocul[,"ExpYear"]), 
                    FUN = "mean")
names(monocul) <- c("CO2", "N", "monospecies", "Year", "Nitrogen")
monocul <- merge(monocul, unique(monoplots[,c(1,3)]))
monocul$monospecies <- factor(monocul$monospecies, 
                              levels = c("Achillea millefolium", "Asclepias tuberosa",
                                         "Solidago rigida", "Amorpha canescens", 
                                         "Lespedeza capitata", "Lupinus perennis", 
                                         "Petalostemum villosum","Agropyron repens", 
                                         "Bromus inermis", "Koeleria cristata", 
                                         "Poa pratensis", "Andropogon gerardi", 
                                         "Bouteloua gracilis", "Schizachyrium scoparium",
                                          "Sorghastrum nutans"))
pal <- c("#0714C9","#0083FF","#00C0DC",#Forb
         "#8C009F","#BB0076","#F47AFF","#4E0241",#legume
         "#DB5700","#971A00","#FFCA79","#FF2400",#C3
         "#ABC31A","#43BB00","#005000","#00BE66")#C4

ggplot(aes(x = Year, y = Nitrogen), data = monocul) +
  geom_line(aes(color = monospecies)) +
  facet_wrap(~N:CO2) + 
  scale_color_manual(values = pal) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x= "ExpYear", y = "Tissue %N", color = "Species")



ggplot(aes(x = Year, y = fit), data = monocul.pred) +
  geom_line(aes(color = monospecies)) +
  facet_wrap(~N:CO2) + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x= "ExpYear", y = "Tissue %N", color = "Species")

