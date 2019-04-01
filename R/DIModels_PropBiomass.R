## Running DI Models but with proportion from that year instead of planted proportion

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

# Tissue N Data 
totdat <- read.csv(here("data", "total_clean.csv"),row.names = 1)
# Monocultures
monoplots <- read.csv(here("data", "monoplots.csv"), row.names = 1)
# Abundance matrix
abund.mat <- read.csv(here("data", "abund_mat.csv"), row.names = 1)

##########################################################################
## dataframe using the actual proportions of the species in the plots
df1 <- merge(totdat[,c(2,6,7)], abund.mat)
df1$Ring <- as.factor(df1$Ring)
df1 <- df1[,-8]
colnames(df1)[9:24] <- paste("P", 1:16, sep = "")
df1 <- df1[,c(1,2,4:24,3)]

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
df1 <- PairwiseInt(df1,7)
# Rename interaction terms to PP1:120
colnames(df1)[25:144] <- paste("PP", 1:120, sep="")
## Compute sum of pairwise interactions
df1$PPsum <- apply(df1[25:144], MARGIN=1, FUN=sum)


df1$l.year <- log(df1$ExpYear)
df1$f.year <- as.factor(df1$ExpYear)
df1$s.year <- df1$ExpYear^2
####################################################################
## Model fitting
## 1. Find best model with time as a main effect
## 2. Find best model with time and environment interactions

#####################################################################

M0 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
            P10 + P11 + P12 + P13 + P14 + P15, random = ~1|Ring, data=df1, method = "ML",
          correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M0)
#################################################################
## M1a: IDENTITY + + AVERAGE PAIRWISE INTERACTION + CO2 + N + CO2:N + LINEAR YEAR. 
#################################################################
M1a <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N + CO2 +
             N:CO2 + ExpYear, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
anova(M1a)						
summary(M1a)
###########################################################
## M1b: IDENTITY + + AVERAGE PAIRWISE INTERACTION + CO2 + N + CO2:N + LOG YEAR. 
###########################################################
M1b <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
anova(M1b)						
summary(M1b)

###########################################################
## M1c: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 + N + CO2:N + FACTOR YEAR. 
###########################################################
M1c <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N + CO2 +
             N:CO2 + f.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
anova(M1c)						
summary(M1c)

###########################################################
## M1d: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N * LINEAR YEAR. 
###########################################################
M1d <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*ExpYear, 
           random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1d)

###########################################################
## M1e: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N * LOG YEAR. 
###########################################################
M1e <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*l.year, 
           random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1e)

###########################################################
## M1f: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N * f.YEAR. 
###########################################################
M1f <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*f.year, 
           random = ~1|Ring, data=df1, method = "ML", 
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1f)

###########################################################
## M1f: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N + QUAD YEAR. 
###########################################################
M1g <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2 + ExpYear*s.year, 
           random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))
summary(M1g)

###########################################################
## M1h: IDENTITY + AVERAGE PAIRWISE INTERACTION + CO2 * N + QUAD YEAR. 
###########################################################
M1h <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N*CO2*ExpYear*s.year, 
           random = ~1|Ring, data=df1, method = "ML",
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
             N:CO2 + ExpYear, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Linear year for interaction and log year for main
M2b <- lme(Nitrogen ~ P1*N*CO2*ExpYear + P2*N*CO2*ExpYear + P3*N*CO2*ExpYear + 
             P4*N*CO2*ExpYear + P5*N*CO2*ExpYear + P6*N*CO2*ExpYear + P7*N*CO2*ExpYear + 
             P8*N*CO2*ExpYear + P9*N*CO2*ExpYear + P10*N*CO2*ExpYear + 
             P11*N*CO2*ExpYear + P12*N*CO2*ExpYear + P13*N*CO2*ExpYear + P14*N*CO2*ExpYear +
             P15*N*CO2*ExpYear + PPsum*N*CO2*ExpYear + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Log year for interaction and linear year main
M2c <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum*N*CO2*l.year + N + CO2 +
             N:CO2 + ExpYear, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Log year for interaction and main
M2d <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum*N*CO2*l.year + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

### Quad year for interaction and main 
M2e <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + P3*N*CO2*ExpYear*s.year + 
             P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year + P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year + P14*N*CO2*ExpYear*s.year +
             P15*N*CO2*ExpYear*s.year + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

## Quad year for interaction log year for main
M2f <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + 
             P3*N*CO2*ExpYear*s.year + P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year +
             P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year +
             P14*N*CO2*ExpYear*s.year + P15*N*CO2*ExpYear*s.year + PPsum*N*CO2*ExpYear*s.year + 
             N + CO2 + N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
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
             N:CO2 + ExpYear, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3b <- lme(Nitrogen ~ P1*N*CO2*ExpYear + P2*N*CO2*ExpYear + P3*N*CO2*ExpYear + 
             P4*N*CO2*ExpYear + P5*N*CO2*ExpYear + P6*N*CO2*ExpYear + P7*N*CO2*ExpYear + 
             P8*N*CO2*ExpYear + P9*N*CO2*ExpYear + P10*N*CO2*ExpYear + 
             P11*N*CO2*ExpYear + P12*N*CO2*ExpYear + P13*N*CO2*ExpYear + P14*N*CO2*ExpYear +
             P15*N*CO2*ExpYear + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3c <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum+ N + CO2 +
             N:CO2 + ExpYear, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3d <- lme(Nitrogen ~ P1*N*CO2*l.year + P2*N*CO2*l.year + P3*N*CO2*l.year + 
             P4*N*CO2*l.year + P5*N*CO2*l.year + P6*N*CO2*l.year + P7*N*CO2*l.year + 
             P8*N*CO2*l.year + P9*N*CO2*l.year + P10*N*CO2*l.year + 
             P11*N*CO2*l.year + P12*N*CO2*l.year + P13*N*CO2*l.year + P14*N*CO2*l.year +
             P15*N*CO2*l.year + PPsum+ N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3e <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + P3*N*CO2*ExpYear*s.year + 
             P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year + P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year + P14*N*CO2*ExpYear*s.year +
             P15*N*CO2*ExpYear*s.year + PPsum + N + CO2 +
             N:CO2 + l.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M3f <- lme(Nitrogen ~ P1*N*CO2*ExpYear*s.year + P2*N*CO2*ExpYear*s.year + P3*N*CO2*ExpYear*s.year + 
             P4*N*CO2*ExpYear*s.year + P5*N*CO2*ExpYear*s.year + P6*N*CO2*ExpYear*s.year + P7*N*CO2*ExpYear*s.year + 
             P8*N*CO2*ExpYear*s.year + P9*N*CO2*ExpYear*s.year + P10*N*CO2*ExpYear*s.year + 
             P11*N*CO2*ExpYear*s.year + P12*N*CO2*ExpYear*s.year + P13*N*CO2*ExpYear*s.year + P14*N*CO2*ExpYear*s.year +
             P15*N*CO2*ExpYear*s.year + PPsum + N + CO2 +
             N:CO2 + ExpYear*s.year, random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

anova(M1f, M3a, M3b, M3c, M3d, M3e, M3f)

#############################################################################
## M4 : Interactions with environment and year on diversity
#############################################################################
M4a <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*ExpYear, 
           random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4b <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*l.year, 
           random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*f.year, 
           random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4d <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
             P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2*ExpYear*s.year, 
           random = ~1|Ring, data=df1, method = "ML",
           correlation = corCAR1(form = ~ 1 | Ring/Plot))

anova(M1f, M4a, M4b, M4c, M4d)
anova(M4c, M1f)

###################################################################
## Moving forward with M4c
##################################################################
M4c1 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2 + PPsum:ExpYear +
              PPsum:N:ExpYear + PPsum:CO2:ExpYear + PPsum:CO2:N:ExpYear + f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c2 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2 + PPsum:l.year +
              PPsum:N:l.year + PPsum:CO2:l.year + PPsum:CO2:N:l.year+ f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c3 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + PPsum*N*CO2 + PPsum:f.year +
              PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c4 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:ExpYear + P2:ExpYear + P3:ExpYear + P4:ExpYear + P5:ExpYear + 
              P6:ExpYear + P7:ExpYear + P8:ExpYear + P9:ExpYear + 
              P10:ExpYear + P11:ExpYear + P12:ExpYear + P13:ExpYear + P14:ExpYear + 
              P15:ExpYear + PPsum*N*CO2 + PPsum:f.year +
              PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c5 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:l.year + P2:l.year + P3:l.year + P4:l.year+ P5:l.year + 
              P6:l.year + P7:l.year + P8:l.year + P9:l.year + 
              P10:l.year + P11:l.year + P12:l.year + P13:l.year + P14:l.year + 
              P15:l.year + PPsum*N*CO2 + PPsum:f.year +
              PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c6 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:ExpYear + P2:ExpYear + P3:ExpYear + P4:ExpYear + P5:ExpYear + 
              P6:ExpYear + P7:ExpYear + P8:ExpYear + P9:ExpYear + 
              P10:ExpYear + P11:ExpYear + P12:ExpYear + P13:ExpYear + P14:ExpYear + 
              P15:ExpYear + P1:s.year + P2:s.year + P3:s.year + P4:s.year + 
              P5:s.year + P6:s.year + P7:s.year + P8:s.year + P9:s.year + 
              P10:s.year + P11:s.year + P12:s.year + P13:s.year + P14:s.year + 
              P15:s.year +  PPsum*N*CO2 + PPsum:f.year +
              PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

M4c7 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:f.year + P2:f.year + P3:f.year + P4:f.year + 
              P6:f.year + P7:f.year + P8:f.year + P9:f.year + 
              P10:f.year + P11:f.year + P12:f.year + P13:f.year + P14:f.year + 
              P15:f.year +  PPsum*N*CO2 + PPsum:f.year +
              PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCAR1(form = ~ 1 | Ring/Plot))

## compare models
anova(M1f, M4c, M4c1, M4c2, M4c3, M4c4, M4c5, M4c6, M4c7)


# M4c7 has the lowest AIC value by 100 AIC points
# This model includes Identity X year factor effects and 
# diveristy x year factor x N x CO2 (no 4-way interaction)

# Comparing different autocorrelation structures
M4c7.1 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
              P10 + P11 + P12 + P13 + P14 + P15 + 
              P1:f.year + P2:f.year + P3:f.year + P4:f.year + 
              P6:f.year + P7:f.year + P8:f.year + P9:f.year + 
              P10:f.year + P11:f.year + P12:f.year + P13:f.year + P14:f.year + 
              P15:f.year +  PPsum*N*CO2 + PPsum:f.year +
              PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
            random = ~1|Ring, data=df1, method = "ML",
            correlation = corCompSymm(form = ~ 1 | Ring/Plot))

M4c7.2 <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
                P10 + P11 + P12 + P13 + P14 + P15 + 
                P1:f.year + P2:f.year + P3:f.year + P4:f.year + 
                P6:f.year + P7:f.year + P8:f.year + P9:f.year + 
                P10:f.year + P11:f.year + P12:f.year + P13:f.year + P14:f.year + 
                P15:f.year +  PPsum*N*CO2 + PPsum:f.year +
                PPsum:N:f.year + PPsum:CO2:f.year + f.year, 
              random = ~1|Ring, data=df1, method = "ML")


anova(M4c7, M4c7.1, M4c7.2)

### Graphs
x <- predictSE.lme(M4c7.1, newdata = df1)
newdat <- cbind(x$fit, df1)
colnames(newdat)[1] <- "fit"
avgpred<- aggregate(newdat[,c(1)], by = list(newdat[,"CO2"], newdat[,"N"], newdat[,"SR"], newdat[,"ExpYear"]), FUN = "mean")
names(avgpred) <- c("CO2", "N", "SR", "Year", "fit")
avgpred1<- aggregate(df1[,c(24)], by = list(df1[,"CO2"], df1[,"N"], df1[,"SR"], df1[,"ExpYear"]), FUN = "mean")
names(avgpred1) <- c("CO2", "N", "SR", "Year", "fit")

cbp1 <- c("#999999","#009E73", "#56B4E9",  "#E69F00")
ggplot(aes(x = Nitrogen, y = x$fit), data = newdat) + 
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
  labs(x = "Experiment Year", y = "Tissue %N", color = "Planted \nRichness") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=14,face="bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"))

ggplot(aes(x = Year, y = fit), data = avgpred) +
  geom_point(aes(color = N:CO2)) +
  geom_line(aes(color = N:CO2)) + 
  facet_wrap(~SR) + 
  geom_point(aes(color = N:CO2, x = Year, y = fit), shape = 2, data = avgpred1) +
  scale_color_manual( values = cbp1) +
  labs(x = "Experiment Year", y = "Tissue %N", color = "Treatment") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=14,face="bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"))

psub <- unique(monoplots$Plot)
monosub <- df1[df1$Plot %in% psub,]
monocul <- merge(monoplots, df1)
monocul$monospecies <- as.character(monocul$monospecies)
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
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=14,face="bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x= "Experiment Year", y = "Tissue %N", color = "Species") 


monocul.pred <- merge(monoplots, newdat)
monocul.pred$monospecies <- factor(monocul.pred$monospecies, 
                                   levels = c("Achillea millefolium", "Asclepias tuberosa",
                                              "Solidago rigida", "Amorpha canescens", 
                                              "Lespedeza capitata", "Lupinus perennis", 
                                              "Petalostemum villosum","Agropyron repens", 
                                              "Bromus inermis", "Koeleria cristata", 
                                              "Poa pratensis", "Andropogon gerardi", 
                                              "Bouteloua gracilis", "Schizachyrium scoparium",
                                              "Sorghastrum nutans"))

ggplot(aes(x = ExpYear, y = fit), data = monocul.pred) +
  geom_line(aes(color = monospecies)) +
  geom_point(aes(x = Year, y = Nitrogen, color = monospecies, alpha = .5), data = monocul)+
  facet_wrap(~N:CO2) +
  scale_color_manual(values = pal) +
  theme_linedraw() +
  ylim(c(0.3,4.4)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=14,face="bold"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x= "Experiment Year", y = "Tissue %N", color = "Species") + 
  guides(alpha = FALSE)


smod <- lme(Nitrogen ~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + 
                    P10 + P11 + P12 + P13 + P14 + P15 + PPsum + N + CO2 + f.year,
                  random = ~1|Ring/Plot, data=df1, method = "ML",
                  correlation = corCAR1(form = ~ 1 | Ring/Plot))

y <- predictSE.lme(smod, df1)
newdat <- cbind(y$fit, df1)
names(newdat)[1] <- "fit"


gr.dat <- data.frame(Bar = NA, TissueN = NA)
Schizplots <- monoplots[monoplots$monospecies == "Schizachyrium scoparium","Plot"]
Schiz.amb <- c("Schizamb", mean(newdat[newdat$Plot %in% Schizplots & newdat$ExpYear == 1 
                                        & newdat$CO2 == "Camb" & newdat$N == "Namb", "fit"]))
Schiz.N <- c("SchizN", mean(newdat[newdat$Plot %in% Schizplots & newdat$ExpYear == 1 
                                    & newdat$CO2 == "Camb" & newdat$N == "Nenrich", "fit"]))
Schiz.CO2 <- c("SchizCO2", mean(newdat[newdat$Plot %in% Schizplots & newdat$ExpYear == 1 
                                        & newdat$CO2 == "Cenrich" & newdat$N == "Namb", "fit"]))
Schiz.CO2N <- c("SchizCO2N", mean(newdat[newdat$Plot %in% Schizplots & newdat$ExpYear == 1 
                                          & newdat$CO2 == "Cenrich" & newdat$N == "Nenrich", "fit"]))
SpeciesEf <- c("SpeciesEF", mean(newdat[!(newdat$Plot %in% Schizplots) & newdat$ExpYear == 1 
                                       & newdat$CO2 == "Camb" & newdat$N == "Namb"
                                       & newdat$SR == 1, "fit"]))
YearEf <- c("YearEF", mean(newdat[newdat$Plot %in% Schizplots & newdat$ExpYear != 1 & 
                                   newdat$CO2 == "Camb" & newdat$N == "Namb", "fit"]))
DivEf <- c("DivEF", mean(newdat[newdat$SR == 16 & newdat$ExpYear == 1 & 
                                 newdat$CO2 == "Camb" & newdat$N == "Namb", "fit"]))
gr.dat <- rbind(gr.dat, Schiz.amb, Schiz.N, Schiz.CO2, Schiz.CO2N, SpeciesEf, YearEf, DivEf) 
                
gr.dat <- gr.dat [-1,]
gr.dat$TissueN <- as.numeric(gr.dat$TissueN)
gr.dat$Bar <- factor(gr.dat$Bar, 
                     levels = c("Schizamb", "SchizN", "SchizCO2",
                                "SchizCO2N", "SpeciesEF", "DivEF", 
                                "YearEF"))
gr.dat$diff <- gr.dat$TissueN - gr.dat$TissueN[1]

ggplot()+
  geom_bar(aes(x=Bar, y = TissueN), data = gr.dat, stat = "identity")+
  labs(x = "", y = "Tissue %N") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+
  geom_bar(aes(x=Bar, y = diff), data = gr.dat[gr.dat$Bar!="Schizamb",], stat = "identity")+
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Difference in Tissue %N") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Average across all years

gr.dat1 <- data.frame(Bar = NA, TissueN = NA)

monoamb <- c("MonoAmb", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                                     newdat$N == "Namb" & newdat$ExpYear ==1, "fit"]))
monoN <- c("MonoN", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                                 newdat$N == "Nenrich"& newdat$ExpYear ==1, "fit"]))
monoCO2 <- c("MonoCO2", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                                     newdat$N == "Namb"& newdat$ExpYear ==1, "fit"]))
monoNCO2 <- c("MonoCO2N", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                                       newdat$N == "Nenrich"& newdat$ExpYear ==1, "fit"]))
polyamb <- c("PolyAmb", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                                     newdat$N == "Namb"& newdat$ExpYear ==1, "fit"]))
year20 <- c("Year20", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                                      newdat$N == "Namb" & newdat$ExpYear ==20, "fit"]))
poly20<- c("Poly20", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                                      newdat$N == "Namb"& newdat$ExpYear ==20, "fit"]))
gr.dat1 <- rbind(gr.dat1, monoamb, monoN, monoCO2, monoNCO2, polyamb,year20, poly20)
gr.dat1 <- gr.dat1 [-1,]
gr.dat1$TissueN <- as.numeric(gr.dat1$TissueN)
gr.dat1$Bar <- factor(gr.dat1$Bar, 
                      levels = c("MonoAmb", "MonoN", "MonoCO2", "MonoCO2N",
                                 "PolyAmb", "Year20", "Poly20"))
ggplot()+
  geom_hline(yintercept = gr.dat1$TissueN[1], linetype = 2, alpha = .5) +
  geom_bar(aes(x=Bar, y = TissueN), data = gr.dat1, stat = "identity")+
  labs(x = "", y = "Tissue %N") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Average for year 1 and year 20

gr.dat2 <- data.frame(Bar = NA, TissueN = NA, se = NA)
se <- function(x) sqrt(var(x)/length(x))

monoamb <- c("MonoAmb", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                                     newdat$N == "Namb" & newdat$ExpYear == 1, "fit"]), 
             se(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                        newdat$N == "Namb", "fit"]))
monoN <- c("MonoN", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                                 newdat$N == "Nenrich", "fit"]),
           se(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                      newdat$N == "Nenrich", "fit"]))
monoCO2 <- c("MonoCO2", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                                     newdat$N == "Namb", "fit"]),
             se(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                        newdat$N == "Namb", "fit"]))
monoNCO2 <- c("MonoCO2N", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                                       newdat$N == "Nenrich", "fit"]),
              se(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                         newdat$N == "Nenrich", "fit"]))
polyamb <- c("PolyAmb", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                                     newdat$N == "Namb", "fit"]),
             se(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                        newdat$N == "Namb", "fit"]))
polyN <- c("PolyN", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                                 newdat$N == "Nenrich", "fit"]),
           se(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                      newdat$N == "Nenrich", "fit"]))
polyCO2 <- c("PolyCO2", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                                     newdat$N == "Namb", "fit"]),
             se(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                        newdat$N == "Namb", "fit"]))      
polyCO2N <- c("PolyCO2N", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                                       newdat$N == "Nenrich", "fit"]),
              se(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                         newdat$N == "Nenrich", "fit"]))
monoamb20 <- c("MonoAmb20", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                                         newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]), 
               se(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                          newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]))
monoN20 <- c("MonoN20", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                                     newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]),
             se(newdat[newdat$SR == 1 & newdat$CO2 == "Camb" & 
                        newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]))
monoCO2.20 <- c("MonoCO2.20", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                                           newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]),
                se(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                           newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]))
monoCO2N20 <- c("MonoCO2N20", mean(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                                           newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]),
                se(newdat[newdat$SR == 1 & newdat$CO2 == "Cenrich" & 
                           newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]))
polyamb20 <- c("PolyAmb20", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                                         newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]),
               se(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                          newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]))
polyN20 <- c("PolyN20", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                                     newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]),
             se(newdat[newdat$SR == 16 & newdat$CO2 == "Camb" & 
                        newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]))
polyCO2.20 <- c("PolyCO2.20", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                                           newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]),
                se(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                           newdat$N == "Namb" & newdat$ExpYear == 20, "fit"]))      
polyCO2N20 <- c("PolyCO2N20", mean(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                                           newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]),
                se(newdat[newdat$SR == 16 & newdat$CO2 == "Cenrich" & 
                           newdat$N == "Nenrich" & newdat$ExpYear == 20, "fit"]))
gr.dat2 <- rbind(gr.dat2, monoamb, monoN, monoCO2, monoNCO2, polyN,
                 polyamb, polyCO2, polyCO2N, monoamb20, monoN20,
                 monoCO2.20, monoCO2N20, polyN20, polyamb20, 
                 polyCO2.20, polyCO2N20)
gr.dat2 <- gr.dat2 [-1,]
gr.dat2$TissueN <- as.numeric(gr.dat2$TissueN)
gr.dat2$se <- as.numeric(gr.dat2$se)
gr.dat2$Bar <- factor(gr.dat2$Bar, 
                      levels = c("MonoAmb", "MonoAmb20", "MonoN", "MonoN20","MonoCO2", 
                                 "MonoCO2.20","MonoCO2N","MonoCO2N20","PolyAmb",
                                 "PolyAmb20", "PolyN", "PolyN20", "PolyCO2", "PolyCO2.20",
                                 "PolyCO2N","PolyCO2N20"))
gr.dat2$col <- c(rep("One",8),rep("Twenty",8))

ggplot()+
  geom_bar(aes(x=Bar, y = TissueN, fill = col), data = gr.dat2, stat = "identity")+
  scale_fill_manual(values = c("red", "blue")) +
  geom_errorbar(aes(x = Bar, ymin = TissueN-se, ymax = TissueN+se, width = .25), data = gr.dat2) +
  labs(x = "", y = "Tissue %N") +
  geom_hline(yintercept = gr.dat2$TissueN[1], linetype = 2) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##################################################
## CALCULATE CWM BASED ON AMBIENT MONOCULTURES ##

CWM.mat <- abund.mat[abund.mat$SR !=1,]
MonoTraits <- totdat[totdat$ExpYear == 1 & totdat$SR ==1 & totdat$N == "Namb" &
                       totdat$CO2 == "Camb", c(2,6)]
MonoTraits <- merge(MonoTraits, monoplots[,c(1,2)])

MonoTraits <- aggregate(MonoTraits$Nitrogen, by = list(MonoTraits$monospecies), FUN = mean)
names(MonoTraits) <- c("monospecies", "TissueN")
rownames(MonoTraits) <- MonoTraits$monospecies
MonoTraits$monospecies <- NULL
CWM.mat <- CWM.mat[CWM.mat$Anemone.cylindrica <.2,]
CWM.mat <- CWM.mat[,-13]
ab.mat <- as.matrix(CWM.mat[,c(9:23)])
CWMs <- ab.mat %*% as.matrix(MonoTraits)
CWM.mat <- cbind(CWMs, CWM.mat)
CWM.mat <- CWM.mat[,c(1:3)]
CWM.mat <- merge(CWM.mat, totdat)

CWM.mat$diff <- CWM.mat$TissueN - CWM.mat$Nitrogen

ggplot(aes(x = Nitrogen, y = TissueN), data = CWM.mat[CWM.mat$ExpYear == 1 | CWM.mat$ExpYear == 20,]) + 
  geom_point(aes(color = N:CO2, shape = as.factor(ExpYear), cex = 1)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(~SR) +
  scale_color_manual(values = cbp1) +
  theme_linedraw() +
  guides(cex = FALSE)
