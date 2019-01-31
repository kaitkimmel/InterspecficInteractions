## Diversity-Interaction Modeling: BioCON Tissue N ##
# Author: Kaitlin Kimmel
# Date: January, 7 2018

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
###########################################################################
# load data

# Data 
# 2017 data
TisN17 <- read.csv(here("data", "TisN17_clean.csv"), row.names = 1)
# dataset with % n and % c from aboveground biomass
TisN <- read.csv(here("data", "TisN_clean.csv"),row.names = 1)
# Experimental design of BioCON
ExpDes <- read.csv(here("data", "BioCONExpDes.csv"))

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
df1$colnums <- seq(24,143, by = 1) #PP1 - PP120 columns

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
df <- TisN17
df <- merge(df, ExpDes, by =c("Plot", "Ring"))
df$Ring <- as.factor(df$Ring)
# Check the SR column = number of planted species
which(df$SR != rowSums(df[,c(7:22)])) # Equals 0 - all good!
# Create mixture column for the unique species combinations
mixtures <- unique.data.frame(df[,c(7:22)])
mixtures$mixture <- seq(1:nrow(mixtures))
df <- inner_join(df, mixtures)

# Get proportion planted
df[c(7:22)] <- df[c(7:22)]/df$SR
# Need to put Tissue N as last column
df <- df[,c(1:5,23,7:22,6)]
# Rename species to P1:16
colnames(df)[7:22] <- paste("P", 1:16, sep = "")

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
df <- PairwiseInt(df,6)
# Rename interaction terms to PP1:120
colnames(df)[24:143] <- paste("PP", 1:120, sep="")

## Compute sum of pairwise interactions
df$PPsum <- apply(df[24:143], MARGIN=1, FUN=sum)

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

# Create dummy variable for ring (not sure why I am doing this...)
for(t in unique(df$Ring)) {
  df[paste("RN",t,sep="")] <- ifelse(df$Ring==t,1,0)}

####################################################################
## Model fitting
#####################################################################
##### MODEL 0: All structure 
#################################################################
## mixture is the unique community compositions in the experiments
M0 <- lm (Nitrogen ~ factor(mixture), data = df)
anova(M0)
summary(M0)

####################################################################################
#### MODEL 1: IDENTITY + N +CO2 + N:CO2. NB CO2 and Ring are aliased. 
#  Ring is omitted here as it enters the random effects model as random and we wish to 
#  use the fixed effects version to feed starting values to the random model later 
#  so we start this pattern here.
##################################################################################

nam1 <- paste("P", 1:15, sep="")		#IDENTITY TERMS
nam2 <- paste("+N+CO2+N:CO2 ", sep="") #Environmental Terms
f1 <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),   paste(nam2)))
M1 <- lm(f1, data=df)
anova(M1)
summary(M1)
AIC(M1)
logLik(M1)
anova(M1, M0)
######################################################################
## M1ran: Modification to allow for random Ring
######################################################################
# nam1 as before	#IDENTITY TERMS
# f1 as before 
M1ran<-lme(f1, data=df, random = ~1|Ring, method = "ML") 
summary(M1ran)
anova(M1ran)
AIC(M1ran)
logLik(M1ran)

##################################################################
## M2a: IDENTITY + PAIRWISE INTERACTION + C02+N+N:CO2
##################################################################
nam3 <- paste("PP", 1:120, sep="")
f2 <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"), paste(nam2, collapse = "+"),
                       paste("+"),paste(nam3, collapse = "+")))
M2a <-lm(f2, data=df)
anova(M2a)
summary(M2a) # NA values for some pairwise interactions. 
logLik(M2a)
anova(M2a, M0)

#################################################################
## M2b: IDENTITY + PAIRWISE INTERACTION + C02+N+N:CO2 + Theta
## NB: Does not work because of NA values estimated in model 2a
#################################################################

df$vN <- as.numeric(df$N)-1
df$vCO2 <- as.numeric(df$CO2)-1
df$vCO2N <- df$vCO2*df$vN
parm1 <- paste("b", 1:15, sep="", paste("*P",1:15,sep=""))
parm2_CO2 <- paste("e1", sep="", paste("*vCO2",sep=""))
parm2_N <- paste("e2", sep="", paste("*vN",sep=""))
parm2_CO2N <- paste("e3", sep="", paste("*vCO2N",sep=""))
parm3 <- paste("d", 1:36, sep="", paste("*PP",1:36,sep=""))		
f2b <- as.formula(paste("Nitrogen ~ ", paste(parm1, collapse= "+"),paste("+"),
                        paste(parm2_CO2),paste("+"),paste(parm2_N),
                        paste("+"), paste(parm2_CO2N),paste("+"),
                        paste(parm3, paste("**theta"), collapse= "+")))
svs <- as.vector(c(coef(M2a),1))
svs <- svs[-1]
names(svs) <- c(paste("b", 1:15, sep=""), paste("e", 1:2, sep=""), 
                paste("d", 1:120, sep=""), "e3", "theta")
svs
M2b <- nls(f2b, data=df, start=svs, algorithm = "port", trace=TRUE)

#################################################################
##  M3a: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION. 
#################################################################

f3 <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),paste(nam2, collapse="+"), 
                       paste("+"),paste("PPsum")))
M3a <- lm(f3, data=df)
anova(M3a)						
summary(M3a)
AIC(M3a)
logLik(M3a)
anova(M3a, M0)
###########################################################
#####   M3ar:  M3a + Ring as random
###########################################################

M3ar <- lme(f3,random=~1|Ring,  data=df)
anova(M3ar)						
summary(M3ar)
AIC(M3ar)
logLik(M3ar)
qqnorm(residuals(M3ar))
plot(M3ar)

#################################################################
##  M3b: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION + THETA 
#################################################################
parm1 <- paste("b", 1:15, sep="", paste("*P",1:15,sep=""))
parm2_CO2 <- paste("e1", sep="", paste("*vCO2",sep=""))
parm2_N <- paste("e2", sep="", paste("*vN",sep=""))
parm2_CO2N <- paste("e3", sep="", paste("*vCO2N",sep=""))
parm3_Dav <- paste("dav", paste("*PP",1:120,sep=""))	
f3b <- as.formula(paste("Nitrogen ~ ",
                        paste(parm1, collapse= "+"),	paste("+"),
                        paste(parm2_CO2),paste("+"),paste(parm2_N),paste("+"),paste(parm2_CO2N),paste("+"),
                        paste(parm3_Dav, paste("**theta"), collapse= "+")))
svs_Dav <- as.vector(c(coef(M3),1))
svs_Dav <- svs_Dav[-1]
names(svs_Dav) <- c(paste("b", 1:15, sep=""),paste("e", 1:2, sep=""),  "dav", "e3", "theta")
svs_Dav
M3b<-nls(f3b, data=df,start=svs_Dav, algorithm = "port", trace=TRUE)
residuals(M3b)
summary(M3b)
AIC(M3b, M3)
logLik(M3b)
logLik(M3)
qqnorm(M3b)
plot(M3b)
res <- summary(M3b)$residuals
RMS <- sqrt(sum(res^2)/310)

#################################################################
## M3c: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION 
## Using Theta estimated from model M3b to create av pairwise 
## interaction variable. 
#################################################################
theta=coef(M3b)[20] # for TPPsum model
df$vN <- as.numeric(df$N)-1
df$vCO2 <- as.numeric(df$CO2)-1
df$vCO2N <- df$vCO2*df$vN
df[163:282]<-df[24:143]^theta
colnames(df)[163:282] <- paste("TPP", 1:120, sep="")
df$TPPsum <- apply(df[163:282], MARGIN=1, FUN=sum)

nam2a <- paste("+vCO2+vN+vCO2N ", sep="")
f3c <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),paste(nam2a, collapse="+"), 
                        paste("+"),("TPPsum")))
M3c <- lm(f3c, data=df)
anova(M3c)						
summary(M3c)
AIC(M3c)
logLik(M3c)
qqnorm(residuals(M3c))
lrtest(M3, M3c) 
anova(M0, M3c)
# model3c is significantly different than M3 - so theta is different from one
#################################################################
## M3cr: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION 
## Using Theta estimated from model M3b to create av pairwise 
## interaction variable and random effect of ring
#################################################################
M3cr <- lme(f3c, random = ~1|Ring, data=df, method = "ML")
anova(M3cr)						
summary(M3cr)
AIC(M3cr)
logLik(M3cr)
qqnorm(residuals(M3cr))
plot(M3cr)



##################################################################
## M4a: IDENTITY + CO2 + N + N:CO2 + FUNCTIONAL GROUPS
##################################################################
f4 <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),
                        paste(nam2, collapse= "+"), paste("+"),
                        paste("PPwfg1 + PPwfg2 + PPwfg3 + PPwfg4 + PPbfg12 + PPbfg13 + PPbfg14 + PPbfg23 + PPbfg24 + PPbfg34")))
M4a <- lm(f4, data=df)
anova(M4a)
AIC(M4a)
logLik(M4a)
summary(M4a)
anova(M4a, M0)
####################################################################
## M4b: IDENTITY + FUNCTIONAL GROUP INTERACTIONS + N + CO2 + N:CO2 + THETA 
####################################################################
parm1 <- paste("b", 1:15, sep="", paste("*P",1:9,sep=""))
parm2_CO2 <- paste("e1", sep="", paste("*vCO2",sep=""))
parm2_N <- paste("e2", sep="", paste("*vN",sep=""))
parm2_CO2N <- paste("e3", sep="", paste("*vCO2N",sep=""))
parm5 <- colnames(df[df1$colnums[which(df1$Within == 1)]]) # within FG1
parm6 <- colnames(df[df1$colnums[which(df1$Within == 2)]]) # within FG2
parm7 <- colnames(df[df1$colnums[which(df1$Within == 3)]]) # within FG3
parm8 <- colnames(df[df1$colnums[which(df1$Within == 4)]]) # within FG4
parm9 <- colnames(df[df1$colnums[which(df1$Between == 12)]]) # between FG1&2
parm10 <- colnames(df[df1$colnums[which(df1$Between == 13)]]) # between FG1&3
parm11 <- colnames(df[df1$colnums[which(df1$Between == 14)]]) # between FG1&4
parm12 <- colnames(df[df1$colnums[which(df1$Between == 23)]]) # between FG2&3
parm13 <- colnames(df[df1$colnums[which(df1$Between == 24)]]) # between FG2&4
parm14 <- colnames(df[df1$colnums[which(df1$Between == 34)]]) # between FG3&4


f4b <- as.formula(paste("Nitrogen ~ ", paste(parm1, collapse= "+"),paste("+"),
                        paste(parm2_CO2, collapse= "+"),paste("+"),
                        paste(parm2_N, collapse= "+"),paste("+"),
                        paste(parm2_CO2N, collapse= "+"),paste("+"),
                        paste(paste("dwfg1*"), parm5, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dwfg2*"), parm6, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dwfg3*"), parm7, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dwfg4*"), parm8, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dbfg12*"), parm9, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dbfg13*"), parm10, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dbfg14*"), parm11, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dbfg23*"), parm12, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dbfg24*"), parm13, paste("**theta"), collapse= "+"),paste("+"),
                        paste(paste("dbfg34*"), parm14, paste("**theta"), collapse= "+")       ))

svs <- as.vector(c(coef(M4),1))[-1]
names(svs) <- c(paste("b", 1:15, sep=""), paste("e", 1:2, sep=""), 
                "dwfg1", "dwfg2", "dwfg3", "dwfg4", 
                "dbfg12", "dbfg13", "dbfg14", "dbfg23", "dbfg24","dbfg34", 
                "e3","theta")
svs
M4b<-nls(f4b, data=df, start=svs, trace=TRUE, algorithm = "port")
summary(M4b)
logLik(M4)
logLik(M4b)

########################################################################
## MODEL 4cr: IDENTITY +vCO2+vN+vCO2*vN + FUNCTIONAL GROUP INTERACTIONS 
# with theta from M3b and Ring random
########################################################################
theta=coef(M3b)[20] # for TPPsum model
df$vN <- as.numeric(df$N)-1
df$vCO2 <- as.numeric(df$CO2)-1
df$vCO2N <- df$vCO2*df$vN

df[163:282]<-df[24:143]^theta
colnames(df)[163:282] <- paste("TPP", 1:120, sep="")
df$TPPsum <- apply(df[163:282], MARGIN=1, FUN=sum)
df1$Tcolnums <- df1$colnums + 139
df$TPPwfg1 <- apply(df[df1$Tcolnums[which(df1$Within == 1)]],MARGIN=1,FUN=sum)
df$TPPwfg2 <- apply(df[df1$Tcolnums[which(df1$Within == 2)]],MARGIN=1,FUN=sum)
df$TPPwfg3 <- apply(df[df1$Tcolnums[which(df1$Within == 3)]],MARGIN=1,FUN=sum)
df$TPPwfg4 <- apply(df[df1$Tcolnums[which(df1$Within == 4)]],MARGIN=1,FUN=sum)
df$TPPbfg12 <- apply(df[df1$Tcolnums[which(df1$Between == 12)]],MARGIN=1,FUN=sum)
df$TPPbfg13 <- apply(df[df1$Tcolnums[which(df1$Between == 13)]],MARGIN=1,FUN=sum)
df$TPPbfg14 <- apply(df[df1$Tcolnums[which(df1$Between == 14)]],MARGIN=1,FUN=sum)
df$TPPbfg23 <- apply(df[df1$Tcolnums[which(df1$Between == 23)]],MARGIN=1,FUN=sum)
df$TPPbfg24 <- apply(df[df1$Tcolnums[which(df1$Between == 24)]],MARGIN=1,FUN=sum)
df$TPPbfg34 <- apply(df[df1$Tcolnums[which(df1$Between == 34)]],MARGIN=1,FUN=sum)

f4c <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),paste(nam2a,"+"),
                  paste("TPPwfg1+ TPPwfg2 + TPPwfg3 + TPPwfg4 + TPPbfg12 + TPPbfg13 +
                        TPPbfg14 + TPPbfg23 + TPPbfg24 + TPPbfg34")))
M4c <- lm(f4c, data = df)
anova(M4c)
anova(M4c, M0)
lrtest(M4c, M4a)
M4cr<-lme(f4c,random=~1|Ring, data=df, method="ML")
summary(M4cr)
logLik(M4cr)
anova(M4cr, M3cr, M1ran)
AICc(M4cr, M3cr, M1ran)
####################################################################
## M4dr: IDENTITY + CO2 +N +CO2:N + FUNCTIONAL GROUP +
## Environmental interactions with identity and diversity
## Using Theta estimated from model M3b to create av pairwise 
## interaction variable and random effect of ring
##################################################################
nam2b <- paste("+vCO2+vN ", sep="")
f4d <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),paste(nam2b,"+"),
                        paste(nam1, collapse = "*vCO2 *vN +"),
                        paste("*vCO2*vN +"),
                        paste("TPPwfg1*vCO2*vN + TPPwfg2*vCO2*vN+ TPPwfg3 *vCO2*vN+ 
                              TPPwfg4*vCO2*vN + TPPbfg12*vCO2*vN + TPPbfg13*vCO2*vN +
                              TPPbfg14*vCO2*vN + TPPbfg23*vCO2*vN + TPPbfg24*vCO2*vN
                              + TPPbfg34*vCO2*vN")))
M4d <- lm(f4d, data = df)
summary(M4d)
M4dr <- lme(f4d, random = ~1|Ring, data = df, method = "ML")
summary(M4dr)


####################################################################
## M3dr: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION +
## Environmental interactions with identity and diversity
## Using Theta estimated from model M3b to create av pairwise 
## interaction variable and random effect of ring
##################################################################
df$vN <- as.numeric(df$N)-1
df$vCO2 <- as.numeric(df$CO2)-1
df$vCO2N <- df$vCO2*df$vN
nam2b <- paste("+vCO2+vN ", sep="")
f3d <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse = "+"), paste(nam2b, collapse="+"), 
                        paste("+"),("TPPsum"), ("+"), paste(nam1, collapse = "*vCO2 *vN +"),
                        paste("*vCO2*vN +"), 
                        paste("TPPsum*vCO2*vN")))
M3d <- lm(f3d, data = df)
summary(M3d)
M3dr <- lme(f3d, random = ~1|Ring, data=df, method = "ML")
summary(M3dr)

#####################################################################
## M3er: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION +
## Environmental interactions with identity and diversity effects 
## No three way interaction
## Using Theta estimated from model M3b to create av pairwise 
## interaction variable and random effect of ring
##################################################################
f3e <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse = "+"), paste(nam2a, collapse="+"), 
                        paste("+"),("TPPsum"), ("+"), paste(nam1, collapse = ":vCO2 +"),
                        paste(":vCO2 +"), paste(nam1, collapse = ":vN +"),
                        paste(":vN +"), paste("TPPsum:vCO2 + "), paste("TPPsum:vN")))
M3e <- lm(f3e, data = df)
summary(M3e)
M3er <- lme(f3e, random = ~1|Ring, data = df, method = "ML")
summary(M3er)

#####################################################################
## M3fr: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION +
## Environmental interactions with identity effects 
## No three way interaction
## Using Theta estimated from model M3b to create av pairwise 
## interaction variable and random effect of ring
##################################################################
f3f <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse = "+"), paste(nam2a, collapse="+"), 
                        paste("+"),("TPPsum"), ("+"), paste(nam1, collapse = ":vCO2 +"),
                        paste(":vCO2 +"), paste(nam1, collapse = ":vN +"), paste(":vN")))
M3f <- lm(f3f, data = df)
summary(M3f)
M3fr <- lme(f3f, random = ~1|Ring, data = df, method = "ML")
summary(M3fr)
anova(M3fr)

#####################################################################
## M3gr: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION +
## Environmental interactions with average diversity effects 
## No three way interaction
## Using Theta estimated from model M3b to create av pairwise 
## interaction variable and random effect of ring
##################################################################
f3g <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse = "+"), paste(nam2a, collapse="+"), 
                        paste("+"),("TPPsum"), ("+"), ("TPPsum:vCO2 +"), ("TPPsum:vN")))
M3g <- lm(f3g, data = df)
summary(M3g)
M3gr <- lme(f3g, random = ~1|Ring, data = df, method = "ML")
summary(M3gr)
anova(M3fr)

anova(M3gr, M3cr, M3fr, M3er)



## Forest's function

pfun <- function(response.columns,trts) {
  par(mfrow=c(1,1), mar=c(2,2,1,1), oma=c(2,2,2,0))
  for(jj in response.columns) {
    j1 <- sapply(seq(0,2,by=0.01),function(theta) {
      ndf$TPPsum <-rowSums(as.data.frame((120*ndf[24:143])^theta/120))
      form <- as.formula(paste(names(ndf)[jj], "~ P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+P11+P12+P13+P14+P15+TPPsum"))
      fit <- lm(form, data=ndf)
      cbind(AIC(fit),theta)})
    jmin <- j1[2,][j1[1,]==min(j1[1,])]
    
    plot(j1[2,],j1[1,],ylab="AIC",xlab="theta", type="l", main=paste(names(ndf)[jj],jmin))
  }
  
  mtext("AIC",2,0,outer=T)
  mtext("theta",1,0.5,outer=T)
  mtext(trts,3,0.5,outer=T)
}
ndf <- df[df$CO2=="Camb" & df$N=="Namb",]
pfun(23,"Amb_Amb")
ndf <- df[df$CO2=="Camb" & df$N=="Nenrich",]
pfun(23,"Amb_Nenr")
ndf <- df[df$CO2=="Cenrich" & df$N=="Namb",]
pfun(23,"Cenrich_Amb")
ndf <- df[df$CO2=="Cenrich" & df$N=="Namb",]
pfun(23,"Cenrich_Nenrich")
