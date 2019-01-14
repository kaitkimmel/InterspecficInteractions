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
###########################################################################
# load data

# Data 
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
df1$colnums <- seq(24,143, by = 1)

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

# Subset TisN for last year of experiment, plot, ring, SR, %N, CO2, & N treatments
df <- TisN[TisN$ExpYear == max(TisN$ExpYear), c(3:7,16)]
df <- merge(df, ExpDes, by = c("Plot", "Ring"))
df$Ring <- as.factor(df$Ring)
# Check the SR column = number of planted species
which(df$SR != rowSums(df[,c(7:22)])) # Equals 0 - all good!
# Create mixture column for the unique species combinations
mixtures <- unique.data.frame(df[,c(7:22)])
mixtures$mixture <- seq(1:119)
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
sum(test)
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

#### MODEL 1: IDENTITY + N +CO2 + N:CO2. NB CO2 and Ring are aliased. 
#  Ring is omitted here as it enters the random effects model as random and we wish to 
#  use the fixed effects version to feed starting values to the random model later 
#  so we start this pattern here.

nam1 <- paste("P", 1:15, sep="")		#IDENTITY TERMS
nam2 <- paste("+N+CO2+N:CO2 ", sep="")
f1 <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),   paste(nam2)))
M1 <- lm(f1, data=df)
anova(M1)
summary(M1)
AIC(M1)
logLik(M1)

######################################################################
## M1ran: Modification to allow for random Ring
######################################################################
# nam1 as before	#IDENTITY TERMS
# f1 as before 
M1ran<-lme(f1, data=df, random = ~1|Ring) 
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

#################################################################
## M2b: IDENTITY + PAIRWISE INTERACTION + C02+N+N:CO2 + Theta
## NB: Does not work because of NA values estimated in model 2a
#################################################################

df$vN=as.numeric(df$N)-1
df$vCO2=as.numeric(df$CO2)-1
df$vCO2N=df$vCO2*df$vN
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
M2b <- nls(f2b, data=df, start=svs, trace=TRUE)

#################################################################
##  M3: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION. 
#################################################################

f3 <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),paste(nam2, collapse="+"), 
                       paste("+"),("PPsum")))
M3 <- lm(f3, data=df)
anova(M3)						
summary(M3)
AIC(M3)
logLik(M3)

#################################################################
##  M3: IDENTITY + CO2 +N +CO2:N + AVERAGE PAIRWISE INTERACTION + THETA 
#################################################################

##################################################################
## M4: IDENTITY + CO2 + N + N:CO2 + FUNCTIONAL GROUPS
##################################################################
f4 <- as.formula(paste("Nitrogen ~ ", paste(nam1, collapse= "+"),
                        paste(nam2, collapse= "+"), paste("+"),
                        paste("PPwfg1 + PPwfg2 + PPwfg3 + PPbfg12 + PPbfg13 + PPbfg23")))
M4 <- lm(f4, data=df)
anova(M4)
AIC(M4)
logLik(M4)


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
