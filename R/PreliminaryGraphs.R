setwd("~/Documents/Research/Interspecific_Interactions")
PercN <- read.csv("~/Documents/Research/Interspecific_Interactions/data/PercN.csv")


PercN <- subset(PercN, PercN$Experiment == "M")
library(ggplot2)

ggplot(PercN, aes(x = as.factor(CountOfSpecies), y = Nitrogen...)) + 
  geom_boxplot(aes(color = Nitrogen.Treatment), 
               position = position_dodge(w = .7)) + 
  labs(x = "SR", y = "Nitrogen", color = "NTrt") +
  theme_classic()

ggplot(PercN, aes(x = as.factor(CountOfSpecies), y = Nitrogen...)) + 
  geom_boxplot(aes(color = CO2.Treatment)) + 
  labs(x = "SR", y = "Nitrogen") +
  theme_classic()

ggplot(PercN, aes(x = as.factor(CountOfSpecies), y = Nitrogen...)) + 
  geom_boxplot(aes(color = CO2.Treatment:Nitrogen.Treatment)) + 
  labs(x = "SR", y = "Nitrogen", color = "C:N") +
  theme_classic()
