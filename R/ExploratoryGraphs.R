## Exploratory graphs

#libraries
library(ggplot2)
library(here)

# load in cleaned data
TisN <- read.csv(here("data", "TisN_clean.csv"))

# Graphs Tissue N throught time by treatment
# Species Richness
ggplot(data = TisN, aes(x = ExpYear, y = Nitrogen)) +
  geom_smooth(method = lm, aes(color = as.factor(SR)))+
  facet_grid(CO2Trt~NTrt)+
  ylab("Nitrogen")+
  theme_linedraw()
# Functional group richness
ggplot(data = TisN, aes(x = ExpYear, y = Nitrogen)) +
  geom_smooth(method = lm, aes(color = as.factor(FGR)))+
  facet_grid(CO2Trt~NTrt)+
  ylab("Nitrogen")+
  theme_linedraw()
# Functional group and species richness
ggplot(data = TisN, aes(x = ExpYear, y = Nitrogen)) +
  geom_smooth(method = lm, aes(color = as.factor(FGR), linetype = as.factor(SR)))+
  facet_grid(CO2Trt~NTrt)+
  ylab("Nitrogen")+
  theme_linedraw()
# Combined
ggplot(data = TisN, aes(x = ExpYear, y = Nitrogen)) +
  geom_smooth(method = lm, aes(color = NTrt:CO2Trt))+
  ylab("Nitrogen")+
  theme_linedraw()
