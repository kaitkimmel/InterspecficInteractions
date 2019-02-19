# Data cleaning: 

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


