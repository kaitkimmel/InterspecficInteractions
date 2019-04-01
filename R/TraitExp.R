library(here)
library(ggplot2)
library(ggpubr)
library(dplyr)

# load data
mono.dat <- read.csv(here("data", "monocultureTraitsBiomass.csv"), row.names = 1)
mono.dat$SeedWt <- as.numeric(mono.dat$SeedWt)
# ambient plots only to start
amb.dat <- mono.dat[mono.dat$N == "Namb" & mono.dat$CO2 == "Camb",]

cor.mat <- cor(amb.dat[,c(10:28)], use = "pairwise.complete.obs")

# Plots of trait value v biomass
pdf(file = here("Figures","Traitplots.pdf"), onefile = TRUE)
for(i in 10:27){
  plot.list <- list()
  xval<- names(amb.dat)[i]
  for(j in unique(amb.dat$ExpYear))local({
    j = j
    dat <- amb.dat[amb.dat$ExpYear == j,]
    p1 <- ggplot(aes(x = dat[,i], y = TotalBiomass), data = dat) + 
      geom_point() + 
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = xval, y = "Biomass", title = j) +
      theme_linedraw()
    plot.list[[j]] <<- p1
    })
    
    p2 <- ggplot(aes(x = amb.dat[,i], y = TotalBiomass), data = amb.dat) + 
      geom_point() + 
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = xval, y = "Biomass", title = "All Years") +
      theme_linedraw()
    plot.list[[j+1]] <- p2
    
    print(ggarrange(plotlist = plot.list))
    
}

dev.off()

avg.dat <- amb.dat[,c(3,10:28)] %>%
  group_by(monospecies) %>% 
  summarise_all(funs(mean(., na.rm = TRUE)))
avg.dat <- as.data.frame(avg.dat)


plot.list2 <- list()
for(i in 2:19)local({
  i = i
  xval<- names(avg.dat)[i]
    p1 <- ggplot(aes(x = avg.dat[,i], y = TotalBiomass), data = avg.dat) + 
      geom_point() + 
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = xval, y = "Biomass") +
      theme_linedraw()
    plot.list2[[i-1]] <<- p1
  })

pdf(file = here("Figures","AvgTraitPlots.pdf"), onefile = TRUE, width = 15, height = 10)
print(ggarrange(plotlist = plot.list2))
dev.off()

