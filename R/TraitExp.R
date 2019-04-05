library(here)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(MuMIn)

# set.seed
set.seed(4234)
# load data
mono.dat <- read.csv(here("data", "monocultureTraitsBiomass.csv"), row.names = 1)
mono.dat$SeedWt <- as.numeric(mono.dat$SeedWt)
CDRsp <- read.csv(here("data", "CDRSPDat.csv"))

# ambient plots only to start
amb.dat <- mono.dat[mono.dat$N == "Namb" & mono.dat$CO2 == "Camb",]
amb.dat <- amb.dat[-which(is.na(amb.dat$ExpYear)),]

cor.mat <- cor(amb.dat[,c(10:30)], use = "pairwise.complete.obs")

# Plots of trait value v biomass
pdf(file = here("Figures","Traitplots.pdf"), onefile = TRUE, width = 15, height = 12)

for(i in 10:29){
  plot.list <- list()
  xval<- names(amb.dat)[i]
  for(j in 1:20)local({
    j = j
    dat <- amb.dat[amb.dat$ExpYear == j,]
    p1 <- ggplot(aes(x = dat[,i], y = TotalBiomass), data = dat) + 
      geom_point(aes(color = monospecies)) + 
      geom_smooth(se = FALSE, method = "lm") +
      labs(x = xval, y = "Biomass", title = j) +
      theme_linedraw()
    plot.list[[j]] <<- p1
    })
    
    p2 <- ggplot(aes(x = amb.dat[,i], y = TotalBiomass), data = amb.dat) + 
      geom_point(aes(color = monospecies)) + 
      geom_smooth(se = FALSE, method = "lm") +
      labs(x = xval, y = "Biomass", title = "All Years") +
      theme_linedraw()
    plot.list[[j+1]] <- p2
    
    print(ggarrange(plotlist = plot.list, common.legend = TRUE))
    
}

dev.off()

avg.dat <- amb.dat[,c(1,10:30)] %>%
  group_by(monospecies) %>% 
  summarise_all(funs(mean(., na.rm = TRUE)))
avg.dat <- as.data.frame(avg.dat)


plot.list2 <- list()
for(i in 2:19)local({
  i = i
  xval<- names(avg.dat)[i]
    p1 <- ggplot(aes(x = avg.dat[,i], y = TotalBiomass), data = avg.dat) + 
      geom_point(aes(color = monospecies)) + 
      geom_smooth(se = FALSE, method = "lm") +
      labs(x = xval, y = "Biomass") +
      theme_linedraw()
    plot.list2[[i-1]] <<- p1
  })

pdf(file = here("Figures","AvgTraitPlot.pdf"), onefile = TRUE, width = 15, height = 12)
print(ggarrange(plotlist = plot.list2, common.legend = TRUE))
dev.off()

## Stats
df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("trait", "intercept", "coef", "coef.pval", "RSE", "R2"))
for(i in 10:(length(amb.dat)-1)){
  mod <- lm(TotalBiomass ~ amb.dat[,i], data = amb.dat)
  df[i-9,]<- c(names(amb.dat)[i],summary(mod)$coefficients[1,1], 
               summary(mod)$coefficients[2,1], 
               summary(mod)$coefficients[2,4],
               sqrt(deviance(mod)/df.residual(mod)),
               summary(mod)$r.squared)

}
df$intercept <- as.numeric(df$intercept)
df$coef <- as.numeric(df$coef)
df$coef.pval <- as.numeric(df$coef.pval)
df$RSE <- as.numeric(df$RSE)
df$R2 <- as.numeric(df$R2)

## PCA of traits
pca.mat <- avg.dat
rownames(pca.mat) <- pca.mat[,1]
pca.mat  <- pca.mat[,-c(1,21)]
for(i in 1:length(pca.mat)){
  for(j in 1:nrow(pca.mat)){
    if(is.nan(pca.mat[j,i])){
      pca.mat[j,i] = NA
    }
  }
}

sp.pca <- princomp(pca.mat[,c(1,8:19)], cor = TRUE, scores = TRUE)
avg.dat$PCA1 <- sp.pca$scores[,1]
avg.dat$PCA2 <- sp.pca$scores[,2]

## STATS
mod <- lm(TotalBiomass ~ PCA1, data = avg.dat)
mod1 <- lm(TotalBiomass ~ PCA2, data = avg.dat)
summary(mod)
summary(mod1)

biplot(sp.pca)
summary(sp.pca)

load.mat <- as.data.frame(sp.pca$loadings[,c(1,2)])
mult <- min(
  (max(avg.dat[,"PCA2"]) - min(avg.dat[,"PCA2"])/(max(load.mat[,"Comp.2"])-min(load.mat["Comp.2"]))),
  (max(avg.dat[,"PCA1"]) - min(avg.dat[,"PCA1"])/(max(load.mat[,"Comp.1"])-min(load.mat["Comp.1"])))
)



load.mat$C1S <- load.mat$Comp.1 * mult*.7
load.mat$C2S <- load.mat$Comp.2 * mult*.7
load.mat$ywords <- ifelse(load.mat$C2 > 0, load.mat$C2, load.mat$C2 - .2)
load.mat$xwords <- ifelse(load.mat$C1 > 0, load.mat$C1 + .5, load.mat$C1 - .75)

names(CDRsp)[1] <- "monospecies"
avg.dat <- plyr::join(avg.dat, CDRsp[,c(1,2)])
avg.dat[1,25] <- "F"
gr1 <- ggplot(aes(x = PCA1, y = PCA2), data = avg.dat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = FunctionalGroup)) +
  #geom_text(aes(x = PCA1, y = PCA2, label = monospecies, vjust = 0, color = FunctionalGroup)) +
  geom_text(data = load.mat, aes(x=xwords +.1, y = ywords + .1, label=rownames(load.mat)), 
            size = 4, vjust=0, color="black", alpha = 0.75) +
  geom_segment(data = load.mat, aes(x=0, y=0, xend=C1S, yend=C2S), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  theme_pubr() +
  labs (x = "PCA1 (62.03%)", y = "PCA2 (15.65%)")

pdf(file = here("Figures","PCAPlot.pdf"), onefile = TRUE, width = 5, height = 5)
print(gr1)
dev.off()


## Attempt at doing PCA through time - but different traits have different coverage 
pl <- vector("list")
ExpYear <- seq(1,19, by = 1)
for (i in 1:19)local({
  i = i
  holddf <- amb.dat[amb.dat$ExpYear == i,]
  holddf <- holddf[,c(1,10:30)] %>%
    group_by(monospecies) %>% 
    summarise_all(funs(mean(., na.rm = TRUE)))
  holddf <- as.data.frame(holddf)
  pca.mat <- holddf
  rownames(pca.mat) <- pca.mat[,1]
  pca.mat  <- pca.mat[,-c(1,21)]
  for(i in 1:length(pca.mat)){
    for(j in 1:nrow(pca.mat)){
      if(is.nan(pca.mat[j,i])){
        pca.mat[j,i] = NA
      }
    }
  }
  sp.pca <- princomp(pca.mat[,c(7:13,15:17)], cor = TRUE, scores = TRUE)
  holddf$PCA1 <- sp.pca$scores[,1]
  holddf$PCA2 <- sp.pca$scores[,2]
  
  load.mat <- as.data.frame(sp.pca$loadings[,c(1,2)])
  mult <- min(
    (max(holddf[,"PCA2"]) - min(holddf[,"PCA2"])/(max(load.mat[,"Comp.2"])-min(load.mat["Comp.2"]))),
    (max(holddf[,"PCA1"]) - min(holddf[,"PCA1"])/(max(load.mat[,"Comp.1"])-min(load.mat["Comp.1"])))
  )
  
  load.mat$C1S <- load.mat$Comp.1 * mult*.7
  load.mat$C2S <- load.mat$Comp.2 * mult*.7
  load.mat$ywords <- ifelse(load.mat$C2 > 0, load.mat$C2, load.mat$C2 - .2)
  load.mat$xwords <- ifelse(load.mat$C1 > 0, load.mat$C1 + .5, load.mat$C1 - .75)
  
  names(CDRsp)[1] <- "monospecies"
  holddf <- plyr::join(holddf, CDRsp[,c(1,2)])
  holddf[1,24] <- "F"
  gr1 <- ggplot(aes(x = PCA1, y = PCA2), data = holddf) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(aes(color = FunctionalGroup)) +
    #geom_text(aes(x = PCA1, y = PCA2, label = monospecies, vjust = 0)) +
    geom_text(data = load.mat, aes(x=xwords +.1, y = ywords + .1, label=rownames(load.mat)), 
              size = 4, vjust=0, color="black", alpha = 0.75) +
    geom_segment(data = load.mat, aes(x=0, y=0, xend=C1S, yend=C2S), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    theme_pubr() +
    labs (x = "PCA1", y = "PCA2")
  print(gr1)
  pl[[i]] <- gr1
  
})



pdf(file = here("Figures","PCAPlot.pdf"), onefile = TRUE, width = 5, height = 5)
print(ggarrange(plotlist = pl, common.legend = TRUE))
dev.off()
