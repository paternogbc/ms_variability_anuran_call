### Plot tree and populations into the map

### Packages:-------------------------------------------------------------------
library(phytools)
library(mapdata)
library(sensiPhy)
library(tidyverse)

### Start-----------------------------------------------------------------------
rm(list = ls())

# Own functions -----------------------------------------------------------
source("R/zzz_functions.R")

### Load data-------------------------------------------------------------------
mat <- read.csv("data/raw/raw_data_anuran_call_parameters.csv", h = T, sep = ",")

### Phylogeny entry:
phy <- read.tree("data/raw/amph_2014.tre")

# Prepare data for ploting
# Get individuals level data
  mat %>% group_by(sp, ID_ind) %>% 
  select(LAT, LONG) %>% 
  summarise(lat = unique(LAT), long = unique(LONG)) %>% 
  ungroup() %>% 
  select(sp, lat, long) %>% 
  as.data.frame() -> mat.ind

sps <- as.character(mat.ind$sp)
bra <- as.matrix(mat.ind[, 2:3])
row.names(bra) <- sps

# Match data and phy
mat %>% group_by(ID_sp) %>% 
  summarise(sp = unique(sp)) %>% 
  as.data.frame() -> species
rownames(species) <- species$sp

tree <- sensiPhy::match_dataphy(formula = sp ~ 1, data = species, phy = phy)$phy


# Plot phylo to map------------------------------------------------------------ 
obj <- phylo.to.map(tree, bra, database = "worldHires",rotate = TRUE,
                  regions = "Brazil", plot = FALSE, type = "phylogram")

cols<-setNames((rainbow(n = Ntip(tree)))[15:1], tree$tip.label)

pdf(width = 9, height = 7, file = "outputs/figures/supp/SFigure_phylo_to_map.pdf")

plot(obj, direction="rightwards", colors = cols, cex.points=c(0,1.5), rotate = TRUE,
     pts = F, fsize = .8, ylim=c(-40, 5), xlim=c(-80, -30),
     mar=c(5.1,4.1,2.1,2.1))

# add lat + long lines
long.lines<-seq(-80, 30, by = 10)
lat.lines<-seq(-40, 10, by = 10)
axis(1, at=long.lines, labels = paste(long.lines, "°", sep = ""))
axis(2, at=lat.lines, labels = paste(lat.lines, "°", sep = ""))
## add gridlines
for(i in 1:length(lat.lines))
  lines(c(par()$usr[1],160),rep(lat.lines[i],2),lty="dotted",
        col="#00000040")
for(i in 1:length(long.lines))
  lines(rep(long.lines[i],2),c(par()$usr[3],-10),lty="dotted",
        col="#00000040")
title(xlab="Longitude")
title(ylab="Latitude ")        

dev.off()
