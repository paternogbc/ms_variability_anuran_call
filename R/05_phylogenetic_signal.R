### Phylogenetyic signal of acoustic parameters

### Packages:-------------------------------------------------------------------
library(tidyverse)
library(picante)

### Start-----------------------------------------------------------------------
rm(list = ls())

# Own functions -----------------------------------------------------------
source("R/zzz_functions.R")

### Load data-------------------------------------------------------------------
mat <- read.csv("data/raw/raw_data_anuran_call_parameters.csv", h = T, sep = ",")
str(mat)

### Phylogeny entry:
phy <- read.tree("data/raw/amph_2014.tre")
str(phy)

# Transforming ID`s intos factors:
mat$ID_can <- as.factor(mat$ID_can)
mat$ID_ind <- as.factor(mat$ID_ind)
mat$ID_pop <- as.factor(mat$ID_pop)
mat$ID_sp <- as.factor(mat$ID_sp)
str(mat)

# Summary for species mean:
mat.sum <- summarise(group_by(mat,sp),dur=mean(dur),
                     df = mean(df), ci = mean(ci, na.rm = T), pn = mean(pn),
                     pr = mean(pr), pl = mean(pl), BS = mean(BS, na.rm = T), 
                     TE = mean(TE, na.rm = T))
rownames(mat.sum) <- mat.sum$sp
mat.sum <- as.data.frame(mat.sum)

# Pruned phylogeny:
tree.drop  <- drop.tip(phy,rownames(mat.sum))              
tree <- (drop.tip(phy,tree.drop$tip.label))
tree$node.label <- makeLabel(tree)$node.label 

### Checking for absent species in data:
sp.match <- match(tree$tip.label, rownames(mat.sum))
### Species absent from phylogeny:
as.character(mat.sum[-sp.match,]$sp)

# Phylogenetic signal-----------------------------------------------------------
k.signal <- multiPhylosignal(x = mat.sum[, 2:8], phy = tree, reps = 10000)
k.signal <- as.matrix(k.signal)

k.signal <- rownames_to_column(data.frame(k.signal), var = "parameter")

# Save table
write_csv(x = k.signal, path = "outputs/tables/Table_phylogenetic_signal.csv")
