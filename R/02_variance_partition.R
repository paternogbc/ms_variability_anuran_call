### Variance partition

### Packages:-------------------------------------------------------------------
library(tidyverse)
library(reshape2)
library(cowplot)

### Start-----------------------------------------------------------------------
rm(list = ls())

# Own functions -----------------------------------------------------------
source("R/zzz_functions.R")

### Load data-------------------------------------------------------------------
mat <- read.csv("data/raw/raw_data_anuran_call_parameters.csv", h = T, sep = ",")
str(mat)

# Transforming ID`s intos factors:
mat$ID_can <- as.factor(mat$ID_can)
mat$ID_ind <- as.factor(mat$ID_ind)
mat$ID_pop <- as.factor(mat$ID_pop)
mat$ID_sp <- as.factor(mat$ID_sp)
str(mat)

### Nested Anova----------------------------------------------------------------
# dur
mod <- lm(dur ~ sp/pop/ind/can , mat)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.dur <- c(ss.sp,ss.pop,ss.ind,ss.can)

# df
mod <- lm(df ~ sp/pop/ind/can,mat)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.df <- c(ss.sp,ss.pop,ss.ind,ss.can)

# ci
mod <- lm(ci ~ sp/pop/ind/can,mat)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.ci <- c(ss.sp,ss.pop,ss.ind,ss.can)

# pn
mod <- lm(pn ~ sp/pop/ind/can,mat)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.pn <- c(ss.sp,ss.pop,ss.ind,ss.can)

# pr
mod <- lm(pr ~ sp/pop/ind/can,mat)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.pr <- c(ss.sp,ss.pop,ss.ind,ss.can)

# pl
mod <- lm(pl ~ sp/pop/ind/can,mat)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.pl <- c(ss.sp,ss.pop,ss.ind,ss.can)

### Organizing data:
level <- rep(c("sp","pop","ind","can"),each=1,times=6)
param <- rep(c("dur","df","ci","pn","pr","pl"),each=4)
variance <- c(var.explained.dur,var.explained.df,var.explained.ci,
              var.explained.pn,var.explained.pr,var.explained.pl)
nested <- data.frame(level=level,param=param,variance=variance)
nested$param <- factor(nested$param, levels =  c("df","pn","pl","pr","dur","ci"))

### Save temp data
saveRDS(nested, file = "outputs/temp/data_nested_variance.RDs'")

# Load data----
nested <- readRDS(file = "outputs/temp/data_nested_variance.RDs'")

g1 <-
    ggplot(nested, aes(x = param, y = variance,fill = level))+
    geom_bar(stat="identity") + 
    scale_fill_grey(name = "", labels = c("call", "individual", "population", "species")) +
    scale_x_discrete(labels = c("Dominant frequency",  "Number of pulses", "Pulse duration", "Pulse rate",
                                "Call duration", "Call interval")) +
    theme_classic(base_size = 8) +
    theme(legend.position = "top" , legend.key.size = unit(3, "mm"), 
          legend.text = element_text(size = 5.5),
          legend.spacing = unit(.1, "mm")) +
    xlab("Acustic parameter")+
    ylab("Variance explained (%)") +
    coord_flip(); g1

ggsave(filename = "outputs/figures/Figure_variance.png", plot = g1, width = 80, height = 60, units = "mm")
ggsave(filename = "outputs/figures/Figure_variance.pdf", plot = g1, width = 80, height = 60, units = "mm")
