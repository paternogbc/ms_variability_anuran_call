### Sensitivity analysis (species with low sample size)
### Multi-level call variation

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

# Summary of data points across species 
mat %>% 
    group_by(sp) %>% 
    summarise(N_pop = n_distinct(ID_pop),
              N_ind = n_distinct(ID_ind),
              N_can = n_distinct(ID_can)) %>% 
    arrange(N_pop, N_ind, N_can) -> summary_data_points
summary_data_points

# save table
write_csv(x = summary_data_points, 
          path = "outputs/tables/supp/sampling_effort_by_species.csv")

# Remove all species with less than 5 sampled individuals--------------------
# Species to remove
summary_data_points %>% 
    filter(N_ind < 5) %>% 
    pull(sp) %>% as.character() -> sp_remove_1

mat1 <- mat %>% 
    filter(!sp %in% sp_remove_1) %>% 
    droplevels()

# remove NA's
mat1 <- mat1 %>% drop_na(ci)
nrow(mat1)

### 1. Multi-level analysis of variation----------------------------------------
### 1.1 Removing species with less than 5 individuals---------------------------
### Inter-specific variation:---------------------------------------------------
sap <- mat1
table <- NULL
table.cv.5 <- NULL

for (j in 1:1000){
    for (i in 1:5){
        # Select one species at random
        s_sp <- sample(levels(sap$ID_sp),1)
        sp.sap <- droplevels(filter(sap,ID_sp==s_sp))
        # Select one call at random
        s_ca <- sample(levels(sp.sap$ID_can),1)
        ca.sap <- droplevels(filter(sp.sap,ID_can==s_ca))
        # Store parameters selected (one random call from one random species)
        table <- rbind(table,ca.sap)
        # Redifines sap without choosen species
        sap <- droplevels(filter(sap,ID_sp!= s_sp))
    }
    
    ### Calculate CV for each parameter (among 5 calls from different species):
    cv.dur <- cv(table$dur)
    cv.df <- cv(table$df)
    cv.ci <- cv(table$ci)
    cv.pn <- cv(table$pn)
    cv.pr <- cv(table$pr)
    cv.pl <- cv(table$pl)
    
    # Store CV`s
    x <- data.frame(level = "inter_sp" ,cv.dur, cv.df, cv.ci, cv.pn, cv.pr, cv.pl)
    table.cv.5 <- rbind(table.cv.5, x)
    table <- NULL
    sap <- mat1
}

### Intra-specific variation:---------------------------------------------------
sap <- mat1
table <- NULL
table.cv.4 <- NULL

for (j in 1:1000){  
    N.ind <- 1
    # Exclude species with less then 5 individuals
    while(  N.ind   < 5) { 
        
        # Select one species at random
        s_sp <- sample(levels(sap$ID_sp),1)
        N.ind <- length(levels(droplevels(subset(sap,ID_sp==s_sp))$ID_ind))
        sp.sap <- droplevels(filter(sap,ID_sp==s_sp))
    }
    
    # Select five individual from the selected species
    s_in <- sample(levels(sp.sap$ID_ind),5)
    
    # Selects one random call from each individual
    for (k in s_in){
        s_ca <- sample(subset(sp.sap,ID_ind == k)$ID_can,1)
        ca.sap <- sp.sap[sp.sap$ID_can == s_ca,]
        table <- rbind(table,ca.sap)
    }
    ### Calculate CV for each parameter (among 5 individuals from the same
    # species)
    cv.dur <- cv(table$dur)
    cv.df <- cv(table$df)
    cv.ci <- cv(table$ci)
    cv.pn <- cv(table$pn)
    cv.pr <- cv(table$pr)
    cv.pl <- cv(table$pl)
    
    x <- data.frame(level = "intra_sp" ,cv.dur, cv.df, cv.ci, cv.pn, cv.pr, cv.pl)
    table.cv.4 <- rbind(table.cv.4, x)
    table <- NULL
    sap <- mat1
    
}

### Inter-population------------------------------------------------------------
# Subset for Phyllomedusa nordestina
sap <- droplevels(filter(mat,sp=="Phyllomedusa_nordestina" ))
table <- NULL
table.cv.3 <- NULL

for (j in 1:1000){  
    
    # Selects 5 random populations
    s_pop <- sample(levels(sap$ID_pop),5)
    
    # Selects one individual and one call from each population selected
    for (k in s_pop){
        s_in <- sample(levels(droplevels(subset(sap,ID_pop==k))$ID_ind),1)
        s_ca <- sample(subset(sap,ID_ind == s_in)$ID_can,1)
        ca.sap <- sap[sap$ID_can == s_ca,]
        table <- rbind(table,ca.sap)
    }
    
    ### Calculate CV for each parameter (among 5 individuals from different populations
    # of Phyllomedusa nordestina
    
    cv.dur <- cv(table$dur)
    cv.df <- cv(table$df)
    cv.ci <- cv(table$ci)
    cv.pn <- cv(table$pn)
    cv.pr <- cv(table$pr)
    cv.pl <- cv(table$pl)
    
    x <- data.frame(level = "inter_pop" ,cv.dur, cv.df, cv.ci, cv.pn, cv.pr, cv.pl)
    table.cv.3 <- rbind(table.cv.3, x)
    table <- NULL
    sap <- droplevels(filter(mat,sp=="Phyllomedusa_nordestina" ))
    
}

### Intra-population------------------------------------------------------------
# Subset for Phyllomedusa nordestina & Phyllomedusa_hypochondrialis
sap <- droplevels(filter(mat,sp=="Phyllomedusa_nordestina" ))
table <- NULL
table.cv.2 <- NULL

# Exclude population with less then 5 individuals
for (j in 1:1000){  
    N.ind <- 1
    while(  N.ind   < 5) { 
        # Select one population from Phyllomedusa nordestina
        s_pop <- sample(levels(sap$ID_pop),1)
        N.ind <- length(levels(droplevels(subset(sap,ID_pop==s_pop))$ID_ind))
    }
    
    # Selects five random individuals from the selected population
    s_in <- sample(levels(droplevels(subset(sap,ID_pop==s_pop))$ID_ind),5)
    
    # Selects one random call from each individual
    for (k in s_in){
        s_ca <- sample(subset(sap,ID_ind == k)$ID_can,1)
        ca.sap <- sap[sap$ID_can == s_ca,]
        table <- rbind(table,ca.sap)
    }
    
    ### Calculate CV for each parameter (among 5 individuals from the same populations
    # of Phyllomedusa nordestina
    cv.dur <- cv(table$dur)
    cv.df <- cv(table$df)
    cv.ci <- cv(table$ci)
    cv.pn <- cv(table$pn)
    cv.pr <- cv(table$pr)
    cv.pl <- cv(table$pl)
    
    x <- data.frame(level = "intra_pop" ,cv.dur, cv.df, cv.ci, cv.pn, cv.pr, cv.pl)
    table.cv.2 <- rbind(table.cv.2, x)
    table <- NULL
    sap <-  droplevels(filter(mat,sp=="Phyllomedusa_nordestina" ))
    
}

### Intra-individual variation-------
sap <- mat1
table <- NULL
table.cv.1 <- NULL

for (j in 1:1000){  
    
    # Exclude species with less then 5 individuals
    N.ind <- 1
    while(  N.ind   < 5) { 
        
        # Select one species at random
        s_sp <- sample(levels(sap$ID_sp),1)
        N.ind <- length(levels(droplevels(subset(sap,ID_sp==s_sp))$ID_ind))
        sp.sap <- droplevels(filter(sap,ID_sp==s_sp))
    }
    
    # Exclude individuals with less then 5 calls
    N.can <- 1
    while(  N.can   < 5) { 
        # Sample on random individual:
        s_in <- sample(levels(sp.sap$ID_ind),1)
        N.can <- length(levels(droplevels(subset(sap,ID_ind==s_in))$ID_can))  
    }
    
    # Sample 5 calls at random from the selected individual
    s_ca <- sample(subset(sp.sap,ID_ind == s_in)$ID_can,5)
    
    for (i in 1:5){
        ca.sap <- sp.sap[sp.sap$ID_can == s_ca[i],]
        table <- rbind(table,ca.sap)
    }
    
    ### Calculate CV for each parameter (among 5 calls from the same individual from 
    # any population or species in the dataset
    cv.dur <- cv(table$dur)
    cv.df <- cv(table$df)
    cv.ci <- cv(table$ci)
    cv.pn <- cv(table$pn)
    cv.pr <- cv(table$pr)
    cv.pl <- cv(table$pl)
    
    x <- data.frame(level = "intra_ind" ,cv.dur, cv.df, cv.ci, cv.pn, cv.pr, cv.pl)
    table.cv.1 <- rbind(table.cv.1, x)
    table <- NULL
    sap <- mat1
    
}

### Integrating results:
table.cv <- rbind(table.cv.1,
                  table.cv.2,
                  table.cv.3,
                  table.cv.4,
                  table.cv.5)

str(table.cv)

table.cv.melt <- melt(table.cv)

str(table.cv)
str(table.cv.melt)

### Re-name factors level:
levels(table.cv.melt$variable) <- c(
    "Call duration","Dominant frequency","Call interval",
    "Number of pulses","Pulse rate","Pulse duration")

# Re-order factors levels:
table.cv.melt$variable <- factor(table.cv.melt$variable,  c(
    "Dominant frequency","Number of pulses","Pulse duration",
    "Pulse rate", "Call duration", "Call interval"))

saveRDS(object = table.cv.melt, file = "outputs/temp/sensi_data_1000x_cv_less_than_5_ind.Rds")

table.cv.melt <- readRDS(file = "outputs/temp/sensi_data_1000x_cv_less_than_5_ind.Rds")

# CV Plot sensi 1 -----------------------------------------------------------------
g1 <- ggplot(table.cv.melt, aes(y=value,x=level,fill=level))+
    geom_jitter(colour = "gray", alpha =.3, size = .05, width = .1)+
    geom_boxplot(alpha = .7, outlier.colour = NA, size = .3)+
    geom_hline(yintercept = 100, colour = "red", linetype = "dotted")+
    scale_x_discrete(labels = c("df", "pn", "pl", "pr", "cd", "ci")) +
    ylab("Coefficient of variantion (%)")+
    xlab("")+
    facet_grid(~variable) +
    theme_classic(base_size = 9) +
    scale_y_continuous(breaks=seq(0,250,25),limits = c(0, 250))+
    theme(legend.position = c(0.08, 0.8),
          legend.background=element_rect(fill="transparent"),
          axis.text.x=element_blank(),
          panel.background=element_rect(fill = "white", colour = "black"),
          legend.text = element_text(size = 6),
          legend.key.size = unit(4, "mm"),
          strip.background = element_rect(
              color = "black", size = 0.5))+
    scale_fill_brewer(palette = 8, labels=c("Intra-individual",
                                            "Intra-population*",
                                            "Inter-population*",
                                            "Intra-species",
                                            "Inter-species")) +
    labs(title = "Sensitivity analysis: removing species with low sample size",
         subtitle = "Removed species: P. burmeisteri, P. megacephala, P. vaillantii and P. azurea"); g1

ggsave(filename = "outputs/figures/supp/SFigure_sensi_cv_less_than_5_ind_cv.pdf", plot = g1, 
       width = 180, height = 90, units = 'mm')
ggsave(filename = "outputs/figures/supp/SFigure_sensi_cv_less_than_5_ind_cv.png", plot = g1, 
       width = 180, height = 90, units = 'mm')

# 2. Variance partition---------------------------------------------------------
### Nested Anova----------------------------------------------------------------
# dur-------------
mod <- lm(dur ~ sp/pop/ind/can , mat1)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.dur <- c(ss.sp,ss.pop,ss.ind,ss.can)

# df--------------
mod <- lm(df ~ sp/pop/ind/can,mat1)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.df <- c(ss.sp,ss.pop,ss.ind,ss.can)

# ci---------------
mod <- lm(ci ~ sp/pop/ind/can,mat1)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.ci <- c(ss.sp,ss.pop,ss.ind,ss.can)

# pn-----------------
mod <- lm(pn ~ sp/pop/ind/can,mat1)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.pn <- c(ss.sp,ss.pop,ss.ind,ss.can)

# pr--------------
mod <- lm(pr ~ sp/pop/ind/can,mat1)
sst <- sum(anova(mod)[2][1:4,1])
ss.sp <- anova(mod)[2][1,1]/sst
ss.pop <- anova(mod)[2][2,1]/sst
ss.ind <- anova(mod)[2][3,1]/sst
ss.can <- anova(mod)[2][4,1]/sst
var.explained.pr <- c(ss.sp,ss.pop,ss.ind,ss.can)

# pl--------------
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
saveRDS(nested, file = "outputs/temp/sensi_data_nested_variance.RDs'")

# Load data
nested <- readRDS(file = "outputs/temp/sensi_data_nested_variance.RDs'")

# Variance partition plot-------------------------------------------------------
g2 <-
    ggplot(nested, aes(x = param, y = variance,fill = level))+
    geom_bar(stat="identity") + 
    scale_fill_grey(name = "", labels = c("call", "individual", "population", "species")) +
    scale_x_discrete(labels = c("Dominant frequency",  "Number of pulses", "Pulse duration", "Pulse rate",
                                "Call duration", "Call interval")) +
    xlab("Acoustic parameter")+
    ylab("Variance explained (%)") +
    coord_flip() +
    labs(title = "Sensitivity analysis: removing species with low sample size",
         subtitle = "Removed species: P. burmeisteri, P. megacephala, P. vaillantii and P. azurea") +
    theme_classic(base_size = 8) +
    theme(legend.position = "top" , legend.key.size = unit(3, "mm"), 
          legend.text = element_text(size = 5.5),
          legend.spacing = unit(.1, "mm"),
          title = element_text(size = 4));g2

ggsave(filename = "outputs/figures/supp/SFigure_sensi_variance_partition.png", plot = g2, width = 80, height = 60, units = "mm")
ggsave(filename = "outputs/figures/supp/SFigure__sensi_variance_partition.pdf", plot = g2, width = 80, height = 60, units = "mm")
