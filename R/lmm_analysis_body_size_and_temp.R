# Mixed models Body size and Temperature

### Packages:-------------------------------------------------------------------
library(cowplot)
library(tidyverse)
library(broom)
library(lme4)
library(DHARMa)
library(sjPlot)
library(car)
library(performance)

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

### Average per individuals
med.ind <- summarise(group_by(mat, ID_ind, ID_sp, ID_pop),
                     dur = mean(dur, na.rm = T),
                     ci = mean(ci, na.rm = T), df = mean(df, na.rm = T), pr = mean(pr, na.rm = T), pl = mean(pl, na.rm = T),
                     pn = mean(pn, na.rm = T), LAT = mean(LAT, na.rm = T),
                     LONG = mean(LONG, na.rm = T), BS = mean(BS, na.rm = T), TE = mean(TE, na.rm = T) )

# Fit models
# Body size
m1 <- lmer(dur ~ BS + (1 | ID_sp), data = med.ind)
m2 <- lmer(ci ~ BS + (1 | ID_sp), data = med.ind)
m3 <- lmer(df ~ BS + (1 | ID_sp), data = med.ind)
m4 <- lmer(pr ~ BS + (1 | ID_sp), data = med.ind)
m5 <- lmer(pl ~ BS + (1 | ID_sp), data = med.ind)
m6 <- lmer(pn ~ BS + (1 | ID_sp), data = med.ind)

# Temperature
# Body size
m7 <- lmer(dur ~ TE + (1 | ID_sp), data = med.ind)
m8 <- lmer(ci ~ TE + (1 | ID_sp), data = med.ind)
m9 <- lmer(df ~ TE + (1 | ID_sp), data = med.ind)
m10 <- lmer(pr ~ TE + (1 | ID_sp), data = med.ind)
m11 <- lmer(pl ~ TE + (1 | ID_sp), data = med.ind)
m12 <- lmer(pn ~ TE + (1 | ID_sp), data = med.ind)

# Organize tables
mod.list <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12)

at <- lapply(mod.list, sl) 
tall <- bind_rows(at)
tall

tall$parameter <- rep(c("Call duration", "Call interval", "Dominant frequency",
                     "Pulse rate", "Pulse duration", "Number of pulses"), times = 2)

write_csv(x = tall, path = "outputs/tables/Table_lmm_body_size_temperature.csv")