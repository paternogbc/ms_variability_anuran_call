#### Script For Mantel Test (only for Phyllomedusa nordestina)

### Packages:-------------------------------------------------------------------
library(tidyverse)
library(ade4)
library(dplyr)

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

### Data only for Phylomedua nordestina:
mat2 <- mat[complete.cases(mat),]
nordestina <- droplevels(filter(mat2,sp == "Phyllomedusa_nordestina"))
str(nordestina)
head(mat2)

### Average by individuals
med.ind <- summarise(group_by(nordestina,ID_ind),dur = mean(dur),
                     ci = mean(ci), df=mean(df),pr=mean(pr),pl=mean(pl),
                     pn=mean(pn),LAT = mean(LAT),LONG=mean(LONG))
ind.dists <- dist(cbind(med.ind$LONG,med.ind$LAT))

### Average by populations
med.pop <- summarise(group_by(nordestina,ID_pop),dur = mean(dur),
                     ci = mean(ci), df=mean(df),pr=mean(pr),pl=mean(pl),
                     pn=mean(pn),LAT = mean(LAT),LONG=mean(LONG))
pop.dists <- dist(cbind(med.pop$LONG,med.pop$LAT))

as.matrix(ind.dists)
as.matrix(pop.dists)

# Test across individuals--------------------------------------------------------
dur.dists <- dist(med.ind$dur)
df.dists <- dist(med.ind$df)
pr.dists <- dist(med.ind$pr)
pl.dists <- dist(med.ind$pl)
ci.dists <- dist(med.ind$ci)
pn.dists <- dist(med.ind$pn)

# Mantel test
mdur <- mantel.rtest(ind.dists, dur.dists , nrepet = 1000)
mdf  <- mantel.rtest(ind.dists, df.dists , nrepet = 1000)
mpn  <- mantel.rtest(ind.dists, pn.dists , nrepet = 1000)
mpl  <- mantel.rtest(ind.dists, pl.dists , nrepet = 1000)
mpr  <- mantel.rtest(ind.dists, pr.dists , nrepet = 1000)
mci  <- mantel.rtest(ind.dists, ci.dists , nrepet = 1000)

# Mantel table
mantel_table <-
  rbind(ext_mantel(x = mdur, name = "Call duration"),
      ext_mantel(x = mdf, name = "Dominant frequency"),
      ext_mantel(x = mpn, name = "Number of pulses"),
      ext_mantel(x = mpl, name = "Pulse duration"),
      ext_mantel(x = mpr, name = "Pulse rate"),
      ext_mantel(x = mci, name = "Call interval"))
mantel_table$level <- "across individuals"

write_csv(x = mantel_table, path = "outputs/tables/Table_mantel_test_individuals.csv")  

# Test across populations--------------------------------------------------------
dur.dists <- dist(med.pop$dur)
df.dists <- dist(med.pop$df)
pr.dists <- dist(med.pop$pr)
pl.dists <- dist(med.pop$pl)
ci.dists <- dist(med.pop$ci)
pn.dists <- dist(med.pop$pn)

# Mantel test
mdur <- mantel.rtest(pop.dists, dur.dists , nrepet = 1000)
mdf  <- mantel.rtest(pop.dists, df.dists , nrepet = 1000)
mpn  <- mantel.rtest(pop.dists, pn.dists , nrepet = 1000)
mpl  <- mantel.rtest(pop.dists, pl.dists , nrepet = 1000)
mpr  <- mantel.rtest(pop.dists, pr.dists , nrepet = 1000)
mci  <- mantel.rtest(pop.dists, ci.dists , nrepet = 1000)

# Mantel table
mantel_table <-
  rbind(ext_mantel(x = mdur, name = "Call duration"),
        ext_mantel(x = mdf, name = "Dominant frequency"),
        ext_mantel(x = mpn, name = "Number of pulses"),
        ext_mantel(x = mpl, name = "Pulse duration"),
        ext_mantel(x = mpr, name = "Pulse rate"),
        ext_mantel(x = mci, name = "Call interval"))
mantel_table$level <- "across populations"
mantel_table

write_csv(x = mantel_table, path = "outputs/tables/Table_mantel_test_populations.csv")  

