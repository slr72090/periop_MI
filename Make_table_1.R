require(ggplot2)
require(foreign)
require(readstata13)
require(dplyr)
require(tidyr)
require(FactoMineR)
require(PCAmixdata)
require(factoextra)
require(ggfortify)
require(RSQLite)
require(GGally)
require(reshape2)
require(cowplot)
require(corrplot)
require(tidyverse)
require(MASS)
library(proto)
library(pROC)
library(xtable)
select <- dplyr::select

load("./mortality_dispo_models/Table_1_mortality.rda")
load("./predict_MI_model/Table_1_MI.rda")

dat_MI <- dat_by_MI %>% 
  as.data.frame %>% 
  mutate(var = rownames(.)) %>% 
  rename(p.MI = chisq.pval)

dat_mort <- dat_by_mortality %>% 
  as.data.frame %>% 
  mutate(var = rownames(.)) %>% 
  rename(p.mort = chisq.pval)
         
dat <- merge(dat_MI, dat_mort)

dat_tab <- xtable(dat)
print(dat_tab, file="Table_1.txt", include.rownames = F)
