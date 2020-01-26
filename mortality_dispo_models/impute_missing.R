## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

# Set the working directory to the source file directory
setwd("~/Desktop/perip_MI_GITHUB")

## Load package scripts -----------------------------------------------
require(ggplot2)
require(dplyr)
require(tidyr)
require(ggfortify)
require(GGally)
require(reshape2)
require(cowplot)
require(corrplot)
require(tidyverse)
require(MASS)
require(mice)
require(VIM)
select <- dplyr::select

data_for_imputation <- data_formatted %>% select(c( age,
                              race,
                              gender,
                              obesity,
                              smoking,
                              alcoholic, 
                              HTN,
                              HLD,
                              hx_DM,
                              hx_ckd,
                              hx_isch_heart,
                              PAD,
                              valve_dz,
                              hx_chf,
                              hx_VTE,
                              chrnlung,
                              malignancy,
                              anemia,
                              hx_CVA,
                              Afib,
                              #sepsis,
                              #PNA,
                              #DVT,
                              #Bleed,
                              NSTEMI#,
                              #invasive_mgmt,
                             # IABP,
                              #cardiogenic_shock,
                              #severe_MI,
                              #hx_revasc,
                              #high_risk_surgery,
                              #ICF
                              ))

## Step 1: Get a look at missing data
tab_missing <- as.data.frame(md.pattern(data_for_imputation))
#aggr_plot <- aggr(data_for_imputation, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

## Step 2: Impute missing values
init = mice(data_for_imputation, maxit=0) 
meth = init$method
predM = init$predictorMatrix

set.seed(103)
imputed = mice(data_for_imputation, method=meth, predictorMatrix=predM, m=5)
complete_data <- complete(imputed)

# Verify completeness of imputed dataset
sapply(complete_data, function(x) sum(is.na(x)))

missing_vars <- data_for_imputation %>% 
  select_if(function(x) any(is.na(x))) %>% 
  names()

