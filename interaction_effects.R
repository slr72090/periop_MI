## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

# Set the working directory to the source file directory
setwd("~/Desktop/perip_MI_GITHUB")

## Load package scripts -----------------------------------------------
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

#Plotting specs
textSize = 12
save_plots = F
source("plot_themes.R")

## Get data ## ---------------------------------------------------------------------
#Data saving
data_filename = "data_all_raw.rda"
save_data = F
drop_missing = F
impute_missing = T
generate_data = F
output_tables = F

if(generate_data == T){
  # Identify procedures that correspond to non-cardiac surgeries
  year_vec = c(2008:2013)
  n_procedures_vec <- rep(15,length(year_vec)) #rep(15,length(year_vec)) # Number of possible procedures per individual for this year 
  n_dx_vec <- c(15, 25, 25, 25, 25, 25) # Number of possible diagnoses per individual for this year 
  source("process_data.R")
  
  data_formatted <- data_all %>% 
    rename(obesity = cm_obese,
           alcoholic = cm_alcohol,
           HTN = cm_htn_c,
           valve_dz = cm_valve,
           chrnlung = cm_chrnlung,
           anemia = cm_anemdef,
           PAD = cm_perivasc) %>% 
    mutate(age_factor = as.factor(ntile(age,3)),
           age = as.numeric(scale(age)),
           nchronic = as.numeric(scale(nchronic)),
           invasive_mgmt = as.factor(invasive_mgmt),
           high_risk_surgery = as.factor(as.numeric(transplant == 1|thoracic_surgery == 1|vascular == 1)),
           hx_revasc = as.factor(as.numeric(prior_CABG == 1 | prior_PCI ==1))) %>% 
    mutate(RCRI_pt = as.factor(as.numeric(RCRI_pt) + as.numeric(high_risk_surgery == 1))) %>% 
    mutate(severe_MI = as.factor(as.numeric(cardiogenic_shock == 1 | IABP == 1))) %>% 
    select(-c(prior_CABG, prior_PCI, prior_MI, CAD, transplant,thoracic_surgery,vascular,nchronic)) %>% 
    mutate(`RCRI >= 3` = as.factor(as.numeric(RCRI_pt) > 3))
  
  if(impute_missing){
    source("impute_missing.R")
    data_formatted[,missing_vars] <- complete_data[,missing_vars]
  }
  
  data = data_formatted
  
  if(save_data){
    save(data, file = data_filename)
  }
}

if(generate_data == F){
  load(data_filename)
  source("impute_missing.R")
  data_formatted[,missing_vars] <- complete_data[,missing_vars]
  data <- data_formatted 
}
## Visualize correlations ## --------------------------------------------------
cor_data <- data %>% select(-c(ind, 
                               year, 
                               invasive_mgmt, 
                               Afib,
                               sepsis,
                               PNA, 
                               PE,
                               DVT,
                               Bleed,
                               died, 
                               ICF,
                               severe_MI,
                               cardiogenic_shock,
                               IABP,
                               RCRI_pt, 
                               age, 
                               `RCRI >= 3`, 
                               cm_mets)) %>% 
  apply(.,2,as.numeric)
M <- cor(cor_data)
diag(M) <- 0
p_cor <- corrplot(M,
                  method = "shade",
                  type = "lower", 
                  order = "hclust",
                  tl.col = "black",
                  tl.cex = .75)

cor_data_acute <- data %>% select(c(age_factor, 
                                    NSTEMI,
                                    severe_MI,
                                    Bleed, 
                                    PE,
                                    DVT,
                                    Afib,
                                    PNA,
                                    sepsis)) %>% 
  mutate(STEMI = 1- as.numeric(NSTEMI)) %>% 
  apply(.,2,as.numeric) 

M_acute <- cor(cor_data_acute)
diag(M_acute) <- 0
p_cor_acute <- corrplot(M_acute,
                        method = "shade",
                        type = "lower",
                        diag=T,
                        order = "hclust",
                        tl.col = "black",
                        tl.cex = .75)


## Interaction effects for mortality model ## ---------------------------------------------------------------------
out.RCRIplus <- glm(as.numeric(died) ~  age +
                      gender + 
                      race + 
                      year + 
                      Afib + 
                      sepsis +
                      PNA + 
                      PE + 
                      DVT + 
                      Bleed + 
                      hx_chf +
                      hx_CVA + 
                      hx_DM + 
                      hx_ckd + 
                      high_risk_surgery + 
                      invasive_mgmt +
                      NSTEMI + 
                      hx_isch_heart + 
                      severe_MI, 
                    data = data %>% filter(ICF == 0), family = "binomial")

mod_coefs1 <- data.frame(apply(exp(cbind(OR = coef(out.RCRIplus), confint(out.RCRIplus))),2,round,2))
# Demographic factors
mod2 <- update(out.RCRIplus, . ~ . + age*race + age*gender + age*hx_DM)
mod_coefs2 <- data.frame(apply(exp(cbind(OR = coef(mod2), confint(mod2))),2,round,2))
# History of ischemic heart disease
mod3 <- update(out.RCRIplus, .~. + hx_isch_heart*age + hx_isch_heart*gender + hx_isch_heart*NSTEMI + hx_isch_heart*high_risk_surgery)
mod_coefs3 <- data.frame(apply(exp(cbind(OR = coef(mod3), confint(mod3))),2,round,2))
# Features of MI and peri-operative conditions 
mod4 <- update(out.RCRIplus, . ~ . + NSTEMI*Bleed + NSTEMI*age*gender*invasive_mgmt )
mod_coefs4 <- data.frame(apply(exp(cbind(OR = coef(mod4), confint(mod4))),2,round,2))

#Aggretgate
mod5 <- update(out.RCRIplus, . ~ . + age*race + 
                 age*gender + 
                 age*hx_DM + 
                 hx_isch_heart*age + 
                 hx_isch_heart*gender + 
                 hx_isch_heart*NSTEMI + 
                 hx_isch_heart*high_risk_surgery + 
                 NSTEMI*Bleed + 
                 NSTEMI*age  
                 )

mod_coefs5 <- data.frame(apply(exp(cbind(OR = coef(mod5), confint(mod5))),2,round,2))

interaction_mod_mortality_BIC <- data.frame(Model = c(1:5),
                                   BIC = c(BIC(out.RCRIplus),
                                           BIC(mod2),
                                           BIC(mod3),
                                           BIC(mod4),
                                           BIC(mod5))) %>% 
                                     mutate(delta = BIC - min(BIC))

## Interaction effects for discharge to ICF model ## ---------------------------------------------------------------------
data$ICF <- as.numeric(data$ICF==1)
var_names_full <- c("age", "race", "gender", "obesity","smoking", "alcoholic",
                    "HTN", "HLD", "hx_DM", "hx_ckd", "hx_isch_heart", "hx_revasc", "PAD", "valve_dz",
                    "hx_chf","hx_VTE", "chrnlung","malignancy","anemia","hx_CVA", "high_risk_surgery",
                    "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI")
data_full_dispo <- data %>% #filter(died == 0) %>% 
  select(c(year, ICF, died, invasive_mgmt, var_names_full)) 
out.full.dispo <- glm(ICF ~. , data = data_full_dispo, family = "binomial")
mod_coefs1.dispo <- data.frame(apply(exp(cbind(OR = coef(out.full.dispo), confint(out.full.dispo))),2,round,2))

# Demographic factors 
mod2.dispo <- update(out.full.dispo, . ~ . + age*race + age*smoking + age*alcoholic + age*obesity + age*gender + age*hx_DM)
mod_coefs2.dispo <- data.frame(apply(exp(cbind(OR = coef(mod2.dispo), confint(mod2.dispo))),2,round,2))

# History of ischemic heart disease
mod3.dispo <- update(out.full.dispo, .~. + hx_isch_heart*gender + 
                       hx_isch_heart*age +
                       hx_isch_heart*gender + 
                       hx_isch_heart*NSTEMI + 
                       hx_isch_heart*HLD + 
                       hx_isch_heart*HTN + 
                       hx_isch_heart*high_risk_surgery)
mod_coefs3.dispo <- data.frame(apply(exp(cbind(OR = coef(mod3.dispo), confint(mod3.dispo))),2,round,2))

# Features of MI and peri-operative conditions 
mod4.dispo <- update(out.full.dispo, . ~ . + NSTEMI*Bleed + 
                       NSTEMI*age + 
                       NSTEMI*hx_isch_heart + 
                       NSTEMI*age*gender*invasive_mgmt)
mod_coefs4.dispo <- data.frame(apply(exp(cbind(OR = coef(mod4.dispo), confint(mod4.dispo))),2,round,2))

#Aggretgate
mod5.dispo <- update(out.full.dispo, . ~ . + age*race + 
                 age*smoking + 
                 age*alcoholic + 
                 age*obesity + 
                 age*gender + 
                 age*hx_DM + hx_isch_heart*gender + 
                 hx_isch_heart*age +
                 hx_isch_heart*gender + 
                 hx_isch_heart*NSTEMI + 
                 hx_isch_heart*HLD + 
                 hx_isch_heart*HTN + 
                 hx_isch_heart*high_risk_surgery + 
                 NSTEMI*Bleed + 
                 NSTEMI*age + 
                 NSTEMI*hx_isch_heart + 
                 NSTEMI*age*gender*invasive_mgmt
)

mod_coefs5.dispo <- data.frame(apply(exp(cbind(OR = coef(mod5.dispo), confint(mod5.dispo))),2,round,2))

interaction_mod_dispo_BIC <- data.frame(Model = c(1:5),
                                            BIC = c(BIC(out.full.dispo),
                                                    BIC(mod2.dispo),
                                                    BIC(mod3.dispo),
                                                    BIC(mod4.dispo),
                                                    BIC(mod5.dispo))) %>% 
  mutate(delta = BIC - min(BIC))
