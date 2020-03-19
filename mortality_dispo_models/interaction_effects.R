## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

# Set the working directory to the source file directory
setwd("~/Desktop/perip_MI_GITHUB/mortality_dispo_models/")

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
library(survey)
select <- dplyr::select

#Plotting specs
textSize = 12
save_plots = F
source("../utility_scripts/plot_themes.R")

#Data saving
data_filename = "Data_objects/data_imputed_with_weights.rda" 
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
           high_risk_surgery = as.factor(as.numeric(transplant == 1|thoracic_surgery == 1|vascular == 1|abdominal == 1)),
           high_risk_surgery_2 = as.factor(as.numeric(transplant == 1|thoracic_surgery == 1|vascular == 1)),
           hx_revasc = as.factor(as.numeric(prior_CABG == 1 | prior_PCI ==1))) %>% 
    mutate(RCRI_pt = as.factor(as.numeric(RCRI_pt) + as.numeric(high_risk_surgery == 1))) %>% 
    mutate(severe_MI = as.factor(as.numeric(cardiogenic_shock == 1 | IABP == 1))) %>% 
    select(-c(prior_CABG, prior_PCI, prior_MI, CAD)) %>% #transplant,thoracic_surgery,vascular,)) %>% 
    mutate(`RCRI >= 3` = as.factor(as.numeric(RCRI_pt) > 3))
  
  if(impute_missing){
    source("impute_missing.R")
    data_formatted[,missing_vars] <- complete_data[,missing_vars]
  }
  
  data = data_formatted
  
  if(save_data){
    save(data_formatted, data_all, file = data_filename)
  }
}

if(generate_data == F){
  load(data_filename)
  data <- data_formatted 
  rm(data_formatted)
}

options(survey.lonely.psu="adjust")
#var_names_full <- c(names(X1), names(X2), "nchronic", "cm_mets", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI", "invasive_mgmt")#, "nchronic")
#dat_full <- data %>% select(c(year, died, trendwt, hospid, nis_stratum, var_names_full))
survey_df_full <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = data, nest =T)


## Interaction effects for mortality model ## ---------------------------------------------------------------------
out.RCRIplus <- svyglm(as.numeric(died) ~ 
                         age +
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
                         hx_DM + 
                         hx_ckd + 
                         hx_CVA + 
                         hx_isch_heart+ 
                         #hx_revasc + 
                         high_risk_surgery + 
                         #surgery_type + 
                         invasive_mgmt +
                         NSTEMI + 
                         severe_MI ,
                       design = survey_df_full,
                       family = "binomial")

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
                 hx_revasc*hx_isch_heart +
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

names(mod_coefs5)<- c("OR", "2.5% CI", "97.5% CI")
coef_df <- data.frame(Var = rownames(mod_coefs5),
                      OR = paste0(mod_coefs5$OR, " [", mod_coefs5$"2.5% CI", ", ",  mod_coefs5$"97.5% CI", "]"))

if(output_tables){
  tab_coefs <- xtable(coef_df)
  print(tab_coefs, file="mortality_model_interaction_coefs.txt", include.rownames = F)
}

## Interaction effects for discharge to ICF model ## ---------------------------------------------------------------------
out.full.dispo <- svyglm(ICF ~ year + 
                     age + 
                     race + 
                     gender + 
                     obesity + 
                     smoking + 
                     alcoholic + 
                     HTN + 
                     HLD + 
                     hx_DM + 
                     hx_ckd + 
                     hx_isch_heart + 
                     PAD + 
                     valve_dz + 
                     hx_chf + 
                     hx_VTE + 
                     chrnlung + 
                     malignancy + 
                     cm_mets + 
                     anemia +
                     hx_CVA + 
                     Afib + 
                     sepsis + 
                     PNA + 
                     PE + 
                     DVT + 
                     Bleed + 
                     NSTEMI + 
                     invasive_mgmt + 
                     high_risk_surgery +
                     severe_MI,
                   design = survey_df_full,
                   family = "binomial")


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

names(mod_coefs5.dispo)<- c("OR", "2.5% CI", "97.5% CI")
coef_df <- data.frame(Var = rownames(mod_coefs5.dispo),
                      OR = paste0(mod_coefs5.dispo$OR, " [", mod_coefs5.dispo$"2.5% CI", ", ",  mod_coefs5.dispo$"97.5% CI", "]"))

if(output_tables){
  tab_coefs <- xtable(coef_df)
  print(tab_coefs, file="dispo_model_interaction_coefs.txt", include.rownames = F)
}
