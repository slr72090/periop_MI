## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

# Set the working directory to the source file directory
setwd("~/Desktop/perip_MI_GITHUB/predict_MI_model/")

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
data_filename =  "data_all_with_weights.rda"
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

## Logit Regression analysis ## ------------------------------------------------------------------------------

options(survey.lonely.psu="adjust")
data$nchronic = scale(data$nchronic)
data$MI = as.numeric(data$MI == 1)
survey_df_full <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = data, nest =T)

out.full <- svyglm(as.numeric(MI) ~ year + 
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
                     nchronic + 
                     hx_CVA + 
                     high_risk_surgery,
                   design = survey_df_full,
                   family = "binomial")

confidence_intervals =  confint.default(out.full, adjust = "bonferroni")
mod_coefs <- data.frame(apply(exp(cbind(OR = coef(out.full),confidence_intervals)),2,round,2))
names(mod_coefs)<- c("OR", "2.5% CI", "97.5% CI")
coef_df <- data.frame(Var = rownames(mod_coefs),
                            OR = paste0(mod_coefs$OR, " [", mod_coefs$"2.5% CI", ", ",  mod_coefs$"97.5% CI", "]"))

if(output_tables){
  tab_coefs <- xtable(coef_df)
  print(tab_coefs, file="MI_model_coefs.txt")
}

## Testing model out of sample ## -------------------------------------------

set.seed(364)
sample <- sample(nrow(data),floor(nrow(data)*0.75))
train <- data[sample,]
test <- data[-sample,]

library(pROC)
out.RCRIplus.train <- predict(out.RCRIplus, newdata = train)
out.full.train <- predict(out.full, newdata = train)

test_prob_RCRIplus = predict(out.RCRIplus, newdata = test, type = "response")
test_roc_RCRIplus = roc(test$died ~ test_prob_RCRIplus, plot = TRUE, print.auc = TRUE)

test_prob_full = predict(out.full, newdata = test, type = "response")
test_roc_full = roc(test$died ~ test_prob_full, plot = TRUE, print.auc = TRUE)

roc.test(test_roc_full, test_roc_RCRIplus)
p_roc <- ggroc(list(`Full model` = test_roc_full, `RCRI model` = test_roc_RCRIplus), linetype = 2) +
  plot_themes + 
  labs(color = "Model") + 
  ggtitle("ROC curves")

if(save_plots){
  save_plot("ROC_curves_mortality.pdf", p_roc, base_width = 8, base_height = 4)
}

## Summary statistics, and evaluate dependnece of ROC on percent test data #### -----------------------
AUC_RCRIplus <- ci.auc(test_roc_RCRIplus, method = "bootstrap")
AUC_full <- ci.auc(test_roc_full, method = "bootstrap")

## Prediction by number of RCRI factors ## ------------------------------------------------------------
dat_sub <- data %>% filter(!is.na(died))
pred.test <- predict(out.RCRIplus, 
                      newdata = dat_sub,
                      type = "response")

dat_sub_agg <- dat_sub %>% select(ind, died, NSTEMI, gender, invasive_mgmt, hx_chf, hx_isch_heart, hx_CVA, hx_DM, hx_ckd, high_risk_surgery) %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate_if(is.character, as.numeric) %>% 
  mutate(RCRI = hx_chf + hx_isch_heart + hx_CVA + hx_DM + hx_ckd + high_risk_surgery) %>% 
  mutate(RCRI2 = RCRI) %>% 
  mutate(prob = pred.test) 

dat_sub_agg[dat_sub_agg$RCRI2 > 4,]$RCRI2 <- 4

df <- dat_sub_agg %>% 
  group_by(RCRI2, invasive_mgmt, hx_isch_heart) %>% 
  summarise(med = median(prob),
            q1 = quantile(prob,.25),
            q2 = quantile(prob,.75))

df$invasive_mgmt = as.factor(df$invasive_mgmt)
levels(df$invasive_mgmt)<- c("No invasive management", "Invasive management")
df$hx_isch_heart = as.factor(df$hx_isch_heart)
levels(df$hx_isch_heart)<- c("No prior CAD", "Prior CAD")

p_risk <- ggplot(df, aes(x = RCRI2, y = med)) +  
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = q1, ymax = q2), linetype = 2) + 
  facet_wrap(.~invasive_mgmt + hx_isch_heart, scales = "free_y") +
  xlab("RCRI") + 
  ylab("Median probability (IQR)") +
  plot_themes
  
if(save_plots){
  save_plot("Risk_by_RCRI_and_CAD.pdf", p_risk, base_width = 8, base_height = 6)
}


############################################################-------------------------------------------------------------------------------
## Analysis of disposition to intermediate care facility ## ------------------------------------------------------------------------------
############################################################-------------------------------------------------------------------------------

data$ICF <- as.numeric(data$ICF==1)
sub.PCA.dispo <- data[,(colnames(data) %in% c("ICF", "hospid","nis_stratum", "year", "race", "gender", "invasive_mgmt", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI", unlist(retained_vars)))]
dat_full_dispo <- data %>% filter(died == 0) %>% select(c(year, ICF, trendwt, hospid, nis_stratum, var_names_full)) 
dat_PCA <- data[,(colnames(data) %in% c("hospid", "nis_stratum", "trendwt", "ICF", "year", "race", "gender", "invasive_mgmt", "PNA", "sepsis", "Afib", "PE", "DVT",  "Bleed", "NSTEMI", "severe_MI", "trendwt", unlist(retained_vars)))]

survey_df_full_dispo <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = dat_full_dispo, nest =T)

out.full.dispo <- svyglm(as.numeric(ICF) ~ year + 
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
                   design = survey_df_full_dispo,
                   family = "binomial")

out.RCRIplus.dispo <- svyglm(as.numeric(ICF) ~ 
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
                         high_risk_surgery + 
                         invasive_mgmt +
                         NSTEMI + 
                         severe_MI ,
                       design = survey_df_full_dispo,
                       family = "binomial")

out.null.dispo <- svyglm(as.numeric(ICF) ~  age +
                     gender + 
                     race + 
                     year + 
                     Afib + 
                     sepsis +
                     PNA + 
                     PE + 
                     DVT + 
                     Bleed + 
                     invasive_mgmt +
                     NSTEMI + 
                     severe_MI +
                     NSTEMI, 
                   design = survey_df_full_dispo,
                   family = "binomial")

model_comp_dispo <- data.frame(model = c("Full", "RCRIPlus", "null"), 
                         n_par = c(out.full$rank, out.RCRIplus$rank, out.null$rank),
                         BIC = c(BIC(out.full.dispo, maximal = out.full.dispo)["BIC"], BIC(out.RCRIplus.dispo, maximal = out.RCRIplus.dispo)["BIC"],  BIC(out.null.dispo, maximal = out.null.dispo)["BIC"])) %>% 
  mutate(delta = BIC - min(BIC))

confidence_intervals_dispo =  confint.default(out.full.dispo, adjust = "bonferroni")
mod_coefs_dispo <- data.frame(apply(exp(cbind(OR = coef(out.full.dispo),confidence_intervals_dispo)),2,round,2))
names(mod_coefs_dispo)<- c("OR", "2.5% CI", "97.5% CI")
coef_df_dispo <- data.frame(Var = rownames(mod_coefs_dispo),
                      OR = paste0(mod_coefs_dispo$OR, " [", mod_coefs_dispo$"2.5% CI", ", ",  mod_coefs_dispo$"97.5% CI", "]"))
if(output_tables){
  tab_coefs <- xtable(coef_df_dispo)
  print(tab_coefs, file="dispo_model_coefs.txt", rownames = F)
}

## PCA analysis ## ------------------------
sub.PCA <- data[,(colnames(data) %in% c("died", "year", "race", "gender", "invasive_mgmt", "PNA", "sepsis", "Afib", "PE", "DVT",  "Bleed", "NSTEMI", "severe_MI", "trendwt", unlist(retained_vars)))]
out.PCA <- glm(as.numeric(died) ~., 
               data = sub.PCA,
               family = "binomial")