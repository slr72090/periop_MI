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
data_filename = "Data_objects/data_raw_with_weights_imputed.rda"
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


## Step 1: retain all variables ## -----------------------------------------------------------------------------------------
dat_PCA <- data %>%
  mutate(RCRI_pt = as.numeric(scale(as.numeric(RCRI_pt))))  

dat_PCA <- dat_PCA[, (colnames(dat_PCA) %in% c("age", 
                                               "race",
                                               "gender",
                                               "obesity",
                                               "smoking",
                                               "alcoholic",
                                               "HTN",
                                               "HLD",
                                               "hx_DM",
                                               "hx_ckd",
                                               "hx_isch_heart",
                                               "hx_revasc",
                                               "PAD",
                                               "valve_dz",
                                               "hx_chf",
                                               "hx_VTE",
                                               "chrnlung",
                                               "malignancy",
                                               "cm_mets",
                                               "anemia",
                                               "hx_CVA",
                                               "high_risk_surgery"
                                               ))]

split.data.all <- splitmix(dat_PCA)

split.data <- split.data.all
X1 <- split.data$X.quanti 
X2 <- split.data$X.quali

res.pcamix <- PCAmix(X.quanti=X1, X.quali=X2,rename.level=TRUE, ndim = 5,
                     graph=FALSE)

out.eig <- as.data.frame(res.pcamix$eig) 
out.eig$dim <- c(1:nrow(out.eig))
  
out.ind <- as.data.frame(res.pcamix$ind$coord) %>% 
  mutate(ind = row_number()) %>% 
  left_join(data, by = "ind") %>% 
  drop_na()

out.quant <- res.pcamix$quanti$contrib.pct
out.qual <- res.pcamix$quali$contrib.pct
out.contrib <- rbind(out.quant,out.qual) %>% 
  as.data.frame()
out.contrib$var <- rownames(out.contrib)

out.coord <- as.data.frame(rbind(res.pcamix$quanti.cor, res.pcamix$categ.coord))
out.coord$var <- rownames(out.coord)
out.loadings <- res.pcamix$sqload

out.coord$sig <- as.factor(out.coord$`dim 1` > 0.4 | out.coord$`dim 2` > 0.5 |out.coord$`dim 1` < -0.4 | out.coord$`dim 2` < -0.5 )

## Feature selection -------------------------------------------------------------------------------------------

#Select retained variables
out.contrib$max_contrib <- apply(out.contrib %>% select(-var),1,max) 
retained_vars <- out.contrib %>% filter(max_contrib > 10) %>% select(var)

### Plot PCA output ## -----------------------------------------------------------------------------------------
xlabel <- paste0(" Dim 1, ",round(out.eig$Proportion[1],2),"%")
ylabel <- paste0(" Dim 2, ",round(out.eig$Proportion[3],2),"%")

out.coord.sig <- out.coord %>% filter(sig ==T)
p_loadings = ggplot()+
  #geom_point(data = out.loadings, aes(x=`dim 1`, y= `dim 2`, colour = "#DCDCDC"))+
  geom_segment(data = out.coord, aes(x=0,y=0,xend=`dim 1`,yend=`dim 2`), color = "grey",
               arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(data = out.coord.sig, aes(x=`dim 1`, y=`dim 2`, label=factor(out.coord.sig$var)), color = "black") + 
  xlab(xlabel) + 
  ylab(ylabel) + 
  plot_themes

p_cumvar <- ggplot(out.eig, aes(x = dim, y = Cumulative)) + 
  geom_point() + 
  geom_line(linetype = 2) +
  geom_hline(aes(yintercept = 70, color = "red")) + 
  labs("") +
  xlab("Principal component dimension") + 
  ylab("Percentage of cumulative variance") +
  plot_themes

p_eig <- ggplot(out.eig, aes(x = dim, y = Eigenvalue)) + 
  geom_point() + 
  geom_line(linetype = 2) +
  geom_hline(aes(yintercept = 1, color = "red")) + 
  labs("") +
  xlab("Principal component dimension") + 
  ylab("Eigenvalue") +
  plot_themes + 
  theme(legend.position = "none")

if(save_plots){
  save_plot("eigenvalues.pdf", p_eig, base_width =  8, base_height =  6)
  save_plot(paste0("loadings.pdf"), p_loadings, base_width = 8, base_height = 6)
}

p_ind_RCRI <-  qplot(data = out.ind, x = `dim 1`, y = `dim 2`, colour = `RCRI >= 3`) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = `RCRI >= 3`)) +
  plot_themes

p_ind_year <-  qplot(data = out.ind, x = `dim 1`, y = `dim 2`, colour = year) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = year)) +
  plot_themes

if(save_plots){
  save_plot("RCRI_ind.pdf", p_ind_RCRI, base_width = 6, base_height =4)
  save_plot("year_ind.pdf", p_ind_year, base_width = 8, base_height =6)
  ind_plots <- plot_grid(p_ind_RCRI, p_ind_year, ncol =2)
  save_plot("ind_plots.pdf", ind_plots, base_width = 12, base_height = 4)
}


## Logit Regression analysis ## ------------------------------------------------------------------------------

survey_df <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = data, nest =T)
sub.PCA <- data[,(colnames(data) %in% c("died", "year", "race", "gender", "invasive_mgmt", "PNA", "sepsis", "Afib", "PE", "DVT",  "Bleed", "NSTEMI", "severe_MI", "trendwt", unlist(retained_vars)))]
out.PCA <- glm(as.numeric(died) ~., 
               data = sub.PCA,
               weights = trendwt,
               family = "binomial")

var_names_full <- c(names(X1), names(X2), "cm_mets", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI")#, "nchronic")
data_full <- data %>% select(c(year, died, invasive_mgmt, trendwt, var_names_full))
out.full <- glm(died ~., 
                data = data_full, 
                family = "binomial")

options(survey.lonely.psu="adjust")

out.RCRIplus <- svyglm(as.numeric(died) ~  age +
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
                         hx_isch_heart + 
                         hx_revasc + 
                         high_risk_surgery + 
                         #surgery_type + 
                         invasive_mgmt +
                         hx_revasc + 
                         NSTEMI + 
                         severe_MI ,
                       design = survey_df)

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
                      hx_isch_heart + 
                      high_risk_surgery + 
                      #surgery_type + 
                      invasive_mgmt +
                      NSTEMI + 
                      hx_revasc + 
                      severe_MI , 
                    data = data, 
                    family = "binomial")

out.null <- glm(as.numeric(died) ~  age +
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
                  NSTEMI*Bleed, 
                data = data,
                weights = trendwt,
                family = "binomial")

model_comp <- data.frame(model = c("Full", "PCA", "RCRIPlus", "null"), 
                         n_par = c(out.full$rank, out.PCA$rank, out.RCRIplus$rank, out.null$rank),
                         LL = c(logLik(out.full),logLik(out.PCA), logLik(out.RCRIplus), logLik(out.null)),
                         BIC = c(BIC(out.full), BIC(out.PCA), BIC(out.RCRIplus), BIC(out.null))) %>% 
  mutate(delta = BIC - min(BIC))

confidence_intervals =  confint.default(out.RCRIplus, adjust = "bonferroni")
mod_coefs <- data.frame(apply(exp(cbind(OR = coef(out.RCRIplus),confidence_intervals)),2,round,2))
mod_coefs <- data.frame(apply(exp(cbind(OR = coef(out.RCRIplus), confint(out.RCRIplus))),2,round,2))
names(mod_coefs)<- c("OR", "2.5% CI", "97.5% CI")
if(output_tables){
  tab_coefs <- xtable(mod_coefs)
  print(tab_coefs, file="mortality_model_coefs.txt")
}


## Testing model out of sample ## -------------------------------------------

set.seed(364)
sample <- sample(nrow(data),floor(nrow(data)*0.75))
train <- data[sample,]
test <- data[-sample,]

library(pROC)

out.RCRIplus.train <- glm(as.numeric(died) ~ age +
                            gender + 
                            race + 
                            year + 
                            Afib + 
                            sepsis +
                            PNA + 
                            PE + 
                            DVT + 
                            Bleed + 
                            hx_isch_heart +
                            PAD + 
                            hx_chf +
                            hx_CVA + 
                            hx_DM + 
                            hx_ckd + 
                            high_risk_surgery + 
                            invasive_mgmt +
                            hx_revasc + 
                            NSTEMI + 
                            severe_MI, 
                          data = train, 
                          family = "binomial")

out.full.train <- glm(died ~., data = (train %>% select(c(year, died, invasive_mgmt, severe_MI, var_names_full))), family = "binomial")

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

if(compare_models){
  var_names_RCRI <- c("age", "race", "gender", 
                      "hx_DM", "Afib", "sepsis", 
                      "PNA", "PE", "DVT",
                      "Bleed", "hx_chf", "hx_CVA", 
                      "hx_DM", "hx_ckd", "hx_isch_heart", 
                      "hx_revasc", "high_risk_surgery", 
                      "invasive_mgmt", "NSTEMI", "severe_MI")
  dat_AUC <- data %>% filter(!is.na(died)) %>% 
    select(c(died,var_names_full)) %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate_if(is.character, as.integer)
  
  mat_RCRI <- as.matrix(dat_AUC[,(names(dat_AUC) %in% var_names_RCRI)])
  mat_Full <- as.matrix(dat_AUC[,var_names_full[!(var_names_full %in% var_names_RCRI)]])[,-13]
  
  AUC_compare <- deltaAUC(y = dat_AUC$died, 
                          x = mat_RCRI,
                          z = mat_Full)
}                  
                        
                      #  x = dat_AUC[,(names(dat_AUC) %in% var_names_RCRI)], 
                      #  z = dat_AUC[,var_names_full[!(var_names_full %in% var_names_RCRI)]])
ROC_test_size <- data.frame()
test_percent_vec <- seq(0.5,1.0,0.05)
for(i in c(1:length(test_percent_vec))){
  print(i)
  set.seed(364)
  sample <- sample(nrow(data),floor(nrow(data)*test_percent_vec[i]))
  this_train <- data[sample,]
  this_test <- data[-sample,]
  out.RCRIplus.train <- glm(as.numeric(died) ~ age +
                              gender + 
                              race + 
                              year + 
                              #smoking + 
                              #HLD + 
                              Afib + 
                              sepsis +
                              PNA + 
                              PE + 
                              DVT + 
                              Bleed + 
                              hx_isch_heart +
                              PAD + 
                              hx_chf +
                              hx_CVA + 
                              hx_DM + 
                              hx_ckd + 
                              high_risk_surgery + 
                              invasive_mgmt +
                              NSTEMI + 
                              severe_MI, 
                            data = this_train, 
                            family = "binomial")
  test_prob_RCRIplus_temp = predict(out.RCRIplus.train, newdata = this_test, type = "response")
  test_roc_RCRIplus_temp = roc(this_test$died ~ test_prob_RCRIplus_temp, plot = TRUE, print.auc = TRUE)
  this_df <- data.frame(percent = test_percent_vec[i],
                        AUC = test_roc_RCRIplus_temp$auc
                        )
  ROC_test_size <- rbind(ROC_test_size, this_df)
  
}

## Examining discharge outcomes 

#Select retained variables
out.contrib$max_contrib <- apply(out.contrib %>% select(-var),1,max) 
retained_vars <- out.contrib %>% filter(max_contrib > 9) %>% select(var)

## Prediction by number of RCRI factors ## --------------------------------------
dat_sub <- data %>% filter(!is.na(died) #& 
                           #invasive_mgmt == 0 & 
                           #sepsis == 0 & 
                           #PNA == 0  
                           #Bleed == 0 
                           )

pred.test <- predict(out.RCRIplus, 
                      dat_sub,
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

p_risk <- ggplot(df, aes(x = RCRI, y = med)) +  
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = q1, ymax = q2), linetype = 2) + 
  facet_wrap(.~invasive_mgmt + hx_isch_heart, scales = "free_y") +
  xlab("RCRI") + 
  ylab("Median probability (IQR)") +
  plot_themes
  
if(save_plots){
  save_plot("Risk_by_RCRI_and_CAD.pdf", p_risk, base_width = 8, base_height = 6)
}


## Analysis of dispo outcome ## ------------------------------------------------------------------------------
data$ICF <- as.numeric(data$ICF==1)
sub.PCA.dispo <- data[,(colnames(data) %in% c("ICF", "died", "year", "race", "gender", "invasive_mgmt", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI", unlist(retained_vars)))]
out.PCA.dispo <- glm(as.numeric(ICF) ~., data = sub.PCA.dispo, family = "binomial")
#summary(out.PCA)
#exp(cbind(OR = coef(out.PCA), confint(out.PCA)))

var_names_full <- c(names(X1), names(X2), "cm_mets", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI","invasive_mgmt", "severe_MI")
data_full_dispo <- data %>% select(c(year, died, ICF, invasive_mgmt, var_names_full)) # %>% select(-c(RCRI_pt, Ischemic_stroke, los, ind, age_factor, `RCRI > 3`))
out.full.dispo <- glm(ICF ~. + gender*NSTEMI*hx_revasc , data = data_full_dispo, family = "binomial")

out.RCRIplus.dispo <- glm(as.numeric(ICF) ~  age +
                  gender + 
                  race + 
                  year + 
                  hx_isch_heart +
                  hx_chf +
                  hx_CVA + 
                  hx_DM + 
                  hx_ckd + 
                  high_risk_surgery + 
                  invasive_mgmt +
                  NSTEMI + 
                    severe_MI + 
                    PNA + 
                    sepsis + 
                    Afib + 
                    PE + 
                    DVT + 
                    Bleed +
                    died, 
                data = data, family = "binomial")

model_comp.dispo <- data.frame(model = c("Full", "PCA",  "RCRIPlus"), 
                         BIC = c(BIC(out.full.dispo), 
                                 BIC(out.PCA.dispo), 
                                 BIC(out.RCRIplus.dispo))) %>% 
  mutate(delta = BIC - min(BIC))
mod_coefs.dispo <- data.frame(apply(exp(cbind(OR = coef(out.full.dispo), confint(out.full.dispo))),2,round,2))
names(mod_coefs.dispo)<- c("OR", "2.5% CI", "97.5% CI")
if(output_tables){
  tab_coefs <- xtable(mod_coefs.dispo)
  print(tab_coefs, file="dispo_model_coefs.txt")
}


## Testing model out of sample ## -------------------------------------------
set.seed(123)
sample <- sample(nrow(data),floor(nrow(data)*0.75))
train.dispo <- data[sample,]
test.dispo <- data[-sample,]

library(pROC)

out.RCRIplus.train.dispo <- glm(as.numeric(ICF) ~ age +
                                  gender + 
                                  race + 
                                  year + 
                                  #smoking + 
                                  #HLD + 
                                  Afib + 
                                  sepsis +
                                  PNA + 
                                  PE + 
                                  DVT + 
                                  Bleed + 
                                  hx_isch_heart +
                                  PAD + 
                                  hx_chf +
                                  hx_CVA + 
                                  hx_DM + 
                                  hx_ckd + 
                                  high_risk_surgery + 
                                  invasive_mgmt +
                                  NSTEMI + 
                                  severe_MI +
                                  died, 
                                data = train, 
                                family = "binomial")

out.full.train.dispo <- glm(ICF ~., data = (train %>% select(c(year, died, invasive_mgmt, severe_MI, var_names_full))), family = "binomial")

test_prob_RCRIplus.dispo = predict(out.RCRIplus.train.dispo, newdata = test, type = "response")
test_roc_RCRIplus.dispo = roc(test$ICF~ test_prob_RCRIplus, plot = TRUE, print.auc = TRUE)

test_prob_full.dispo = predict(out.full.train.dispo, newdata = test, type = "response")
test_roc_full.dispo = roc(test$ICF ~ test_prob_full.dispo, plot = TRUE, print.auc = TRUE)

p_roc <- ggroc(list(`Full model` = test_roc_full.dispo, `RCRI model` = test_roc_RCRIplus.dispo), linetype = 2) +
  plot_themes + 
  labs(color = "Model") + 
  ggtitle("ROC curves")

if(save_plots){
  save_plot("ROC_curves_dispo.pdf", p_roc, base_width = 8, base_height = 4)
}