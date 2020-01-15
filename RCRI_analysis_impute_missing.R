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

## Exploratory plots ## ------------------------------------------------------------------------------------------------
dfm <- melt(data_all, id.vars = c("ind","year"))

p_RCRI <- ggplot(dfm %>% filter(variable == "RCRI_pt"), aes(x = value )) + 
  geom_bar(aes(fill = year), stat = "count", position = "dodge", alpha = .5) +
  xlab("RCRI") + ylab("") +
  ggtitle("RCRI distribution by year") + 
  plot_themes# + 
 # facet_grid(year~., scales = "free_y")

if(save_plots){
  save_plot("RCRI_distribution.pdf", p_RCRI, base_width = 8, base_height = 4)
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

## Logit Regression analysis ## ------------------------------------------------------------------------------

sub.PCA <- data[,(colnames(data) %in% c("died", "year", "race", "gender", "invasive_mgmt", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI", unlist(retained_vars)))]
out.PCA <- glm(as.numeric(died) ~., data = sub.PCA, family = "binomial")

var_names_full <- c(names(X1), names(X2), "cm_mets", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI")#, "nchronic")
data_full <- data %>% select(c(year, died, invasive_mgmt, var_names_full))
out.full <- glm(died ~., data = data_full, family = "binomial")

out.RCRI <- glm(as.numeric(died) ~  age +
               gender + 
               race + 
               year + 
               #smoking + 
               #HLD + 
               hx_isch_heart +
               hx_chf +
               hx_CVA + 
               hx_DM + 
               hx_ckd + 
               high_risk_surgery + 
               invasive_mgmt + 
               NSTEMI + 
               severe_MI, 
               data = data, family = "binomial")

out.RCRIplus <- glm(as.numeric(died) ~  age +
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
                      hx_chf +
                      hx_CVA + 
                      hx_DM + 
                      hx_ckd + 
                      high_risk_surgery + 
                      invasive_mgmt +
                      NSTEMI*Bleed + 
                      hx_isch_heart + 
                      severe_MI , 
                    data = data, family = "binomial")

model_comp <- data.frame(model = c("Full", "PCA", "RCRIPlus"), 
                         BIC = c(BIC(out.full), BIC(out.PCA), BIC(out.RCRIplus))) %>% 
  mutate(delta = BIC - min(BIC))

mod_coefs <- data.frame(apply(exp(cbind(OR = coef(out.RCRIplus), confint(out.RCRIplus))),2,round,2))
names(mod_coefs)<- c("OR", "2.5% CI", "97.5% CI")
if(output_tables){
  tab_coefs <- xtable(mod_coefs)
  print(tab_coefs, file="mortality_model_coefs.txt")
}


## Testing model out of sample ## -------------------------------------------

set.seed(364)
sample <- sample(nrow(data),floor(nrow(data)*0.8))
train <- data[sample,]
test <- data[-sample,]

library(pROC)

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

## Examining discharge outcomes 

#Select retained variables
out.contrib$max_contrib <- apply(out.contrib %>% select(-var),1,max) 
retained_vars <- out.contrib %>% filter(max_contrib > 9) %>% select(var)

## Logit Regression analysis ## ------------------------------------------------------------------------------
data$ICF <- as.numeric(data$ICF==1)
sub.PCA.dispo <- data[,(colnames(data) %in% c("ICF", "died", "year", "race", "gender", "invasive_mgmt", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI", unlist(retained_vars)))]
out.PCA.dispo <- glm(as.numeric(ICF) ~., data = sub.PCA.dispo, family = "binomial")
#summary(out.PCA)
#exp(cbind(OR = coef(out.PCA), confint(out.PCA)))

var_names_full <- c(names(X1), names(X2), "cm_mets", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", "NSTEMI", "severe_MI")
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

# out.RCRI <- glm(as.numeric(died) ~ age + 
#                         gender + 
#                         race + 
#                         year + 
#                         smoking +  
#                         HLD + 
#                         hx_ckd + 
#                         hx_isch_heart + 
#                         hx_chf +
#                         high_risk_surgery + 
#                         invasive_mgmt, data = data, family = "binomial")

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
sample <- sample(nrow(data),floor(nrow(data)*0.8))
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