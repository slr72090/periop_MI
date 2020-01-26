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
data_filename = "./predict_MI_model/data_predict_MI_imputed.rda" 
save_data = F
drop_missing = F
impute_missing = T
generate_data = F
output_tables = F
generate_PCA = F
save_models = F

if(generate_data == T){
  # Identify procedures that correspond to non-cardiac surgeries
  year_vec = c(2008:2013)
  n_procedures_vec <- rep(15,length(year_vec)) #rep(15,length(year_vec)) # Number of possible procedures per individual for this year 
  n_dx_vec <- c(15, 25, 25, 25, 25, 25) # Number of possible diagnoses per individual for this year 
  source("process_data_predict_MI.R")
  
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
    select(-c(prior_CABG, prior_PCI, prior_MI, CAD, nchronic)) %>% #transplant,thoracic_surgery,vascular,)) %>% 
    mutate(`RCRI >= 3` = as.factor(as.numeric(RCRI_pt) > 3))
  
  if(impute_missing){
    source("impute_missing.R")
    data_formatted[,missing_vars] <- complete_data[,missing_vars]
  }
  
  data = data_formatted
  
  if(save_data){
    save(data_formatted, file = "data_predict_MI_imputed.rda")
    save(data_all, file = "data_predict_MI_raw.rda")
  }
}
 
if(generate_data == F){
  load(data_filename)
   data <- data_formatted 
   rm(data_formatted)
   data$MI <- as.numeric(data$MI == 1)
}

## Step 1: PCA ## -----------------------------------------------------------------------------------------

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

if(generate_PCA){
  res.pcamix <- PCAmix(X.quanti=X1, X.quali=X2,rename.level=TRUE, ndim = 5,
                     graph=FALSE)

  save(res.pcamix, split.data, X1, X2, file = "PCA_output_full_predict_MI.rda")
}
if(!generate_PCA){
  
}

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
  geom_segment(data = out.coord, aes(x=0,y=0,xend=`dim 1`,yend=`dim 2`), color = "grey",
               arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(data = out.coord.sig, aes(x=`dim 1`, y=`dim 2`, label=factor(out.coord.sig$var)), color = "black") + 
  xlab(xlabel) + 
  ylab(ylabel) + 
  plot_themes

if(save_plots){
  save_plot(paste0("loadings.pdf"), p_loadings, base_width = 8, base_height = 6)
}

out.ind.plot <- out.ind %>% 
  sample_n(.,20000)

p_ind_RCRI <-  qplot(data = out.ind.plot, x = `dim 1`, y = `dim 2`, colour = `RCRI >= 3`) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = `RCRI >= 3`)) +
  plot_themes

p_ind_RCRI <-  ggplot(data = out.ind.plot, aes(x = `dim 1`, y = `dim 2`, colour = `RCRI >= 3`)) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = `RCRI >= 3`)) +
  plot_themes

p_ind_MI <-  ggplot(data = out.ind.plot, aes(x = `dim 1`, y = `dim 2`, colour = MI)) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = MI)) +
  plot_themes

p_ind_MI <-  qplot(data = out.ind.plot, x = `dim 1`, y = `dim 2`, colour = MI) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = MI)) +
  plot_themes


# p_ind_year <-  qplot(data = out.ind, x = `dim 1`, y = `dim 2`, colour = year) +
#   stat_ellipse(geom = "polygon", alpha = .2, aes(fill = year)) +
#   plot_themes

p_ind_year <-  ggplot(data = out.ind.plot, aes( x = `dim 1`, y = `dim 2`, colour = year)) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = year)) +
  plot_themes

if(save_plots){
  ind_plots <- plot_grid(p_ind_RCRI, p_ind_year, ncol =2)
  save_plot("ind_plots.pdf", ind_plots, base_width = 12, base_height = 4)
}


## Logit Regression analysis ## ------------------------------------------------------------------------------
sub.PCA <- data[,(colnames(data) %in% c("MI", "year", "race", "gender", "Afib", unlist(retained_vars)))]
#c("MI", "year", "race", "gender", "PNA", "sepsis", "Afib", "PE", "DVT", "Bleed", unlist(retained_vars))
out.PCA <- glm(as.numeric(MI) ~., data = sub.PCA, family = "binomial")

var_names_full <- c(names(X1), names(X2), "Afib")
                    #"PNA", "sepsis", "Afib", "PE", "DVT", "Bleed","cm_mets")
data_full <- data %>% select(c(year, MI, var_names_full))
out.full <- glm(MI ~., data = data_full, family = "binomial")

out.RCRIplus <- glm(as.numeric(MI) ~  age +
                      gender + 
                      race + 
                      year + 
                      Afib + 
                      #sepsis +
                      #PNA + 
                      #PE + 
                      #DVT + 
                      #Bleed + 
                      hx_chf +
                      hx_CVA + 
                      hx_DM + 
                      hx_ckd +  
                      high_risk_surgery + 
                      hx_isch_heart + 
                      hx_revasc, 
                    data = data, family = "binomial")

model_comp <- data.frame(model = c("Full", "PCA", "RCRIPlus"), 
                         BIC = c(BIC(out.full), BIC(out.PCA), BIC(out.RCRIplus))) %>% 
  mutate(delta = BIC - min(BIC))

if(save_models){
  save(out.full, out.PCA, out.RCRIplus, file = "models_predict_MI.rda")
}

mod_coefs <- data.frame(apply(exp(cbind(OR = coef(out.full), confint(out.full))),2,round,2))
names(mod_coefs)<- c("OR", "2.5% CI", "97.5% CI")
if(output_tables){
  tab_coefs <- xtable(mod_coefs)
  print(tab_coefs, file="MI_model_coefs.txt")
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

