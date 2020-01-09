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
select <- dplyr::select

#Plotting specs
textSize = 12
save_plots = F
source("plot_themes.R")

#Data saving
data_output_filename = "data_full_NA.rda"
save_data = T
drop_missing = F
generate_data = F

if(generate_data){
# Identify procedures that correspond to non-cardiac surgeries
year_vec = c(2008:2013)
n_procedures_vec <- rep(15,length(year_vec)) #rep(15,length(year_vec)) # Number of possible procedures per individual for this year 
n_dx_vec <- c(15, 25, 25, 25, 25, 25) # Number of possible diagnoses per individual for this year 
source("process_data.R")

data <- data_all %>% 
  rename(obesity = cm_obese,
         alcoholic = cm_alcohol,
         HTN = cm_htn_c,
         valve_dz = cm_valve,
         chrnlung = cm_chrnlung,
         anemia = cm_anemdef,
         PAD = cm_perivasc,
         liver_dz = cm_liver) %>% 
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

if(save_data){
  save(data, file = data_output_filename)
}
}

if(!generate_data){
  load(data_output_filename)
}
## Missing data ## ---------------------------------------------------------------------------------------------------------------
dfm_missing_demog <- data %>% select(c(year, race, gender, smoking, alcoholic, high_risk_surgery)) %>% 
  melt(., id.vars = "year") %>% 
  group_by(year,variable) %>% 
  summarize(Completeness = (mean(! is.na(value)) * 100),
            n = n()) %>% 
  mutate(var_type = "Demographic variables")

dfm_missing_ASCVD_risk <- data %>% 
  select(year, hx_chf, hx_DM, HTN, obesity, HLD, hx_ckd) %>%  
  melt(., id.vars = "year") %>% 
  group_by(year,variable) %>% 
  summarize(Completeness = (mean(! is.na(value)) * 100),
            n = n()) %>% 
  mutate(var_type = "ASCVD risk factors")

dfm_missing_ASCVD <- data %>% 
  select(year, hx_isch_heart, PAD, valve_dz, hx_CVA) %>%  
  melt(., id.vars = "year") %>% 
  group_by(year,variable) %>% 
  summarize(Completeness = (mean(! is.na(value)) * 100),
            n = n()) %>% 
  mutate(var_type = "Known ASCVD")

dfm_missing_other_risk <- data %>% 
  select(year, chrnlung, cm_mets, malignancy, anemia) %>%  
  melt(., id.vars = "year") %>% 
  group_by(year,variable) %>% 
  summarize(Completeness = (mean(! is.na(value)) * 100),
            n = n()) %>% 
  mutate(var_type = "Other risk factors")

df_missing_data <- rbind(dfm_missing_demog, 
                         dfm_missing_ASCVD_risk, 
                         dfm_missing_ASCVD, 
                         dfm_missing_other_risk
                         ) %>% 
  ungroup() %>% 
  mutate( year = as.numeric(as.character(year)))
  
p_missing <- ggplot(df_missing_data, aes(x = as.numeric(year), y = Completeness)) + 
  geom_line(aes(color = variable)) + 
  facet_wrap(~var_type) + 
  plot_themes + 
  ylim(50,100) + 
  xlab("")

if(save_plots){
  save_plot("Missing_data.pdf", p_missing, base_width = 6, base_height = 6)
}

## Exploratory plots ## ---------------------------------------------------------------------------------------------------------------

dfm_dat_demographic_cat <- data_all %>% select(c(year, race, gender, smoking)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

dfm_dat_demographic_cont <- data_all %>% select(c(year,age)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

dfm_smoking_age <- data_all %>% select(c(year, age_factor, smoking)) %>% 
  drop_na() %>% 
  group_by(year, age_factor) %>% 
  summarize(mean = mean((as.numeric(smoking)-1)),
            sd = sd(smoking))

p_dem_cat <- ggplot(dfm_dat_demographic_cat, aes(x = year, y = mean)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(group = variable, color = variable), linetype = 2) +
  plot_themes + 
  xlab("") + 
  ylab("Mean fraction") + 
  ggtitle("Demographic and behavioral factors")

p_dem_cont <- 
  ggplot(dfm_dat_demographic_cont, aes(x = year, y = mean)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(group = variable, color = variable), linetype = 2) +
  geom_errorbar(aes(color = variable, ymin = mean - sd, ymax = mean + sd)) + 
  plot_themes + 
  ylim(40,85) + 
  xlab("") + 
  ylab("Mean")

p_smoking_age <- ggplot(dfm_smoking_age, aes(x = year, y = mean)) + 
  geom_point(aes(color = age_factor)) + 
  geom_line(aes(group = age_factor, color = age_factor), linetype = 2) +
  #geom_errorbar(aes(color = variable, ymin = mean - sd, ymax = mean + sd)) + 
  plot_themes + 
  xlab("") + 
  ylab("Mean fraction") + 
  labs(color = "Age category") +
  ggtitle("Smoking by age group")

p_dem <- plot_grid(p_dem_cat, p_dem_cont, ncol = 2)
save_plots = F
if(save_plots){
  save_plot("demographics.pdf", p_dem, base_width = 14, base_height = 6)
  save_plot("smoking_by_age.pdf", p_smoking_age, base_width = 8, base_height = 4)
}

dfm_dat_cardiac_risk <- data_all %>% #select(c(cm_chf, cm_dm, cm_perivasc, cm_htn_c, cm_renlfail, HLD, CAD, prior_PCI, prior_CABG, year)) %>% 
  select(year, hx_chf, hx_DM, HTN, obesity, HLD, hx_ckd) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

p_cardiac_risk <- ggplot(dfm_dat_cardiac_risk, aes(x = year, y = mean)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(group = variable, color = variable), linetype = 2) +
  plot_themes + 
  xlab("") +
  labs(color = "Cardiac risk factor") +
  ylab("Mean fraction")+ 
  ggtitle("Cardiac risk factors")

dfm_dat_ASCVD <- data_all %>% #select(c(cm_chf, cm_dm, cm_perivasc, cm_htn_c, cm_renlfail, HLD, CAD, prior_PCI, prior_CABG, year)) %>% 
  select(year, hx_isch_heart, PAD, hx_CVA) %>% 
  mutate(ASCVD = as.numeric(PAD ==1 | hx_isch_heart == 1 | hx_CVA == 1)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

p_ASCVD <- ggplot(dfm_dat_ASCVD, aes(x = year, y = mean)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(group = variable, color = variable), linetype = 2) +
  plot_themes + 
  xlab("") +
  ylab("Mean fraction") + 
  labs(color = "ASCVD type") + 
  ggtitle("ASCVD")

p_cardiac <- plot_grid(p_cardiac_risk, p_ASCVD, ncol = 2)
if(save_plots){
  save_plot("cardiac_risk.pdf", p_cardiac, base_width = 14, base_height = 6)
}

data_all$invasive_mgmt = as.numeric(data_all$invasive_mgmt) -1
dfm_died <- data_all %>% select(c(year, died)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

dfm_inv <- data_all %>% select(c(year, invasive_mgmt)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

p_died <- ggplot(dfm_died, aes(x = year, y = mean)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(group = variable, color = variable), linetype = 2) +
  #geom_errorbar(aes(group = variable, color = variable, ymin = mean -sd, ymax = mean + sd)) + 
  ylim(0,0.2) + 
  plot_themes + 
  ggtitle("Mortality")

p_invasive <- ggplot(dfm_inv, aes(x = year, y = mean)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(group = variable, color = variable), linetype = 2) +
  #geom_errorbar(aes(group = variable, color = variable, ymin = mean -sd, ymax = mean + sd)) + 
  ylim(0,1) + 
  plot_themes + 
  ggtitle("Invasive management")

if(save_plots){
  save_plot("Mortality.pdf", p_died, base_width = 8, base_height = 4)
  save_plot("Invasive_mgmt.pdf", p_invasive, base_width = 8, base_height = 4)
}

dfm_outcomes_cont <- data_all %>% select(c(year, ndx,los)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

p_outcomes <- ggplot(dfm_outcomes, aes(x = year, y = mean)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(group = variable, color = variable), linetype = 2) +
  geom_errorbar(aes(group = variable, color = variable, ymin = mean -sd, ymax = mean + sd)) + 
  facet_grid(variable~., scales = "free") + 
  plot_themes 

