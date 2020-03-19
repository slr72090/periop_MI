## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

# Set the working directory to the source file directory
setwd("~/Desktop/perip_MI_GITHUB/predict_MI_model")

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
select <- dplyr::select

#Plotting specs
textSize = 12
save_plots = F
source("../utility_scripts/plot_themes.R")

#Data saving
dbfilename = "../data_files/NIS.db"
data_filename = "data_raw_with_surgery_type.rda"
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
           #age = as.numeric(scale(age)),
           #nchronic = as.numeric(scale(nchronic)),
           #invasive_mgmt = as.factor(invasive_mgmt),
           high_risk_surgery = as.factor(as.numeric(transplant == 1|thoracic_surgery == 1|vascular == 1|abdominal ==1)),
           hx_revasc = as.factor(as.numeric(prior_CABG == 1 | prior_PCI ==1))) %>% 
    mutate(RCRI_pt = as.factor(as.numeric(RCRI_pt) + as.numeric(high_risk_surgery == 1))) %>% 
    select(-c(prior_CABG, prior_PCI, prior_MI, CAD, transplant,thoracic_surgery,vascular,nchronic)) %>% 
    mutate(`RCRI >= 3` = as.factor(as.numeric(RCRI_pt) > 3))
  
  data <- data %>% rename(high_risk_surgery_2 = high_risk_surgery) %>% 
    mutate(high_risk_surgery = as.factor(as.numeric(high_risk_surgery_2 == 1|abdominal ==1)))
  
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
   
  data_formatted <- data %>% 
    rename(obesity = cm_obese,
           alcoholic = cm_alcohol,
           HTN = cm_htn_c,
           valve_dz = cm_valve,
           chrnlung = cm_chrnlung,
           anemia = cm_anemdef,
           PAD = cm_perivasc) %>% 
    mutate(age_factor = as.factor(ntile(age,3)),
           #age = as.numeric(scale(age)),
           #nchronic = as.numeric(scale(nchronic)),
           #invasive_mgmt = as.factor(invasive_mgmt),
           high_risk_surgery = as.factor(as.numeric(transplant == 1|thoracic_surgery == 1|vascular == 1)),
           hx_revasc = as.factor(as.numeric(prior_CABG == 1 | prior_PCI ==1))) %>% 
    mutate(RCRI_pt = as.factor(as.numeric(RCRI_pt) + as.numeric(high_risk_surgery == 1))) %>% 
    select(-c(prior_CABG, prior_PCI, prior_MI, CAD, transplant,thoracic_surgery,vascular,nchronic)) %>% 
    mutate(`RCRI >= 3` = as.factor(as.numeric(RCRI_pt) > 3))
  
  data <- data_formatted %>% 
    mutate(MI = as.numeric(MI==1))
  rm(data_formatted)
  rm(data_all)
}
# 
# 
# ## Determine surgery type ## -----------------------------------------
# pr_list <- list(c(1:9),c(10:12),
#                 c(13:21), 
#                 c(22:35), 
#                 c(36:42),
#                 c(43:50, 62:63),
#                 c(51:61),
#                 c(66:99), 
#                 c(100:118),
#                 c(119:125,129:132),
#                 c(126:128, 133:141),
#                 c(142:164),
#                 c(165:167),
#                 c(168:175),
#                 176)
# names(pr_list) = c("Neuro", "Endo","Optho","ENT","Thoracic","Cardiac","Vascular","General","GU", "Gynecologic", "Obstetric", "Orthopedic","Breast","Skin/burn", "Transplant")
# 
# pr_df <- stack(pr_list) %>% 
#   as.data.frame() %>% 
#   stats::setNames(c("code","type"))
# 
# data$surgery_type <- sapply(as.numeric(core_df$prccs1), pr_class_fun, code_df = pr_df)

## MI statistics ## ---------------------------------------
n_MI <- data %>% group_by(year) %>% 
  drop_na() %>% 
  summarize(n = sum(MI),
            percent = mean(MI)) %>% 
  mutate(year = as.character(year)) %>% 
    mutate(lab = paste0("n = ",n))

p_MI <- ggplot(n_MI, aes(x= as.numeric(year), y = percent)) +
  geom_point() + 
  geom_line(linetype = 2) + 
  plot_themes + 
  ylim(0,0.01) + geom_text(aes(label = lab, vjust = -1.5)) +
  xlab("") + ylab("Yearly incidence of AMI")

if(save_plots){
  save_plot("n_MI.pdf", p_MI, base_width =8 , base_height = 4 )
}

n_MI["Total",] = c("Aggregate", sum(data$MI, na.rm=T),mean(data$MI, na.rm=T), paste0("n = ", sum(data$MI, na.rm=T)))


dat_by_MI <- 
  data %>% select(-c(contains("RCRI"), died, ICF, year,age,ind, sepsis, PNA, PE, DVT, Bleed, NSTEMI, died, invasive_mgmt, IABP, cardiogenic_shock, ICF, age_factor, surgery_type, high_risk_surgery_2,prccs1))  %>%
  mutate_if(is.factor, as.character) %>% 
  mutate_if(is.character, as.numeric) %>% 
  drop_na() %>% 
  group_by(MI) %>% 
  summarise_all(.,list(~mean(.))) %>% 
  select(-MI) %>% 
  t() %>% 
  as.data.frame()

dat_surg_type <- data %>% 
  filter(MI ==1 ) %>% 
  filter(!is.na(surgery_type) &!is.na(died)) %>% 
  group_by(died,as.factor(surgery_type)) %>% 
  summarise(n=n())%>% 
  ungroup() %>% 
  as.data.frame() %>% 
  group_by(died) %>% 
  mutate(prop = round(n/sum(n),3)) %>% 
  as.data.frame()
  

names(dat_by_MI) <- c("no MI", "MI")
dat_by_MI$chisq.pval <- NA #array(NA,nrow(dat_by_mortality))
for(i in c(1:nrow(dat_by_MI))){
  print(i)
  name <- rownames(dat_by_MI[i,])
  dat_sub <- data %>% select(c(MI,name)) %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate_if(is.character, as.numeric) %>% 
    drop_na()
  dat_by_MI[i,]$chisq.pval <- chisq.test(table(dat_sub))$p.value
}
  
dat_by_MI <- dat_by_MI%>% 
  apply(.,2, round,4)

if(output_tables){
  n_MI_tab <- xtable(n_MI)
  dat_by_MI_tab <- xtable(dat_by_MI)
  print(n_MI_tab, file="n_MI.txt")
  print(dat_by_MI_tab, file="dat_by_MI.txt")
}

## Characterize Missing data ## ---------------------------------------------------------------------------------------------------------------
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
  
p_missing <- ggplot(dfm_missing_data, aes(x = as.numeric(as.character(year)), y = (100-Completeness))) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(color = variable), linetype =2) + 
  #facet_wrap(~var_type) + 
  plot_themes + 
  ylim(0,50) + 
  xlab("") + 
  ylab("Percent missing values")

if(save_plots){
  save_plot("Missing_data_demog_all.pdf", p_missing, base_width = 8, base_height = 6)
}

## Characterize Outcomes ## --------------------------------------------------------------------------------
dfm_outcomes <- data %>% select(year,hx_isch_heart, hx_revasc) %>% 
  apply(.,2,as.numeric) %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  summarize(
            hx_isch_heart = mean(hx_isch_heart, na.rm = T),
            hx_revasc = mean(hx_revasc, na.rm = T)
  ) %>% 
  ungroup()

trend_isch <- data %>% select(year, hx_isch_heart) %>% 
  table() %>% 
  CochranArmitageTest(., alternative = c("increasing"))

trend_revasc <- data %>% select(year, hx_revasc) %>% 
  table() %>% 
  CochranArmitageTest(., alternative = c("increasing"))
                      


                    
## Exploratory plots ## ---------------------------------------------------------------------------------------------------------------

dfm_dat_demographic_cat <- data %>% select(c(year, race, gender, smoking)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

dfm_dat_demographic_cont <- data %>% select(c(year,age)) %>% 
  melt(., id.vars = "year") %>% 
  drop_na() %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(year, variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

dfm_smoking_age <- data %>% select(c(year, age_factor, smoking)) %>% 
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

if(save_plots){
  save_plot("demographics.pdf", p_dem, base_width = 14, base_height = 6)
  save_plot("smoking_by_age.pdf", p_smoking_age, base_width = 8, base_height = 4)
}

dfm_dat_cardiac_risk <- data %>% #select(c(cm_chf, cm_dm, cm_perivasc, cm_htn_c, cm_renlfail, HLD, CAD, prior_PCI, prior_CABG, year)) %>% 
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

dfm_dat_ASCVD <- data %>% #select(c(cm_chf, cm_dm, cm_perivasc, cm_htn_c, cm_renlfail, HLD, CAD, prior_PCI, prior_CABG, year)) %>% 
  select(year, hx_isch_heart, PAD, hx_CVA) %>% 
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



