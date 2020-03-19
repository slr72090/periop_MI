# Survey weights analysis 

# Set the working directory to the source file directory
setwd("~/Desktop/perip_MI_GITHUB/mortality_dispo_models")

## Load package scripts -----------------------------------------------
library(survey)
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
library(tableone)
select <- dplyr::select
options(survey.lonely.psu="adjust")

#Plotting specs
textSize = 12
save_plots = F
source("../utility_scripts/plot_themes.R")

#Data saving
data_filename = "data_all_with_NA_weighted.rda"
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
           invasive_mgmt = as.factor(invasive_mgmt),
           high_risk_surgery = as.factor(as.numeric(transplant == 1|thoracic_surgery == 1|vascular == 1|abdominal == 1)),
           high_risk_surgery_2 = as.factor(as.numeric(transplant == 1|thoracic_surgery == 1|vascular == 1)),
           hx_revasc = as.factor(as.numeric(prior_CABG == 1 | prior_PCI ==1))) %>% 
    mutate(RCRI_pt = as.factor(as.numeric(RCRI_pt) + as.numeric(high_risk_surgery == 1))) %>% 
    mutate(severe_MI = as.factor(as.numeric(cardiogenic_shock == 1 | IABP == 1))) %>% 
    select(-c(prior_CABG, prior_PCI, prior_MI, CAD)) %>% #transplant,thoracic_surgery,vascular,)) %>% 
    mutate(`RCRI >= 3` = as.factor(as.numeric(RCRI_pt) > 3))
  
  save(data_formatted, data_all, file = data_filename)
  
  if(impute_missing){
    source("impute_missing.R")
    data_formatted[,missing_vars] <- complete_data[,missing_vars]
    save(data_formatted, data_all, file = data_filename)
  }
  
  data = data_formatted
  
  if(save_data){
    save(data_formatted, data_all, file = data_filename)
  }
}
if(generate_data == F){
  load(data_filename)
}
# Classify surgery type 
## Define procedure type ## --------------------------------------
pr_list <- list(c(1:9),
                c(10:12),
                c(13:21), 
                c(22:35), 
                c(36:42),
                c(43:50, 62:63),
                c(51:61),
                c(66:99), 
                c(100:118),
                c(119:125,129:132),
                c(126:128, 133:141),
                c(142:164),
                c(165:167),
                c(168:175),
                176)
names(pr_list) = c("Neuro", "Endo","Optho","ENT","Thoracic","Cardiac","Vascular","General","GU", "Gynecologic", "Obstetric", "Orthopedic","Breast","Skin/burn", "Transplant")

pr_class_fun <- function(x,code_df){
  if(x %in% code_df$code){
    pr_type <- code_df %>% filter(code == x) %>% 
      select(type) %>% 
      unlist()
  }
  if(!(x %in% code_df$code)){
    pr_type = NA
  }
  return(as.character(pr_type))
}

pr_df <- stack(pr_list) %>% 
  as.data.frame() %>% 
  stats::setNames(c("code","type"))

data_all$surgery_type <- sapply(as.numeric(as.character(data_all$prccs1)), pr_class_fun, code_df = pr_df)

# Specify the sampling design 
dat_sv <- data_formatted %>% mutate(row = rep(1,nrow(data_all))) 
dat_MI <- data_formatted %>% filter(MI == 1) %>% mutate(row = 1)
survey_df_total <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = dat_sv, nest =T)
survey_df_MI <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = dat_MI, nest =T)

sv_total_row <- svytotal(~row, design = survey_df_total, survey.lonely.psu="adjust", fpc = fpc)
sv_total_MI <- svytotal(~MI, design = survey_df_MI, survey.lonely.psu="adjust", fpc = fpc)
sv_total_MI_died <- svytotal(~died, design = survey_df_MI,na.rm=T,  survey.adjust.domain.lonely=TRUE, fpc = fpc)
sv_total_MI_ICF <- svytotal(~ICF, design = survey_df_MI, na.rm = T, survey.adjust.domain.lonely=TRUE, fpc = fpc)


## Make Table 1 ## -----------------------------

tab1_MI <- svyCreateTableOne(vars = c("race","age", "gender","obesity", "smoking", "alcoholic", "nchronic",
                                         "HTN", "HLD", "hx_DM", "hx_ckd", "hx_isch_heart", "hx_revasc", "PAD",
                                         "valve_dz", "hx_chf", "hx_VTE", "chrnlung", "malignancy","anemia", "cm_mets",
                                         "hx_CVA", "Afib", "sepsis", "PNA", "PE", "DVT", "Bleed", "NSTEMI", "high_risk_surgery", "surgery_type"),
                                strata = "died", data = survey_df_MI,
                                factorVars = c("race","surgery_type"))

tab1_MI_print <- xtable(rbind(print(tab1_MI$CatTable), print(tab1_MI$ContTable)))

print(tab1_MI_print, file = "tab1_MI.txt", include.rownames = T)

tab1_total <-  svyCreateTableOne(vars = c("race","age","gender","obesity", "smoking", "alcoholic", "nchronic",
                                          "HTN", "HLD", "hx_DM", "hx_ckd", "hx_isch_heart", "hx_revasc", "PAD",
                                          "valve_dz", "hx_chf", "hx_VTE", "chrnlung", "malignancy","anemia", "cm_mets",
                                          "hx_CVA", "Afib", "sepsis", "PNA", "PE", "DVT", "Bleed", "NSTEMI", "high_risk_surgery", "surgery_type"),
                                 strata = "MI", data = survey_df_total,
                                 factorVars = c("race", "surgery_type"))

tab1_total_print <- xtable(rbind(print(tab1_total$CatTable), print(tab1_total$ContTable)))

if(save_tables){
  save(tab1_total, tab1_MI, file = "Table1.rda")
}

load("Table1.rda")
tab1_MI_print <- xtable(rbind(print(tab1_MI$CatTable), print(tab1_MI$ContTable)))
table_1 <- cbind(tab1_total_print, tab1_MI_print)
tab_1 <- table_1[,!(names(table_1) == "test")]
print(xtable(tab_1), file = "table_1.txt", include.rownames = T)


## Plotting ### ------------------------------------------------------------------------------------------------------------------
dfm_missing_demog <- data_formatted %>% select(c(year, race, gender)) %>% 
  melt(., id.vars = "year") %>% 
  group_by(year,variable) %>% 
  summarize(Completeness = (mean(! is.na(value)) * 100),
            n = n()) %>% 
  mutate(var_type = "Demographic variables")

df_missing_data <- dfm_missing_demog %>% 
  ungroup() %>% 
  mutate( year = as.numeric(as.character(year)))

p_missing <- ggplot(dfm_missing_demog, aes(x = as.numeric(as.character(year)), y = (Completeness))) + 
  geom_point() + 
  geom_line(linetype =2) + 
  facet_wrap(~variable, scales = "free_y") + 
  plot_themes + 
  ylim(70,100) + 
  xlab("") + 
  ylab("Data completeness (%)")

if(save_plots){
  save_plot("Missing_data_demog_all.pdf", p_missing, base_width = 14, base_height = 6)
}


## Trends ### ------------------------------------------------------------------------------------------------------------------

dat_cardiac <- data_formatted %>% filter(MI == 1) %>% select(c(hospid, nis_stratum, trendwt, year, hx_chf, hx_DM, HTN, obesity, HLD, hx_ckd, NSTEMI, invasive_mgmt, hx_revasc, died))
dat_cardiac$HTN <- as.numeric(dat_cardiac$HTN == 1)
dat_cardiac$HLD <- as.numeric(dat_cardiac$HLD == 1)
dat_cardiac$hx_chf <- as.numeric(dat_cardiac$hx_chf == 1)
dat_cardiac$hx_DM <- as.numeric(dat_cardiac$hx_DM == 1)
dat_cardiac$hx_ckd <- as.numeric(dat_cardiac$hx_ckd == 1)
dat_cardiac$obesity <- as.numeric(dat_cardiac$obesity == 1)
dat_cardiac$hx_revasc <- as.numeric(dat_cardiac$hx_revasc == 1)
dat_cardiac$NSTEMI <- as.numeric(dat_cardiac$NSTEMI == 1)
dat_cardiac$invasive_mgmt <- as.numeric(dat_cardiac$invasive_mgmt == 1)
dat_cardiac$died <- as.numeric(dat_cardiac$died == 1)



survey_cardiac_risk <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = dat_cardiac, nest =T)
dat_chf <- svyby(~hx_chf, ~year, survey_cardiac_risk, svyciprop) 
df_chf <- data.frame(year = c(2008:2013), 
                     var = "Chronic heart failure",
                     mean = dat_chf %>% select(starts_with("hx_")) %>% unlist,
                     SE = dat_chf %>% select(starts_with("se."))%>% unlist
                     )

dat_DM <- svyby(~hx_DM, ~year, survey_cardiac_risk, svyciprop)
df_DM <- data.frame(year = c(2008:2013), 
                     var = "Diabetes mellitus",
                     mean = dat_DM %>% select(starts_with("hx_")) %>% unlist,
                     SE = dat_DM %>% select(starts_with("se."))%>% unlist)


dat_ckd <- svyby(~hx_ckd, ~year, survey_cardiac_risk, svyciprop)
df_ckd <- data.frame(year = c(2008:2013), 
                     var = "Chronic kidney disease",
                     mean = dat_ckd %>% select(starts_with("hx_")) %>% unlist,
                     SE = dat_ckd %>% select(starts_with("se.")) %>% unlist)


dat_HTN <- svyby(~HTN, ~year, survey_cardiac_risk, svyciprop) 
df_HTN <- data.frame(year = c(2008:2013), 
                     var = "Hypertension",
                     mean = dat_HTN %>% select(starts_with("HTN")) %>% unlist,
                     SE = dat_HTN %>% select(starts_with("se.")) %>% unlist
)

dat_HLD <- svyby( ~HLD,~year, survey_cardiac_risk, svyciprop) 
df_HLD <- data.frame(year = c(2008:2013), 
                     var = "Hyperlipidemia",
                     mean = dat_HLD %>% select(starts_with("HLD")) %>% unlist,
                     SE = dat_HLD %>% select(starts_with("se.")) %>% unlist
)

dat_obesity <- svyby( ~obesity, ~year, survey_cardiac_risk, svyciprop) 
df_obesity <- data.frame(year = c(2008:2013), 
                     var = "Obesity",
                     mean = dat_obesity %>% select(starts_with("obesity")) %>% unlist,
                     SE = dat_obesity %>% select(starts_with("se.")) %>% unlist
)

df_cardiac_risk <- rbind(df_chf, df_DM,df_ckd, df_HTN, df_HLD, df_obesity)
p_cardiac_risk <- ggplot(df_cardiac_risk, aes(x = year, y = mean)) + 
  geom_point(aes(color = var)) + 
  geom_line(aes(group = var, color = var), linetype = 2) +
  plot_themes + 
  xlab("") +
  labs(color = "Cardiac risk factor") +
  ylab("Mean fraction")+  
  ggtitle("Cardiac risk factors") 

dat_ASCVD <- data_formatted %>% filter(MI == 1) %>% select(c(hospid, nis_stratum, trendwt, year,PAD, hx_isch_heart, hx_CVA))
dat_ASCVD$PAD <- as.numeric(dat_ASCVD$PAD == 1)
dat_ASCVD$hx_isch_heart <- as.numeric(dat_ASCVD$hx_isch_heart == 1)
dat_ASCVD$hx_CVA <- as.numeric(dat_ASCVD$hx_CVA == 1)

survey_ASCVD <- svydesign(ids = ~hospid, strata = ~nis_stratum, weights = ~trendwt, data = dat_ASCVD, nest =T)
dat_PAD <- svyby(~PAD, ~year, survey_ASCVD, svyciprop) 
df_PAD <- data.frame(year = c(2008:2013), 
                     var = "Peripheral artery disease",
                     mean = dat_PAD %>% select(starts_with("PAD")) %>% unlist,
                     SE = dat_chf %>% select(starts_with("se."))%>% unlist
)

dat_CAD <- svyby(~hx_isch_heart, ~year, survey_ASCVD, svyciprop) 
df_CAD <- data.frame(year = c(2008:2013), 
                     var = "Coronary artery disease",
                     mean = dat_CAD %>% select(starts_with("hx_")) %>% unlist,
                     SE = dat_CAD %>% select(starts_with("se."))%>% unlist
)

dat_CVA <- svyby(~hx_CVA, ~year, survey_ASCVD, svyciprop) 
df_CVA <- data.frame(year = c(2008:2013), 
                     var = "Prior cerebrovascular accident",
                     mean = dat_CVA %>% select(starts_with("hx_")) %>% unlist,
                     SE = dat_CVA %>% select(starts_with("se."))%>% unlist
)
df_ASCVD <- rbind(df_PAD, df_CAD,df_CVA)

p_ASCVD <- ggplot(df_ASCVD, aes(x = year, y = mean)) + 
  geom_point(aes(color = var)) + 
  geom_line(aes(group = var, color = var), linetype = 2) +
  plot_themes + 
  xlab("") +
  labs(color = "Condition") +
  ylab("Mean fraction")+  
  ggtitle("ASCVD") 

if(save_plots){
  save_plot("cardiac_risk.pdf", p_cardiac_risk, base_width = 14, base_height = 6)
  save_plot("ASCVD.pdf", p_ASCVD, base_width = 14, base_height = 6)
  p_cardiac_combined <- plot_grid(p_cardiac_risk, p_ASCVD, ncol = 2)
  save_plot("cardiac_risk.pdf", p_cardiac_combined, base_width = 16, base_height = 6)
}

dat_died <- svyby(~died, ~year, survey_cardiac_risk, svyciprop) 
df_died <- data.frame(year = c(2008:2013), 
                     var = "Mortality rate",
                     mean = dat_died %>% select(starts_with("died")) %>% unlist,
                     SE = dat_died %>% select(starts_with("se."))%>% unlist
)

dat_NSTEMI <- svyby(~NSTEMI, ~year, survey_cardiac_risk, svyciprop) 
df_NSTEMI <- data.frame(year = c(2008:2013), 
                      var = "NSTEMI",
                      mean = dat_NSTEMI %>% select(starts_with("NSTEMI")) %>% unlist,
                      SE = dat_NSTEMI %>% select(starts_with("se."))%>% unlist
)

dat_inv <- svyby(~invasive_mgmt, ~year, survey_cardiac_risk, svyciprop) 
df_inv <- data.frame(year = c(2008:2013), 
                        var = "Invasive management",
                        mean = dat_inv %>% select(starts_with("inv")) %>% unlist,
                        SE = dat_inv %>% select(starts_with("se."))%>% unlist
)

dat_revasc <- svyby(~hx_revasc, ~year, survey_cardiac_risk, svyciprop) 
df_revasc <- data.frame(year = c(2008:2013), 
                        var = "Prior revascularization",
                        mean = dat_revasc %>% select(starts_with("hx")) %>% unlist,
                        SE = dat_revasc %>% select(starts_with("se."))%>% unlist
)

df_MI_feat <- rbind(df_inv, df_NSTEMI)

p_MI_feat <- ggplot(df_MI_feat, aes(x = year, y = mean)) + 
  geom_point(aes(color = var)) + 
  geom_line(aes(group = var, color = var), linetype = 2) +
  plot_themes + 
  xlab("") +
  labs(color = "") +
  ylab("Mean fraction")+  
  ggtitle("Type of AMI and management stratgey") 

p_died <- ggplot(df_died, aes(x = year, y = mean)) + 
  geom_point() + 
  geom_line(linetype = 2) +
  plot_themes + 
  xlab("") +
  labs(color = "") +
  ylim(0,0.2) + 
  ylab("Annual mortality rate")+  
  ggtitle("Mortality") 

if(save_plots){
  save_plot("Mortality.pdf", p_died, base_width = 6, base_height = 4)
  save_plot("Invasive_mgmt.pdf", p_MI_feat + ylim(0,1), base_width = 8, base_height = 6)
}
            