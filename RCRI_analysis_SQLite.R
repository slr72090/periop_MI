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
# Identify procedures that correspond to non-cardiac surgeries
year_vec = c(2008:2013)
n_procedures_vec <- rep(15,length(year_vec)) #rep(15,length(year_vec)) # Number of possible procedures per individual for this year 
n_dx_vec <- c(15, 25, 25, 15, 25, 25) # Number of possible diagnoses per individual for this year 
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
  select(-c(prior_CABG, prior_PCI, prior_MI, CAD, transplant,thoracic_surgery,vascular,nchronic)) %>% 
  mutate(`RCRI >= 3` = as.factor(as.numeric(RCRI_pt) > 3))
  #select(-c(prior_CABG, prior_PCI, transplant, thoracic_surgery, vascular))

## Exploratory plots ## ------------------------------------------------------------------------------------------------
textSize = 12
save_plots = F
source("plot_themes.R")

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

#### Evaluate correlations in covariates ## ---------------------------------------------------------------------------------------------
cor_data <- data %>% select(-c(ind, year, los, invasive_mgmt, died, Ischemic_stroke, RCRI_pt, age, `RCRI >= 3`, cm_mets)) %>% 
  apply(.,2,as.numeric)
M <- cor(cor_data)
diag(M) <- 0
p_cor <- corrplot(M,
                  method = "shade",
                  type = "lower", 
                  order = "hclust",
                  tl.col = "black",
                  tl.cex = .75)

## Step 1: retain all variables ## -----------------------------------------------------------------------------------------
dat_PCA <- data %>% mutate(RCRI_pt = as.numeric(scale(as.numeric(RCRI_pt))))
split.data.all <- splitmix(dat_PCA[,!(colnames(dat_PCA) %in% c("ind",
                                                        "year",
                                                        #"ndx",
                                                        "RCRI_pt", 
                                                        "RCRI >= 3",
                                                        "age_factor",
                                                        "died",
                                                        "los", 
                                                        "MI", 
                                                        "prior_PCI",
                                                        "prior_CABG",
                                                        "cm_mets",
                                                        #"hx_isch_heart",
                                                        #"hx_revasc",
                                                        "Ischemic_stroke", 
                                                        "thoracic_surgery", 
                                                        "invasive_mgmt",
                                                        "transplant",
                                                        "vascular"))]
)


split.data.RCRI <- splitmix(data[,(colnames(data) %in% c("age", "race", "gender", "hx_CVA", "hx_DM", "hx_isch_heart", "hx_ckd", "hx_chf", "high_risk_surgery"))])

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

out.coord$sig <- as.factor(out.coord$`dim 1` > 0.5 | out.coord$`dim 2` > 0.5 |out.coord$`dim 1` < -0.5 | out.coord$`dim 2` < -0.5 )

## Plotting -----------------------------------------------------------------------------------
textSize = 12
source("plot_themes.R")

## Scree plot

p_scree <- ggplot(out.eig, aes(x = dim, y = Eigenvalue)) + 
  geom_point() + 
  geom_line(linetype = 2) +
  geom_hline(yintercept = 1, color = "red", linetype = 2) + 
  xlab("Principal component") +
  plot_themes 

# Loadings plots ## --------------------------------
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

#p_loadings_dim1 <- ggplot()+
  #geom_point(data = out.loadings, aes(x=`dim 1`, y= `dim 2`, colour = "#DCDCDC"))+
#  geom_bar(data = out.coord, aes(x=var,y=`dim 1`), stat = "identity")+ 
#  plot_themes + 
 # theme(axis.text.x = element_text(angle = 90)) + xlab("")

p_squared_loadings = ggplot()+
  #geom_point(data = out.loadings, aes(x=`dim 1`, y= `dim 2`, colour = "#DCDCDC"))+
  geom_segment(data = out.loadings, aes(x=0,y=0,xend=`dim 1`,yend=`dim 2`),
               arrow=arrow(length=unit(0.1,"cm")), color = "black") +
  geom_text(data = out.loadings, aes(x=`dim 1`, y=`dim 2`, label=factor(rownames(out.loadings))),color="red") + 
  xlab(xlabel) + 
  ylab(ylabel) + 
  plot_themes

if(save_plots){
  save_plot("scree_plot.pdf", p_scree, base_width = 6, base_height = 4)
  save_plot(paste0("loadings.pdf"), p_loadings, base_width = 8, base_height = 6)
}

## Individuals plot ## --------------------------------------
p_ind_RCRI <-  qplot(data = out.ind, x = `dim 1`, y = `dim 2`, colour = `RCRI >= 3`) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = `RCRI >= 3`)) +
  plot_themes

p_ind_year <-  qplot(data = out.ind, x = `dim 1`, y = `dim 2`, colour = year) +
  stat_ellipse(geom = "polygon", alpha = .2, aes(fill = year)) +
  plot_themes

p_ind_died <- ggplot(out.ind, aes(x = `dim 1`, y = `dim 2`)) + 
  geom_point(aes(color = factor(died))) + 
  plot_themes

if(save_plots){
  ind_plots <- plot_grid(p_ind_RCRI, p_ind_year, ncol =2)
  save_plot("ind_plots.pdf", ind_plots, base_width = 12, base_height = 4)
}



## Feature selection -------------------------------------------------------------------------------------------

#Select retained variables
out.contrib$max_contrib <- apply(out.contrib %>% select(-var),1,max) 
retained_vars <- out.contrib %>% filter(max_contrib > 10) %>% select(var)

## Logit Regression analysis ## ------------------------------------------------------------------------------

sub.PCA <- data[,(colnames(data) %in% c("died", "year", "race", "gender", "invasive_mgmt", unlist(retained_vars)))]
out.PCA <- glm(as.numeric(died) ~., data = sub.PCA, family = "binomial")
#summary(out.PCA)
#exp(cbind(OR = coef(out.PCA), confint(out.PCA)))

#dat.PCA.vars <- out.ind %>% select(c(contains("dim"), "year", "died"))
#out.PCA.comps <- glm(as.numeric(died) ~., data = dat.PCA.vars, family = "binomial")

var_names_full <- c(names(X1), names(X2))
data_full <- data %>% select(c(year, died, invasive_mgmt, var_names_full))# %>% select(-c(RCRI_pt, Ischemic_stroke, los, ind, age_factor, `RCRI > 3`))
out.full <- glm(died ~., data = data_full, family = "binomial")

# out.1 <- glm(as.numeric(died) ~ age + 
#                as.factor(gender) + 
#                as.factor(race) + 
#                as.factor(obesity) +
#                smoking + 
#                as.factor(HTN) + 
#                as.factor(HLD) + 
#                as.factor(hx_DM) + 
#                as.factor(hx_ckd) + 
#                as.factor(hx_isch_heart) + 
#                as.factor(hx_revasc) + 
#                as.factor(PAD) + 
#                as.factor(hx_chf) + 
#                as.factor(valve_dz) + 
#                as.factor(chrnlung) + 
#                as.factor(malignancy) + 
#                as.factor(anemia) + 
#                as.factor(alcoholic) + 
#                as.factor(high_risk_surgery) + 
#                as.factor(invasive_mgmt) + 
#                as.factor(year), data = data, family = "binomial")

t1 <- Sys.time()
out.step <- out.full %>% 
  stepwise(out.full, 
       direction = "backward",
       criterion = "BIC",
       trace = F)
  #stepAIC(trace = F)
t2 <- Sys.time()
print(t2 - t1)

#exp(cbind(OR = coef(out.step), confint(out.step)))

out.RCRI <- glm(as.numeric(died) ~  age +
               gender + 
               race + 
               year + 
               smoking + 
               HLD+ 
               hx_isch_heart +
               hx_chf +
               hx_CVA + 
               hx_DM + 
               hx_ckd + 
               high_risk_surgery + 
               invasive_mgmt, data = data, family = "binomial")

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

model_comp <- data.frame(model = c("Full", "PCA", "RCRI"), 
                         BIC = c(BIC(out.full), BIC(out.PCA), BIC(out.RCRI))) %>% 
  mutate(delta = BIC - min(BIC))


## Testing model out of sample ## -------------------------------------------

set.seed(364)
sample <- sample(nrow(data),floor(nrow(data)*0.8))
train <- data[sample,]
test <- data[-sample,]


library(pROC)

out.RCRI.train <- glm(as.numeric(died) ~ age + 
                        gender + 
                        race + 
                        year + 
                        smoking +  
                        HLD + 
                        hx_ckd + 
                        hx_isch_heart + 
                        hx_chf +
                        high_risk_surgery + 
                        invasive_mgmt, data = train, family = "binomial")

out.full.train <- glm(died ~., data = (train %>% select(c(year, died, invasive_mgmt, var_names_full))), family = "binomial")

test_prob_RCRI = predict(out.RCRI, newdata = test, type = "response")
test_roc_RCRI = roc(test$died ~ test_prob_RCRI, plot = TRUE, print.auc = TRUE)

test_prob_full = predict(out.full, newdata = test, type = "response")
test_roc_full = roc(test$died ~ test_prob_full, plot = TRUE, print.auc = TRUE)

p_roc <- ggroc(list(`Full model` = test_roc_full, `RCRI model` = test_roc_RCRI), linetype = 2) +
  plot_themes + 
  labs(color = "Model") + 
  ggtitle("ROC curves")

if(save_plots){
  save_plot("ROC_curves.pdf", p_roc, base_width = 8, base_height = 4)
}
               