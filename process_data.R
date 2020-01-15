## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

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

## Functions to identify relevant diagnoses and procedures by ICD9 codes

# Identify diagnoses
dx_fun <- function(string, x){
  return(as.integer(sum(grepl(string, x)) > 0))
}

# Identify procedures
pr_fun <- function(code, x){
  val <- sum(as.integer(x == code), na.rm=T)
  if(val == 0){
    return(val)
  }
  if(val > 0){
    return(1)
  }
}

# Query the database to obtain the relevant dataframe ## ------------------------------------------------------------------------

pr_mapping <- read.csv("./data_files/pc2015.csv") %>% 
  as.data.frame()
names(pr_mapping) <- c("ICD9", "cat")

## Generate string of ICD9 procedure codes that correspond to "major therapeutic interventions"
ICD9_procedure_cats <- pr_mapping %>% 
  filter(cat == 4) %>% 
  select(ICD9) %>% 
  unlist()

ICD9_string <- as.character(ICD9_procedure_cats)
pr_codes <- paste("(", toString(paste(ICD9_string,sep="")), ")", sep="")

# Exlcude primary procedure codes:
pr_ccs_excluded <- c(8, # non-OR nervous system
                     13:21, # optho procedures,
                     29, # dental procedures,
                     41, # non-OR respiratory,
                     43,44,45,46,48,49,50,52, # cardiac
                     62,63, # diagnostic and non-therapeutic cardiac procedures
                     64,65, # bone marrow
                     95, # non-OR GI
                     111, 117, 131, # non-OR GU
                     163, # non-OR MSK
                     174 # non-OR skin/breast,
                     #229 #non-operative removal of foreign body
)

pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")

data_all <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  n_procedures = n_procedures_vec[i]
  n_dx = n_dx_vec[i]
  print(year)
  
  pr_snippet <- " pr1 IN pr_codes"
  for(j in c(2:n_procedures)){
   new_txt <- paste0(" OR pr",j," IN pr_codes")
   pr_snippet <- paste0(pr_snippet, new_txt)
  }
  
  pr_ccs_snippet <- " prccs1 NOT IN pr_ccs_excluded"
  
  dx_snippet <- " dx1 LIKE '410%'"
  for(j in c(2:n_dx)){
    new_txt <- paste0(" OR dx",j," LIKE '410%'")
    dx_snippet <- paste0(dx_snippet, new_txt)
  }
  nis_query <- "SELECT * FROM data_y WHERE (pr_snippet) AND (pr_ccs_snippet) AND (dx_snippet) AND age > 45 AND elective == 1 LIMIT 1000"

  #Sub correct year and procedure codes 
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  nis_query <- gsub("pr_ccs_snippet", pr_ccs_snippet, nis_query)
  nis_query <- gsub("pr_ccs_excluded", pr_ccs_excluded_codes, nis_query)
  nis_query <- gsub("dx_snippet", dx_snippet, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), "NIS.db") 
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  t <- Sys.time()
  core_df <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  cat("nrow is ", nrow(core_df), "\n")
  
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  
  # test_dx <- core_df %>% select(starts_with("dx")) %>% 
  #select(-starts_with("dxccs"))
  
  ## ASCVD risk factors ## --------------------------------------------------------------------------
  
  # Hyperlipidemia 
  core_df$HLD <- apply(core_df %>% select(starts_with("dx")), 
                       1, 
                       dx_fun,  
                       string = "2720|2721|2722|2723|2724")
  
  
  # Tobacco 
  core_df$smoking <- apply(core_df %>% select(starts_with("dx")), 
                           1, 
                           dx_fun,
                           string = "V1582|3051")
  
  # CKD
  core_df$ckd <- apply(core_df %>% select(starts_with("dx")), 
                       1, 
                       dx_fun,
                       string = "^585|7925|^V420|^V451|^V56|9580")
  
  # DM 
  core_df$hx_DM <- apply(core_df %>% select(starts_with("dx")), 
                       1, 
                       dx_fun,
                       string =  "24900|25000|25001|V5391|V6546|24901|24910|24911|24920|24921|24930|24931|24940|24941|24950|24951|24960|24961|24970|24971|24980|24981|24990|24991|^250"
  )
  # Malignancy
 # malig_string <- c("588")
 # for(n in c(140:172,174:195,200:208)){
  # malig_string <- c(malig_string, paste0("|",n))
 # }
  
  #malig_string <- c(588,140:172,174:195,200:208)
  malig_string <- c(11:42)
  
  core_df$malignancy <- apply(core_df %>% 
                                select(starts_with("dxccs")),
                              1, 
                              pr_fun,
                              code = malig_string)
                          
  
  ## Prior VTE
  core_df$hx_VTE <- apply(core_df %>% select(starts_with("dx")), 
                              1, 
                              dx_fun,
                              string = "V1251")
  
  ## ASCVD ## --------------------------------------------------------------------------
  
  core_df$CAD <- apply(core_df %>% select(starts_with("dx")), 
                       1, 
                       dx_fun,
                       string = "^414|^V4581|^V4582")
  
  core_df$PAD <- apply(core_df %>% select(starts_with("dx")), 
                       1, 
                       dx_fun,
                       string = "^440|4439|^557")
  
  core_df$prior_MI <- apply(core_df %>% select(starts_with("dx")), 
                            1, 
                            dx_fun,
                            string = "^412")
  
  core_df$prior_TIA <- apply(core_df %>% select(starts_with("dx")), 
                             1, 
                             dx_fun,
                             string = "V1254")
  
  core_df$hx_CVA <- apply(core_df %>% select(starts_with("dx")), 
                       1, 
                       dx_fun,
                       string = "^438") #V1254|
  
  core_df$prior_PCI <- apply(core_df %>% select(starts_with("dx")), 
                             1, 
                             dx_fun,
                             string = "V4582")
  
  core_df$prior_CABG <- apply(core_df %>% select(starts_with("dx")), 
                              1, 
                              dx_fun,
                              string = "V4581")
  
  core_df$CHF_sys <- apply(core_df %>% select(starts_with("dx")), 
                          1, 
                          dx_fun,
                          string = "^4282")
  
  core_df$CHF_dia <- apply(core_df %>% select(starts_with("dx")), 
                           1, 
                           dx_fun,
                           string = "^4283")
  
  core_df$CHF_com <- apply(core_df %>% select(starts_with("dx")), 
                           1, 
                           dx_fun,
                           string = "^4284")
  
  core_df$CHF_unsp <- apply(core_df %>% select(starts_with("dx")), 
                           1, 
                           dx_fun,
                           string = "^4280")
  
  
  ##  Surgery_type or intervention ## ---------------------------------------------------------------------
  core_df$thoracic_surgery <- sapply(core_df$prccs1,# %>%  select(starts_with("prccs")), 
                                     #1, 
                                     pr_fun,
                                     code = c(36:40,42))
  
  core_df$transplant <- apply(core_df %>% select(starts_with("prccs")), 
                              1, 
                              pr_fun,
                              code = c(105,176))
  
  core_df$vascular <- apply(core_df %>% select(starts_with("prccs")), 
                            1, 
                            pr_fun,
                            code = c(54:56))
  
  core_df$invasive_mgmt <- apply(core_df %>% select(starts_with("prccs")), 
                                 1, 
                                 pr_fun,
                                 code = c(44,45,46,47)) 
  
  ## Use of balloon pump
  core_df$IABP <- apply(core_df %>% select(starts_with("pr")) %>% 
                                             select(-starts_with("prccs")), 
                                 1, 
                                 pr_fun,
                                 code = 3761) 

## MINS/VISION post-op covariates
  core_df$Afib <- apply(core_df %>% select(starts_with("dx")), 
                        1, 
                        dx_fun,  
                        string = "^4273")
  
  core_df$sepsis <- apply(core_df %>% select(starts_with("dx")), 
                          1, 
                          dx_fun,  
                          string = "^9959|78552")
  
  core_df$PNA <-  apply(core_df %>% select(starts_with("dx")), 
                        1, 
                        dx_fun,  
                        string = "^480|^481|^482|^483|^484|^485|^486|^4870")
  
  core_df$PE <- apply(core_df %>% select(starts_with("dx")), 
                      1, 
                      dx_fun,  
                      string = "^4151")
  
  core_df$Bleed <- apply(core_df %>% select(starts_with("dx")), 
                      1, 
                      dx_fun,  
                      string = "^2851|^99811") # Acute post-hemorrhage anemia, hemorrhage complicating procedure 
  
  core_df$DVT <- apply(core_df %>% select(starts_with("dx")), 
                         1, 
                         dx_fun,  
                         string = "^4534|^4538|4539") # Acute DVT
  
  ##  Outcomes ## ---------------------------------------------------------------------
  
  ## NSTEMI vs STEMI 
  core_df$NSTEMI <- apply(core_df %>% select(starts_with("dx")), 
                                   1, 
                                   dx_fun,  
                                   string = "^4107")
  
  core_df$cardiogenic_shock <- apply(core_df %>% select(starts_with("dx")), 
                                     1, 
                                     dx_fun,  
                                     string = "78551")
  
  ## Other outcomes
  
  core_df$Ischemic_stroke <- apply(core_df %>% select(starts_with("dx")), 
                                   1, 
                                   dx_fun,  
                                   string = "^433|^434|^436|^4371")
  
  ## Dispo outcomes 
  
  core_df$ICF <- as.integer(core_df$dispuniform == 5)
  
  core_df <- core_df %>% mutate(
    # hx_DM = as.integer( cm_dmcx == 1), #cm_dm == 1 |
    hx_isch_heart = as.integer(prior_MI==1 | CAD == 1),
    hx_ckd = as.integer(ckd ==1 | cm_renlfail == 1),
    hx_chf = as.integer(CHF_dia == 1 | CHF_sys == 1 | CHF_unsp == 1 | CHF_com == 1)) %>% 
    mutate(RCRI_pt = select(., c(hx_CVA, hx_DM, hx_isch_heart, hx_ckd, hx_chf)) %>% apply(1, sum, na.rm=TRUE)) %>% 
    select(-c(prior_TIA,
              cm_dm, 
              cm_renlfail,
              cm_dmcx, 
              ckd, 
              cm_chf))
  
  data <-core_df %>% 
    mutate(age = as.numeric(age),
      nchronic = as.numeric(nchronic), # of chronic diagnoses
      died = as.numeric(died), # died during hospital admission
      los = as.numeric(los), # length of stay 
      race = as.integer(race == 2), # binary race variable: non-hispanic african american 
      gender = as.integer(female == 0),
      malignancy = as.integer(malignancy)
      ) %>% # male gender 
    mutate_if(is.integer, as.factor) %>%  
    select(c(age, 
             race, 
             gender,
             nchronic,
             cm_obese,
             smoking,
             cm_alcohol, 
             cm_htn_c, 
             HLD,
             hx_DM, 
             #cm_dm,
             #cm_dmcx,
             hx_ckd, 
             #cm_liver,
             #cm_renlfail,
             CAD,
             prior_MI,
             hx_isch_heart,
             prior_PCI,
             prior_CABG,
             cm_perivasc,
             cm_valve,
             hx_chf,
             # CHF_unsp,
             # CHF_sys,
             # CHF_dia,
             # CHF_com,
             hx_VTE,
             cm_chrnlung, 
             cm_alcohol,
             malignancy,
             cm_anemdef,
             malignancy,
             cm_mets,
             hx_CVA, 
             thoracic_surgery,
             transplant,
             vascular, 
             RCRI_pt,
             #Ischemic_stroke,
             Afib,
             sepsis,
             PNA,
             PE,
             DVT,
             Bleed,
             NSTEMI,
             died,
             invasive_mgmt,
             IABP,
             cardiogenic_shock,
             ICF,
             year))  
  
  if(drop_missing){
    data <- data %>% 
    drop_na()
  }## For now, remove missing data for analyses
  data_all <- rbind(data_all, data)
}

data_all$ind <- c(1:nrow(data_all))