## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

# Set the working directory to the source file directory
setwd("~/Desktop/perip_MI_GITHUB/mortality_dispo_models")

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
source("../utility_scripts/plot_themes.R")
dbfilename = "../data_files/NIS.db"

# Identify procedures that correspond to non-cardiac surgeries
year_vec = c(2008:2013)
n_procedures_vec <- rep(15,length(year_vec)) #rep(15,length(year_vec)) # Number of possible procedures per individual for this year 
n_dx_vec <- c(15, 25, 25, 25, 25, 25)

pr_mapping <- read.csv(".././data_files/pc2015.csv") %>% 
  as.data.frame()
names(pr_mapping) <- c("ICD9", "cat")

## Generate string of ICD9 procedure codes that correspond to "major therapeutic interventions"
ICD9_procedure_cats <- pr_mapping %>% 
  filter(cat == 4) %>% 
  select(ICD9) %>% 
  unlist()

ICD9_string <- as.character(ICD9_procedure_cats)
pr_codes <- paste("(", toString(paste(ICD9_string,sep="")), ")", sep="")

## Query 1: n discharges 2008 - 2013 ## ---------------------------------------------------------------------

count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  nis_query <- "SELECT COUNT(*) FROM data_y "
  nis_query <- gsub("y", toString(year), nis_query)
# Perform the query
t <- Sys.time()
db <- dbConnect(RSQLite::SQLite(), dbfilename) 
q <- RSQLite::dbSendQuery(db, nis_query)
# Get the data
count_df_sub <- RSQLite::dbFetch(q, n=-1)
t2 <- Sys.time()
# Clear the query and disconnect from database to prevent memory leaks
RSQLite::dbClearResult(q)
dbDisconnect(db)
time_elapsed = t2 - t
print(time_elapsed)
count_df <- rbind(count_df, data.frame(yr = year, 
                                       count  = as.numeric(count_df_sub$`COUNT(*)`)))
}

## Query 2: n discharges 2008 - 2013, age > 45 and major surgery ## ---------------------------------------------------------------------
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND age > 45"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  t <- Sys.time()
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}

## Query 3: n discharges 2008 - 2013, age > 45 and major surgery, elective admissions ## ---------------------------------------------------------------------
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND age > 45 AND elective == 1"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  t <- Sys.time()
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}


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
                     174, # non-OR skin/breast,
                     229 #non-operative removal of foreign body
)

pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")
pr_ccs_snippet <- " prccs1 NOT IN pr_ccs_excluded"


## Query 4: n discharges 2008 - 2013, age > 45 and major surgery, elective admissions: cardiac procedures ## ---------------------------------------------------------------------
# Exlcude primary procedure codes:
pr_ccs_excluded <- c(43,44,45,46,48,49,50,52)
pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  pr_ccs_snippet <- " prccs1 IN pr_ccs_excluded"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND (pr_ccs_snippet) AND age > 45 AND elective == 1"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  nis_query <- gsub("pr_ccs_snippet", pr_ccs_snippet, nis_query)
  nis_query <- gsub("pr_ccs_excluded", pr_ccs_excluded_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  t <- Sys.time()
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}

## Query 5: n discharges 2008 - 2013, age > 45 and major surgery, elective admissions: bone marrow procedures ## ---------------------------------------------------------------------
# Exlcude primary procedure codes:
pr_ccs_excluded <- c(64,65)
pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  pr_ccs_snippet <- " prccs1 IN pr_ccs_excluded"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND (pr_ccs_snippet) AND age > 45 AND elective == 1"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  nis_query <- gsub("pr_ccs_snippet", pr_ccs_snippet, nis_query)
  nis_query <- gsub("pr_ccs_excluded", pr_ccs_excluded_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  t <- Sys.time()
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}

## Query 6: n discharges 2008 - 2013, age > 45 and major surgery, elective admissions: optho procedures ## ---------------------------------------------------------------------
# Exlcude primary procedure codes:
pr_ccs_excluded <- c(13:21,29)
pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  pr_ccs_snippet <- " prccs1 IN pr_ccs_excluded"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND (pr_ccs_snippet) AND age > 45 AND elective == 1"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  nis_query <- gsub("pr_ccs_snippet", pr_ccs_snippet, nis_query)
  nis_query <- gsub("pr_ccs_excluded", pr_ccs_excluded_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  t <- Sys.time()
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}

## Query 7: n discharges 2008 - 2013, age > 45 and major surgery, elective admissions: dental procedures ## ---------------------------------------------------------------------
# Exlcude primary procedure codes:
pr_ccs_excluded <- c(211 )
pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  pr_ccs_snippet <- " prccs1 IN pr_ccs_excluded"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND (pr_ccs_snippet) AND age > 45 AND elective == 1"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  nis_query <- gsub("pr_ccs_snippet", pr_ccs_snippet, nis_query)
  nis_query <- gsub("pr_ccs_excluded", pr_ccs_excluded_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  t <- Sys.time()
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}



## Query 8: n discharges 2008 - 2013, age > 45 and major surgery, elective admissions: radiation therapy ## ---------------------------------------------------------------------
# Exlcude primary procedure codes:
pr_ccs_excluded <- c(211 )
pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  pr_ccs_snippet <- " prccs1 IN pr_ccs_excluded"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND (pr_ccs_snippet) AND age > 45 AND elective == 1"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  nis_query <- gsub("pr_ccs_snippet", pr_ccs_snippet, nis_query)
  nis_query <- gsub("pr_ccs_excluded", pr_ccs_excluded_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  t <- Sys.time()
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}



## Query 9: n discharges 2008 - 2013, age > 45 and major surgery, elective admissions: non-OR procedures ## ---------------------------------------------------------------------
# Exlcude primary procedure codes:
pr_ccs_excluded <- c(8,
                     41, 
                     62,63, # diagnostic and non-therapeutic cardiac procedures
                     95, # non-OR GI
                     111, 117, 131, # non-OR GU
                     163, # non-OR MSK
                     174, # non-OR skin/breast,
                     229 #non-operative removal of foreign body
)

pr_ccs_excluded_codes <- paste("(", toString(paste(pr_ccs_excluded,sep="")), ")", sep="")
count_df <- data.frame()
for(i in c(1:length(year_vec))){
  year = year_vec[i]
  print(year)
  pr_snippet <- " pr1 IN pr_codes"
  pr_ccs_snippet <- " prccs1 IN pr_ccs_excluded"
  nis_query <- "SELECT COUNT(*) FROM data_y WHERE (pr_snippet) AND (pr_ccs_snippet) AND age > 45 AND elective == 1"
  nis_query <- gsub("y", toString(year), nis_query)
  nis_query <- gsub("pr_snippet", pr_snippet, nis_query)
  nis_query <- gsub("pr_codes", pr_codes, nis_query)
  nis_query <- gsub("pr_ccs_snippet", pr_ccs_snippet, nis_query)
  nis_query <- gsub("pr_ccs_excluded", pr_ccs_excluded_codes, nis_query)
  
  # Perform the query
  db <- dbConnect(RSQLite::SQLite(), dbfilename) 
  t <- Sys.time()
  q <- RSQLite::dbSendQuery(db, nis_query)
  # Get the data
  count_df_sub <- RSQLite::dbFetch(q, n=-1)
  t2 <- Sys.time()
  # Clear the query and disconnect from database to prevent memory leaks
  RSQLite::dbClearResult(q)
  dbDisconnect(db)
  time_elapsed = t2 - t
  print(time_elapsed)
  count_df <- rbind(count_df, data.frame(yr = year, 
                                         count  = as.numeric(count_df_sub$`COUNT(*)`)))
}


