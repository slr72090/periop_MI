## Exploratory data analysis; Perioperative MI 
## Sylvia Ranjeva 

## Load package scripts -----------------------------------------------
require(ggplot2)
require(foreign)
require(readstata13)
require(dplyr)
require(FactoMineR)
source("dta_to_sqlite.R")

## Specify data ----------------------------
this_year = 2010
this_db = "NIS.db"
this_file = "./data_files/Core_Hospital_2010.dta"


dta_to_sqlite(dbfilename = this_db, data_year = this_year, filename = this_file)