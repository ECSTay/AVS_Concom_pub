##Thuy's concom/infant data - for Model A

library(tidyverse)
library(stringr)
library(data.table)
library(here)

readRenviron(here(".Renviron"))
RDS_PATH <- Sys.getenv("RDS_PATH")

read_rds(RDS_PATH, file = "/Analyses/Concomitant vaccination/Infant_NIP_MenB_concom_modelling.rda")

load("Z:/Analyses/Concomitant vaccination/Infant_NIP_MenB_concom_modelling.rda")
str(infant)

dat <- as.data.table(infant)
colnames(dat) <- tolower(colnames(dat))
dat[is.na(dat)] <- 0
dat <- dat[!(is.na(any_event)) &
                           !(is.na(atsi)) &
                           !(is.na(sex)) & 
                           !(is.na(group)) &
                           !(is.na(pmh)), ]




table(dat$group)
#Concomitant vaccination    Seperate vaccination 
#8749                    2352 
dat$group <- str_replace_all(dat$group, c("Concomitant vaccination" = "1", "Separate vaccination" = "2"))
dat$sex <- str_replace_all(dat$sex, c("Female" = "1", "Male" = "2"))
setnames(dat, "atsi", "indig")
dat[, indig := indig + 1]
dat$schedule <- str_replace_all(dat$schedule, c("2 mths" = "1", "4 mths" = "2",
                                                 "6 mths" = "3", "12 mths" = "4"))





#changing to factors
dat[,
           `:=`(
             sex = factor(sex, levels = c("1", "2")), #Female = 1
             indig = factor(indig, levels = c("1", "2")),#Non_aboriginal = 1
             group = factor(group, levels = c("1", "2")),#Concom = 1, Separate = 2, 
             schedule = as.factor(schedule)
             
           )]

