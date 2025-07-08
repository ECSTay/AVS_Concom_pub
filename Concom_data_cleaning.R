##Thuy's concom/infant data - for Model A

library(tidyverse)
library(stringr)
library(data.table)
library(here)

load("Z:/Analyses/Concomitant vaccination/Infant_NIP_MenB_concom_modelling.rda")#connect to USyd RDS)
str(infant)

dat <- as.data.table(infant)
colnames(dat) <- tolower(colnames(dat))

table(dat$group)
#Concomitant vaccination    Seperate vaccination 
#8749                    2352 
dat$group <- str_replace_all(dat$group, c("Concomitant vaccination" = "1", "Seperate vaccination" = "2"))
dat$group <- as.integer(dat$group)
dat$schedule <- str_replace_all(dat$schedule, c("2 months" = "1", "4 months" = "2",
                                                "6 months" = "3", "12 months" = "4"))
dat$schedule <- as.integer(dat$schedule)
dat$sex <- str_replace_all(dat$sex, c("Female" = "1", "Male" = "0"))
dat$sex <- as.numeric(dat$sex)
setnames(dat, "atsi", "indig")
dat$any_event <- as.integer(dat$any_event)






