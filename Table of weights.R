##table of weights
library(data.table)
library(ggplot2)
library(tidyverse)
library(stringr)
library(posterior)

#load in the relevant posterior
draws <- readRDS(file ="C:/Users/ETay/Documents/postr_concom_AEFI.rds")

dat_C <- readRDS("C:/Users/ETay/Documents/dat_C_AEFI.rds")

#load in dat
dat <- fread(file = "C:/Users/ETay/Documents/NIP_MenB_dat.csv")

sex_weights <- dat[, .(mn = mean(sex)), by = schedule]
indig_weights <- dat[, .(mn = mean(indig)), by = schedule]
state_weights <- data.table(mn = as.vector(prop.table(table(dat$schedule,dat$clinic_state), margin = 1)[,-2]),
                            schedule = rep(1:4, 7),
                            state = rep(1:7, each = 4))
clinic_weights <- data.table(mn = as.vector(prop.table(table(dat$schedule,dat$clinic_type), margin = 1)[,-2]),
                             schedule = rep(1:4, 2),
                             clinic = rep(1:2, each = 4))
comorb_weights <- dat[, .(mn = mean(pmh)), by = schedule]

sex <- as.data.frame(t(sex_weights[,-1]))
indig <- as.data.frame(t(indig_weights[,-1]))

state <- reshape(state_weights, idvar = "schedule", timevar = "state", direction = "wide")
clinic <- reshape(clinic__weights, idvar = "schedule", timevar = "state", direction = "wide")
co_morb <- as.data.frame(t(indig_weights[,-1]))