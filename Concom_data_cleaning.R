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

write.csv(dat, file = "C:/Users/ETay/Documents/Work documents/AVS work/Thuy_concom/dat_modelA.csv", row.names = FALSE)


##########################
#sim code from the SAP
N_A <- 1000
N_strat_A <- 2
N_sched_A <- 1
s_A <- sample(1:N_strat_A, size = N_A, replace = TRUE)
t_A <- rep(1:N_sched_A, N_A/N_sched_A)
w_A <- sample(c(0,1), size = N_A, replace = TRUE)
x_A <- sample(c(0,1), size = N_A, replace = TRUE)
z_A <- sample(c(0,1), size = N_A, replace = TRUE)
mu_A <- matrix(qlogis(0.3), nrow = N_strat_A, ncol = N_sched_A)
beta_A <- gamma_A <- delta_A <- 0
p_A <- sapply(1:N_A, function(i) plogis(mu_A[s_A[i], t_A[i]] +
                                          w_A[i]*beta_A +
                                          x_A[i]*gamma_A +
                                          z_A[i]*delta_A))
y_A <- rbinom(N_A, 1, p_A)
a <- qlogis(0.3)
b <- 2

dat_A <- list(N = N_A,
              N_strat = N_strat_A,
              N_sched = N_sched_A,
              s = s_A,
              t = t_A,
              w = w_A,
              x = x_A,
              z = z_A,
              y = y_A,
              a = a,
              b = b)
str(dat_A)





