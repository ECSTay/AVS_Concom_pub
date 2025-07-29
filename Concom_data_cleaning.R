##Model A prob of reporting at least one MA following either(1) Concom and (2) NIP or Men B separately
#estimating P(at least one MA)
library(tidyverse)
library(stringr)
library(data.table)
library(here)

load("Z:/Analyses/Concomitant vaccination/Infant_NIP_MenB_concom_modelling_updated.rda")#connect to USyd RDS)
str(infant)

dat <- as.data.table(infant)
colnames(dat) <- tolower(colnames(dat))
#333
dat <- dat[!(is.na(any_event)) &
                           !(is.na(schedule)) &
                           !(is.na(atsi)) &
                           !(is.na(sex)) & 
                           !(is.na(group)) &
                           !(is.na(pmh)), ]
dat$vax_time_diff[is.na(dat$vax_time_diff)] <- 0
dat <- dat[dat$vax_time_diff != "367",]

table(dat$medical_attention)

# 0     1 
# 10852   236 
#use Model A for medical attention - no responses for 2 MA

table(dat$medical_attention, dat$schedule)

#     12 months 2 months 4 months 6 months
# 0      2993     2734     3747     1377
# 1        89       37       81       29

table(dat$group)

table(dat$medical_attention, dat$group)

#         Concomitant vaccination Seperate vaccination
# 0                    8560                 2291
# 1                     178                   58


dat$group <- str_replace_all(dat$group, c("Concomitant vaccination" = "1", "Seperate vaccination" = "2"))
dat$group <- as.integer(dat$group)
dat$schedule <- str_replace_all(dat$schedule, c("2 months" = "1", "4 months" = "2",
                                                "6 months" = "3","11" = "4"))
dat$schedule <- as.integer(dat$schedule)
dat$sex <- str_replace_all(dat$sex, c("Female" = "1", "Male" = "0"))
dat$sex <- as.numeric(dat$sex)
setnames(dat, "atsi", "indig")

#dat$any_event <- as.integer(dat$any_event)
#dat$impact <- as.integer(dat$impact)
#at least one any_event
#at <- dat %>%
# mutate(any_event = case_when(any_event == "2" ~ 1, .default = any_event))

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
##########################################

#MOdel C data sim from SAP
N_C <- 1000
N_strat_C <- 2
N_sched_C <- 1
s_C <- sample(1:N_strat_C, size = N_C, replace = TRUE)
t_C <- rep(1:N_sched_C, N_C/N_sched_C)
w_C <- sample(c(0,1), size = N_C, replace = TRUE)
x_C <- sample(c(0,1), size = N_C, replace = TRUE)
z_C <- sample(c(0,1), size = N_C, replace = TRUE)
mu_C <- array(NA, dim = c(N_strat_C, N_sched_C, 2))
mu_C[,,1] <- qlogis(0.3)
mu_C[2,,2] <- qlogis(0.05)
beta1_C <- beta2_C <- 0
gamma1_C <- gamma2_C <- 0
delta1_C <- delta2_C <- 0
p_C <- eta_C <- matrix(0, nrow = N_C, ncol = 3)
for(i in 1:N_C){
  eta_C[i, 2] <- mu_C[s_C[i], t_C[i], 1] +
    w_C[i]*beta1_C +
    x_C[i]*gamma1_C +
    z_C[i]*delta1_C
  if(s_C[i] == 1){
    eta_C[i, 3] <- -Inf
  } else {
    eta_C[i, 3] <- mu_C[s_C[i], t_C[i], 2] +
      w_C[i]*beta2_C +
      x_C[i]*gamma2_C +
      z_C[i]*delta2_C
  }
  p_C[i,] <- exp(eta_C[i,])/sum(exp(eta_C[i,]))
}
y_C <- sapply(1:N_C, function(i) sample(1:3, size = 1, prob = p_C[i,]))
a <- qlogis(0.3)
b <- d <- 2
c <- qlogis(0.05)
dat_C <- list(N = N_C,
              N_strat = N_strat_C,
              N_sched = N_sched_C,
              s = s_C, #strategy/vax seqeunce - int
              t = t_C, #schedule - int
              w = w_C, #sex
              x = x_C, #indig
              z = z_C, #comorbidity
              y = y_C, #outcome - int
              a = a,
              b = b,
              c = c,
              d = d)

str(dat_C)
##########################################################################################################
##Model C prob of reporting at least one AEFI/IMPACT  following (1)Concom, (2) NIP first or (3) Men B first

library(tidyverse)
library(stringr)
library(data.table)
library(here)

load("Z:/Analyses/Concomitant vaccination/Infant_NIP_MenB_concom_modelling_updated.rda")#connect to USyd RDS)
str(infant)

dat <- as.data.table(infant)
colnames(dat) <- tolower(colnames(dat))
#11101
dat <- dat[!(is.na(any_event)) &
             !(is.na(schedule)) &
             !(is.na(atsi)) &
             !(is.na(sex)) & 
             !(is.na(group)) &
             #!is.na(fever)) & 4 records missing a value for fever
             !(is.na(pmh)), ]#11088
dat$vax_time_diff[is.na(dat$vax_time_diff)] <- 0

dat <- dat[dat$vax_time_diff != "367",] #11087

dat$group <- str_replace_all(dat$group, c("Concomitant vaccination" = "1", "Seperate vaccination" = "2"))
dat$group <- as.integer(dat$group)

dat$vax_sequence[is.na(dat$vax_sequence)] <- "Concomitant vaccination"
dat$vax_sequence <- str_replace_all(dat$vax_sequence, c("Concomitant vaccination" = "1", "NIP first" = "2", "MenB first" = "3"))
dat$vax_sequence <- as.integer(dat$vax_sequence)
table(dat$vax_sequence)

#1    2    3 
#8738 2208  141 

dat$any_event <- as.integer(dat$any_event)
table(dat$any_event)
#   0    1    2 
# 6375 4339  371 

dat$impact <- as.integer(dat$impact)
table(dat$impact)
#   0     1     2 
# 10632   449     7 

dat$schedule <- str_replace_all(dat$schedule, c("2 months" = "1", "4 months" = "2",
                                                "6 months" = "3","11" = "4"))
dat$schedule <- as.integer(dat$schedule)
table(dat$vax_sequence, dat$schedule)

#     1    2    3    4
# 1 2118 3113 1156 2351
# 2  653  630  242  683
# 3    0   85    8   48

dat$sex <- str_replace_all(dat$sex, c("Female" = "1", "Male" = "0"))
dat$sex <- as.numeric(dat$sex)
setnames(dat, "atsi", "indig")

write.csv(dat, file = "C:/Users/ETay/Documents/Work documents/AVS work/Thuy_concom/dat_modelC.csv", row.names = FALSE)
