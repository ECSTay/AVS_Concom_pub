##############Model A analysis script - prob of reporting at least one AEFI/MA following either NIP or Men B alone AND concom

library(cmdstanr)
library(posterior)
library(ggplot2)

load("Z:/Analyses/Concomitant vaccination/Infant_NIP_MenB_concom_dat.rda")#load in the cleaned data for analysis

N_A                         ## number of responders
N_strat_A                   ## number of strategies
N_sched_A                   ## number of schedules

s_A <- dat$group             ## vaccine strategy, 1 = "Concomitant vaccination", 2 = "Separate"
t_A <- dat$schedule          ## schedule - 1 = 2 months, 2 = 4 months, 3 = 6 months, 4 = 12 months
w_A <- dat$sex               ## sex - 0 = "Male", 1 = "Female"
x_A <- dat$indig             ## Indigenous status -0 = Non-indig, 1 = Aboriginal and Torres Strait Islander
z_A <- dat$pmh               ## comorbidity - 0 = No, 1 = Yes
y_A <- dat$any_event         ## outcome - AEFI, MA or ?Fever - 0 = No, 1 = Yes


a <-  qlogis(0.3)          ## prior distribution mean - depends on the schedule and the vaccine strategy?
b <-                       ## prior distribution standard deviation


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




model_A <- cmdstan_model(write_stan_file(readLines("C:/Users/etay/Documents/Work documents/AVS work/Thuy_concom/model_A.stan")))###

fit_A <- model_A$sample(data = dat_A, 
                  chains = 8, parallel_chains = 8)###
postr <- posterior::as_draws_matrix(fit_A$draws())

##SAP code
model_A <- cmdstan_model("C:/Users/etay/Documents/Work documents/AVS work/Thuy_concom/model_A.stan")
drp <- capture.output({fit_A <- model_A$sample(dat_A, chains = 8, parallel_chains = 8)}) # fits the model and shows the output
draws_A <- as_draws_matrix(fit_A$draws(c("mu")))#posterior

dat_vis_A <- data.frame(x = plogis(as.vector(draws_A)),
                        Parameter = rep(c("Concomitant", "Separate"), each = 8000))

summary_fun_A <- function(strat, dat){
  mn <- format(round(mean(dat[dat$Parameter == strat,]$x), 2), nsmall = 2)
  cr <- format(round(quantile(dat[dat$Parameter == strat,]$x, c(0.025, 0.975)), 2), nsmall = 2)
  paste0(strat, " strategy: ", mn, " (", cr[1], ", ", cr[2], ")")
}

summary_fun_A("Concomitant", dat_vis_A)
summary_fun_A("Separate", dat_vis_A)
