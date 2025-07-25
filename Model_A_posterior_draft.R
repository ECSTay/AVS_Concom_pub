##############Model A analysis script - prob of reporting at least one MA following either Men B first, NIP first OR concom

library(cmdstanr)
library(posterior)
library(ggplot2)



N_A  = nrow(dat)                     ## number of responders
N_strat_A  = 2                      ## number of strategies #3
N_sched_A  = 4                       ## number of schedules

#s_A <- dat$vax_sequence              ## vaccine strategy - "Concomitant vaccination" = "1", "NIP first" = "2", "MenB first" = "3"
s_A <- dat$group                     ## vaccine strategy - "Concomitant vaccination" = "1", "Separate" = "2"
t_A <- dat$schedule                  ## schedule - 1 = 2 months, 2 = 4 months, 3 = 6 months, 4 = 12 months
w_A <- dat$sex                       ## sex - 0 = "Male", 1 = "Female"
x_A <- dat$indig                     ## Indigenous status -0 = Non-indig, 1 = Aboriginal and Torres Strait Islander
z_A <- dat$pmh                       ## comorbidity - 0 = No, 1 = Yes
y_A <- dat$medical_attention         ## outcome - MA - 0 = No, 1 = Yes


a <-  qlogis(0.3)           ## prior distribution mean - depends on the schedule and the vaccine strategy?
b <-   2                    ## prior distribution standard deviation


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




model_A <- cmdstan_model(write_stan_file(readLines("C:/Users/ETay/Documents/Work documents/AVS work/Thuy_concom/AVS_Concom_pub/model_A.stan")))###

fit_A <- model_A$sample(data = dat_A, 
                  chains = 4, parallel_chains = 4)###
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
