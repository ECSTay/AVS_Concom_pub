##############Model A analysis script - prob of reporting at least one AEFI/MA following either NIP or Men B alone AND concom

library(cmdstanr)
library(posterior)
library(ggplot2)

conds <- c("pmh_heart_dis", "pmh_bp", "pmh_diabetes", "pmh_lung",
           "pmh_obesity", "pmh_kidney", "pmh_liver", "pmh_cancer", "pmh_blood_cancer", 
           "pmh_chemo_rad", "pmh_transplant_organ", "pmh_transplant_bone", "pmh_neuro", "pmh_inflam", "pmh_immunodef")
cond_names <- c("Heart disease",
                "Poorly controlled blood pressure", "Diabetes", "Chronic lung disease",
                paste0("Obesity (BMI ", "\U2265", "40)"), "Chronic kidney failure", "Chronic liver disease",
                "Cancer (Solid organ)", "Haematologial cancer", "Receiving chemotherapy or radiotherapy",
                "Organ transplant recipient", "Bone marrow transplant recipient", "Neurological condition",
                "Chronic inflammatory condition", "Immunodeficiency")

#load in data

N                         ## number of responders
N_strat                   ## number of strategies
N_sched                   ## number of schedules

s <- dat$strat            ## vaccine strategy
t <- dat$sched            ## schedule - 2mths, 4mths, 6mths, 12 mths
w <- dat$sex              ## sex
x <- dat$indig            ## Indigenous status
z <- dat$conds            ## comorbidity
y <- dat$AEFI             ## outcome - AEFI, MA or ?Fever


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
