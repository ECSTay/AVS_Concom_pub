##############Model C analysis script - prob of reporting at least one AEFI/MA following either NIP or Men B alone AND concom

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
z <- dat$PMH            ## comorbidity
y <- dat$AEFI             ## outcome - AEFI, MA or ?Fever


a <-  qlogis(0.3)          ## prior distribution mean for one event - depends on the schedule and the vaccine strategy?
b <-                       ## prior distribution standard deviation for two events
c <- qlogis(0.05)          ## prior distribution mean for two events
d <-                       ## prior distribution for standard deviation for two events

  dat_C <- list(N = N_C,
                N_strat = N_strat_C,
                N_sched = N_sched_C,
                s = s_C,
                t = t_C,
                w = w_C,
                x = x_C,
                z = z_C,
                y = y_C,
                a = a,
                b = b,
                c = c,
                d = d)



model_C <- cmdstan_model(write_stan_file(readLines("C:/Users/etay/Documents/Work documents/AVS work/Thuy_concom/model_C.stan")))###

fit_C <- model_C$sample(data = dat_C, 
                  chains = 8, parallel_chains = 8)###

postr <- posterior::as_draws_matrix(fit_C$draws())

##SAP code
model_C <- cmdstan_model("C:/Users/etay/Documents/Work documents/AVS work/Thuy_concom/model_C.stan")

drp <- capture.output({fit_C <- model_C$sample(dat_C, chains = 8, parallel_chains = 8)})
draws_C <- as_draws_matrix(fit_C$draws(c("mu")))
dat_vis_C <- data.frame(x = plogis(as.vector(draws_C)),
                        strategy = rep(rep(c("Concomitant", "Separate"), each = 8000), 2),
                        Parameter = rep(c("P(k = 1)", "P(k = 2)"), each = 8000*2))
dat_vis_C <- dat_vis_C[!(dat_vis_C$strategy == "Concomitant" & dat_vis_C$Parameter == "P(k = 2)"),]

summary_fun_C <- function(strat, k, dat){
  mn <- format(round(mean(dat[dat$strategy == strat & dat$Parameter == paste0("P(k = ", k, ")"),]$x), 2), nsmall = 2)
  cr <- format(round(quantile(dat[dat$strategy == strat & dat$Parameter == paste0("P(k = ", k, ")"),]$x, c(0.025, 0.975)), 2), nsmall = 2)
  paste0(strat, " strategy P(k = ", k, "): ", mn, " (", cr[1], ", ", cr[2], ")")
}

summary_fun_C("Concomitant", 1, dat_vis_C)
summary_fun_C("Separate", 1, dat_vis_C)
summary_fun_C("Separate", 2, dat_vis_C)
