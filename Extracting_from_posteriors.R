###Visualising and tabulating the posteriors - marginalised estimates
library(data.table)
library(ggplot2)
library(tidyverse)
library(stringr)
library(posterior)
postr <- readRDS(file = "C:/Users/ETay/Documents/Concom_AEFI_posterior_SA_eps.rds")
postr <- posterior::as_draws_matrix(fit_C$draws())    
#load in the relevant posterior
draws_full <- readRDS(file ="C:/Users/ETay/Documents/postr_concom_AEFI.rds")



#load in dat
dat <- fread(file = "C:/Users/ETay/Documents/NIP_MenB_dat.csv")
dat$res <- as.integer(factor(dat$uid_person, levels = unique(dat$uid_person)))
dat$vax_sequence[dat$vax_sequence == "2"] <- 0
dat$vax_sequence <- as.integer(dat$vax_sequence)


#marginalise - using option 2.	The event probability (e.g., one/two MAs) for an average person 
#from the population of survey responders receiving their X month schedule.

draws_marg <- array(NA, dim = c(nrow(draws_full), 2, 4, 2), 
                    dimnames = list(draw = 1:nrow(draws_full), 
                                    strategy = c("concom", "separate"), 
                                    schedule = c("2 months", "4 months", "6 months", "12 months"), 
                                    parameter = c("one", "two")))
sex_weights <- dat[, .(mn = mean(sex)), by = schedule]
indig_weights <- dat[, .(mn = mean(indig)), by = schedule]
state_weights <- data.table(mn = as.vector(prop.table(table(dat$schedule,dat$clinic_state), margin = 1)[,-2]),
                            schedule = rep(1:4, 7),
                            state = rep(1:7, each = 4))
clinic_weights <- data.table(mn = as.vector(prop.table(table(dat$schedule,dat$clinic_type), margin = 1)[,-2]),
                             schedule = rep(1:4, 2),
                             state = rep(1:2, each = 4))
comorb_weights <- dat[, .(mn = mean(pmh)), by = schedule]

# sex_weights <- dat[, .(mn = mean(sex))]
# indig_weights <- dat[, .(mn = mean(indig))]
# state_weights <- data.table(mn = as.vector(prop.table(table(dat$clinic_state))[-2]),
#                             state = 1:7)
# clinic_weights <- data.table(mn = as.vector(prop.table(table(dat$clinic_type))[-2]),
#                              state = 1:2)
# comorb_weights <- dat[, .(mn = mean(pmh))]


for(strat in 1:2){
  for(sched in 1:4){
    for(par in 1:2){
      mu_par <- draws_full[,paste0("mu[", sched, ",", par, "]")]
      alpha_par <- (strat == 2)*draws_full[,paste0("alpha[", sched, "]")]
      sex_par <- draws_full[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
      indig_par <- draws_full[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
      state_par <- as.vector(state_weights[schedule==sched]$mn%*%t(draws_full[,paste0("rho[", 1:7, ",", par, "]")]))
      clinic_par <- as.vector(clinic_weights[schedule==sched]$mn%*%t(draws_full[,paste0("tau[", 1:2, ",", par, "]")]))
      comorb_par <- draws_full[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
      if(strat == 2 & par == 2){
        draws_marg[, strat, sched, par] <- -Inf
      } else {
        draws_marg[, strat, sched, par] <- as.vector(mu_par + alpha_par + sex_par + indig_par + state_par + clinic_par + comorb_par)  
      }
    }
  }
}

dat_vis <- data.table(x = as.vector(draws_marg),
                      samp = rep(1:8000,2*4*2),
                      strategy = rep(rep(c("Separate", "Concomitant"), each = 8000), 4*2),
                      schedule = rep(rep(c("2 months", "4 months", "6 months", "12 months"), each = 8000*2), 2),
                      Parameter = rep(c("eta_1", "eta_2"), each = 8000*2*4))
                      
                      #Parameter = rep(c("P(k = 1)", "P(k = 2)"), each = 8000*2*4))
dat_vis <- dcast.data.table(dat_vis, samp + strategy + schedule ~ Parameter, value.var = 'x')
dat_vis[, `:=`(
  p1 = exp(eta_1)/(1 + exp(eta_1) + exp(eta_2)),
  p2 = exp(eta_2)/(1 + exp(eta_1) + exp(eta_2))
)]
dat_vis[, p12 := p1 + p2]
dat_vis <- melt.data.table(dat_vis, id.vars = c('samp', 'strategy', 'schedule'), measure.vars = "p12", variable.name = 'Parameter', value.name = 'x')
#dat_vis <- melt.data.table(dat_vis, id.vars = c('samp', 'strategy', 'schedule'), measure.vars = c('P(k = 1)', 'P(k = 2)', 'P(k >= 1)'), variable.name = 'Parameter', value.name = 'x')
#dat_vis <- dat_vis[!(strategy == 'Concomitant' & Parameter %in% c('P(k = 2)', 'P(k >= 1)'))]
#dat_vis <- dat_vis[!(strategy == 'Separate' & Parameter %in% c('P(k = 1)', 'P(k = 2)'))]

dat_vis[strategy == "Separate", mean(x), by = .(schedule, Parameter)]
dat_vis[strategy == "Concomitant", mean(x), by = .(schedule, Parameter)]

#contrasts - separate strategy at  least one event postr - concomitant 1 event postr 
#then summarise as another column in Tabble 3

#AEFI
ggplot(dat_vis, aes(x = x, colour = strategy)) +
  facet_grid(strategy~forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("Concomitant" = "#D55E00", "Separate" = "#009E73"),
                      labels = c("One report of AEFI", "At least one report of AEFI"))

# ggplot(dat_vis, aes(x = x, colour = Parameter)) +
#   facet_grid(strategy~forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")) +
#   geom_density() +
#   xlab("Probability") +
#   ylab("Density") +
#   labs(color = NULL) +
#   theme(legend.position = "bottom") +
#   scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k >= 1)" = "#009E73"),
#                       labels = c("One report of AEFI", "At least one report of AEFI"))

ggsave("AEFI_3.png", dpi = 400, width = 6, height = 5, units = "in")

## Figure 3 code

## probabilities and 95% credible intervals by schedule

summary_fun_C <- function(strat, sch, dat_vis){
  if(strat == "Concomitant"){
    dist <- dat_vis[strategy == "Concomitant" & schedule == sch & Parameter == "P(k = 1)",]$x
    mn <- format(round(mean(dist), 2), nsmall = 2)
    cr <- format(round(quantile(dist, c(0.025, 0.975)), 2), nsmall = 2)
    paste0("Concomitant strategy P(k = 1): ", mn, " (", cr[1], ", ", cr[2], ")")
  } else {
    dist1 <- dat_vis[strategy == "Separate" & schedule == sch & Parameter == "P(k = 1)",]$x
    dist2 <- dat_vis[strategy == "Separate" & schedule == sch & Parameter == "P(k = 2)",]$x
    dist_atleast <- dat_vis[strategy == "Separate" & schedule == sch & Parameter == "P(k = 1)",]$x + dat_vis[strategy == "Separate" & schedule == sch & Parameter == "P(k = 2)",]$x
    mn1 <- format(round(mean(dist1), 2), nsmall = 2)
    cr1 <- format(round(quantile(dist1, c(0.025, 0.975)), 2), nsmall = 2)
    mn2 <- format(round(mean(dist2), 2), nsmall = 2)
    cr2 <- format(round(quantile(dist2, c(0.025, 0.975)), 2), nsmall = 2)
    mn_atleast <- format(round(mean(dist_atleast), 2), nsmall = 2)
    cr_atleast <- format(round(quantile(dist_atleast, c(0.025, 0.975)), 2), nsmall = 2)
    paste0("Separate strategy P(k = 1): ", mn1, " (", cr1[1], ", ", cr1[2], ") ",
           "Separate strategy P(k = 2): ", mn2, " (", cr2[1], ", ", cr2[2], ") ",
           "Separate strategy P(k >= 1): ", mn_atleast, " (", cr_atleast[1], ", ", cr_atleast[2], ")")
  }
}

results1 <- list()
results2 <- list()
schedules <- c("2 months","4 months", "6 months", "12 months")

for (sch in schedules) {
  results1[sch] <- summary_fun_C("Concomitant", sch, dat_vis)
  results2[sch] <- summary_fun_C("Separate", sch, dat_vis)
}
print(results1)
print(results2)



sched <- 3 # 6 months
par <- 2 
strat <- 1
      mu_par <- draws_full[,paste0("mu[", sched, ",", par, "]")]
       alpha_par <- (strat == 2)*draws_full[,paste0("alpha[", sched, "]")]
       sex_par <- draws_full[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
       indig_par <- draws_full[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
       state_par <- as.vector(state_weights[schedule==sched]$mn%*%t(draws_full[,paste0("rho[", 1:7, ",", par, "]")]))
       clinic_par <- as.vector(clinic_weights[schedule==sched]$mn%*%t(draws_full[,paste0("tau[", 1:2, ",", par, "]")]))
       comorb_par <- draws_full[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
 summary(mu_par)


