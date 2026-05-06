##extracting from ...eps

library(data.table)
library(ggplot2)
library(tidyverse)
library(stringr)
library(posterior)
postr <- readRDS(file = "C:/Users/ETay/Documents/Concom_AEFI_posterior_SA_eps.rds")
postr <- posterior::as_draws_matrix(postr) 

draws_full <- postr[,str_detect(colnames(postr), "mu|alpha|beta|delta|rho|tau|gamma|sigma")]


draws_marg <- array(NA, dim = c(nrow(draws_full), 2, 4, 2), 
                     dimnames = list(draw = 1:nrow(draws_full), 
                     strategy = c("concom", "separate"), 
                     schedule = c("2 months", "4 months", "6 months", "12 months"), 
                     parameter = c("one", "two")))
# sex_weights <- dat[, .(mn = mean(sex)), by = schedule]
# indig_weights <- dat[, .(mn = mean(indig)), by = schedule]
# state_weights <- data.table(mn = as.vector(prop.table(table(dat$schedule,dat$clinic_state), margin = 1)[,-2]),
#                                            schedule = rep(1:4, 7),
#                                            state = rep(1:7, each = 4))
# clinic_weights <- data.table(mn = as.vector(prop.table(table(dat$schedule,dat$clinic_type), margin = 1)[,-2]),
#                                            schedule = rep(1:4, 2),
#                                            state = rep(1:2, each = 4))
# comorb_weights <- dat[, .(mn = mean(pmh)), by = schedule]
for(strat in 1:2){
     for(sched in 1:4){
         for(par in 1:2){
             mu_par <- draws_full[,paste0("mu[", sched, ",", par, "]")]
             alpha_par <- (strat == 2)*draws_full[,paste0("alpha[", sched, "]")]
             #sex_par <- draws_full[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
             #indig_par <- draws_full[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
             #state_par <- as.vector(state_weights[schedule==sched]$mn%*%t(draws_full[,paste0("rho[", 1:7, ",", par, "]")]))
            # clinic_par <- as.vector(clinic_weights[schedule==sched]$mn%*%t(draws_full[,paste0("tau[", 1:2, ",", par, "]")]))
             #comorb_par <- draws_full[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
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
dat_vis


dat_vis <- dcast.data.table(dat_vis, samp + strategy + schedule ~ Parameter, value.var = 'x')

dat_vis[, `:=`(
    p1 = exp(eta_1)/(1 + exp(eta_1) + exp(eta_2)),
     p2 = exp(eta_2)/(1 + exp(eta_1) + exp(eta_2))
   )]

dat_vis[, p12 := p1 + p2]
dat_vis[,.(mean(p1),mean(p2)), by = .(strategy, schedule)][order(schedule)]

dat_vis <- melt.data.table(dat_vis, id.vars = c('samp', 'strategy', 'schedule'), measure.vars = "p12", variable.name = 'Parameter', value.name = 'x')

