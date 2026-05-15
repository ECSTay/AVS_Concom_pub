###Visualising and tabulating the posteriors - marginalised estimates
library(data.table)
library(ggplot2)
library(ggdist)
library(distributional)
library(tidyverse)
library(stringr)
library(posterior)

#load in the relevant posterior
draws <- readRDS(file ="C:/Users/ETay/Documents/postr_concom_AEFI.rds")
#draws <- readRDS(file ="C:/Users/ETay/Documents/postr_concom_ma_SA.rds")# SA with N(0),10) priors

#load in dat

dat <- readRDS("C:/Users/ETay/Documents/dat_C_AEFI.rds")
#marginalise - using option 2.	The event probability (e.g., one/two MAs) for an average person 
#from the population of survey responders receiving their X month schedule.
#extracting_posteriors <- function(draws, dat){
  ndraws <- nrow(draws)
  p_post <- array(data = NA_integer_, dim = c(ndraws, dat$N_R))
  for(r in 1:dat$N_R){
    mu1 <- draws[1:ndraws,paste0("mu[", dat$t[r], ",1]")]
    mu2 <- draws[1:ndraws,paste0("mu[", dat$t[r], ",2]")]
    alpha <- draws[1:ndraws,paste0("alpha[", dat$t[r], "]")]
    epsilon1 <- rnorm(ndraws, sd = draws[1:ndraws,"sigma_epsilon[1]"])
    epsilon2 <- rnorm(ndraws, sd = draws[1:ndraws,"sigma_epsilon[2]"])
    beta1 <- draws[1:ndraws,"beta[1]"]
    beta2 <- draws[1:ndraws,"beta[2]"]
    gamma1 <- draws[1:ndraws,"gamma[1]"]
    gamma2 <- draws[1:ndraws,"gamma[2]"]
    rho1 <- draws[1:ndraws,paste0("rho[", 1:dat$N_state, ",1]")]
    rho2 <- draws[1:ndraws,paste0("rho[", 1:dat$N_state, ",2]")]
    tau1 <- draws[1:ndraws,paste0("tau[", 1:dat$N_clinic, ",1]")]
    tau2 <- draws[1:ndraws,paste0("tau[", 1:dat$N_clinic, ",2]")]
    delta1 <- draws[1:ndraws,"delta[1]"]
    delta2 <- draws[1:ndraws,"delta[2]"]
    eta <- array(NA, dim = c(ndraws, 3))
    eta[,1] <- 0
    eta[,2] <- mu1 + dat$s[r]*alpha + epsilon1 + dat$w[r]*beta1 + dat$x[r]*gamma1 + as.vector(dat$q[r,]%*%t(rho1)) + as.vector(dat$c[r,]%*%t(tau1)) + dat$z[r]*delta1
    if(dat$s[r] == 0){
      eta[,3] <- mu2 + epsilon2 + dat$w[r]*beta2 + dat$x[r]*gamma2 + as.vector(dat$q[r,]%*%t(rho2)) + as.vector(dat$c[r,]%*%t(tau2)) + dat$z[r]*delta2
    } else {
      eta[,3] <- - Inf
    }
    p <- exp(eta - apply(eta, 1, function(x) max(x) + log(sum(exp(x - max(x))))))
    p_post[,r] <- apply(p, 1, function(p_sim) sample.int(3, size = 1, prob = p_sim) - 1)
  }
  
  dat_vis <- data.table(sim = rep(1:ndraws, dat$N_R),#can add eg sex
                       y   = as.vector(p_post),
                       s   = factor(rep(dat$s, each = ndraws), labels = c("Separate", "Concomitant")),
                       t   = factor(rep(dat$t, each = ndraws), labels = c("2 months", "4 months", "6 months", "12 months")))
  dat_vis <- dat_vis[, .(`0` = mean(y == 0), `1` = mean(y == 1), `2` = mean(y == 2)), by = .(sim, s, t)]##avge over the covariates
  
  dat_vis <- melt.data.table(dat_vis, id.vars = c("sim", "s", "t"), measure.vars = c("0", "1", "2"), variable.name = "Outcome")
  dat_vis <- dat_vis[!(s == "Concomitant" & Outcome == "2"), .(s,t,Outcome, value)]
  dat_vis <- dat_vis[, .(value = list(value)), by = .(s, t, Outcome)]
  
  dat_raw <- data.table(s = factor(rep(c("Separate", "Concomitant"), each = 4*3), levels = c("Separate", "Concomitant")),
                        t = factor(rep(rep(c("2 months", "4 months", "6 months", "12 months"), each = 3), 2), levels = c("2 months", "4 months", "6 months", "12 months")),
                        Outcome = factor(rep(0:2, 2*4)),
                        y_mean = c(mean(dat$y[dat$s == 0 & dat$t == 1] == 1), mean(dat$y[dat$s == 0 & dat$t == 1] == 2), mean(dat$y[dat$s == 0 & dat$t == 1] == 3),
                                   mean(dat$y[dat$s == 0 & dat$t == 2] == 1), mean(dat$y[dat$s == 0 & dat$t == 2] == 2), mean(dat$y[dat$s == 0 & dat$t == 2] == 3),
                                   mean(dat$y[dat$s == 0 & dat$t == 3] == 1), mean(dat$y[dat$s == 0 & dat$t == 3] == 2), mean(dat$y[dat$s == 0 & dat$t == 3] == 3),
                                   mean(dat$y[dat$s == 0 & dat$t == 4] == 1), mean(dat$y[dat$s == 0 & dat$t == 4] == 2), mean(dat$y[dat$s == 0 & dat$t == 4] == 3),
                                   mean(dat$y[dat$s == 1 & dat$t == 1] == 1), mean(dat$y[dat$s == 1 & dat$t == 1] == 2), mean(dat$y[dat$s == 1 & dat$t == 1] == 3),
                                   mean(dat$y[dat$s == 1 & dat$t == 2] == 1), mean(dat$y[dat$s == 1 & dat$t == 2] == 2), mean(dat$y[dat$s == 1 & dat$t == 2] == 3),
                                   mean(dat$y[dat$s == 1 & dat$t == 3] == 1), mean(dat$y[dat$s == 1 & dat$t == 3] == 2), mean(dat$y[dat$s == 1 & dat$t == 3] == 3),
                                   mean(dat$y[dat$s == 1 & dat$t == 4] == 1), mean(dat$y[dat$s == 1 & dat$t == 4] == 2), mean(dat$y[dat$s == 1 & dat$t == 4] == 3)))
  
  dat_raw <- dat_raw[!(s == "Concomitant" & Outcome == "2")]
  dat_vis <- merge(dat_vis,dat_raw,by = c("s","t","Outcome"))
  dat_vis
  
  quantile1 <- vector()
  for(i in 1:ndraws){
  quantile1[i] <- quantile(dat_vis$value[[i]], prob = 0.025)
  }
  
  
  
  
  
  p <- ggplot(dat_vis, aes(y = Outcome, xdist = dist_sample(value))) +
    facet_grid(t ~ s) +
    stat_histinterval(aes(fill = after_stat(level)), .width = c(.50, .95, 1)) +
    stat_pointinterval(.width = c(.50, .95, 1)) +
    scale_fill_brewer(breaks = c(1, .95, .50),name = "Proportion") +
    geom_point(aes(x = y_mean), colour = "red", shape = 17) +
    xlab("Probability") +
    ylab("Number of reported adverse events")
  return(p)


########################################################################

draws_marg <- array(NA, dim = c(nrow(draws), 2, 4, 2), 
                    dimnames = list(draw = 1:nrow(draws), 
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


for(strat in 1:2){
  for(sched in 1:4){
    for(par in 1:2){
      mu_par <- draws[,paste0("mu[", sched, ",", par, "]")]
      alpha_par <- (strat == 2)*draws[,paste0("alpha[", sched, "]")]
      sex_par <- draws[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
      indig_par <- draws[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
      state_par <- as.vector(state_weights[schedule==sched]$mn%*%t(draws[,paste0("rho[", 1:7, ",", par, "]")]))
      clinic_par <- as.vector(clinic_weights[schedule==sched]$mn%*%t(draws[,paste0("tau[", 1:2, ",", par, "]")]))
      comorb_par <- draws[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
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


dat_vis[, `:=`(
  p1 = exp(eta[,2] - apply(eta, 1, function(x) max(x) + log(sum(exp(x - max(x)))))),
  p2 = exp(eta[,3] - apply(eta, 1, function(x) max(x) + log(sum(exp(x - max(x))))))
)]
dat_vis[, p12 := p1 + p2]
dat_vis[,.(mean(p1),mean(p2)), by = .(strategy, schedule)][order(schedule)]
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



## probabilities and 95% credible intervals by schedule for Table 3 - add a contrast column

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
      mu_par <- draws[,paste0("mu[", sched, ",", par, "]")]
       alpha_par <- (strat == 2)*draws[,paste0("alpha[", sched, "]")]
       sex_par <- draws[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
       indig_par <- draws[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
       state_par <- as.vector(state_weights[schedule==sched]$mn%*%t(draws[,paste0("rho[", 1:7, ",", par, "]")]))
       clinic_par <- as.vector(clinic_weights[schedule==sched]$mn%*%t(draws[,paste0("tau[", 1:2, ",", par, "]")]))
       comorb_par <- draws[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
 summary(mu_par)


