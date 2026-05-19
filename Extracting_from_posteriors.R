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

#load in STAN dat

dat <- readRDS("C:/Users/ETay/Documents/dat_C_AEFI.rds")
#marginalise - using option 2.	The event probability (e.g., one/two MAs) for an average person 
#from the population of survey responders receiving their X month schedule.

ndraws <- nrow(draws)
p_post <- array(data = NA, dim = c(ndraws, 3, dat$N_R))
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
    p_post[,,r] <- exp(eta - apply(eta, 1, function(x) max(x) + log(sum(exp(x - max(x))))))
    
}
  
dat_vis <- data.table(sim = rep(1:ndraws, 3*dat$N_R),#can add eg sex
                       p   = as.vector(p_post),
                       Outcome = rep(rep(c("0", "1", "2"), each = ndraws), dat$N_R),
                       s   = factor(rep(dat$s, each = ndraws*3), labels = c("Separate", "Concomitant")),
                       t   = factor(rep(dat$t, each = ndraws*3), labels = c("2 months", "4 months", "6 months", "12 months")))
  dat_vis <- dat_vis[Outcome != "0", .(p_marg = mean(p)), by = .(sim, Outcome, s, t)]##avge over the covariates
  dat_vis <- dcast.data.table(dat_vis, sim + s + t ~ Outcome, value.var = "p_marg")
  dat_vis[, `1+` := `1` + `2`]
  dat_vis <- melt.data.table(dat_vis, id.vars = c("sim", "s", "t"), measure.vars = c("1", "2", "1+"), variable.name = "Outcome", value.name = "p_marg")
  dat_vis[, Outcome := factor(Outcome, levels = c("1", "2", "1+"), labels = c("P(1 event)", "P(2 events)", "P(at least one event)"))]
  
  ggplot(dat_vis[Outcome == "P(at least one event)"], aes(x = p_marg, colour = s)) +
    facet_grid(s~forcats::fct_relevel(t, "2 months", "4 months", "6 months", "12 months")) +
    geom_density() +
    xlab("Probability") +
    ylab("Density") +
    labs(color = NULL) +
    theme(legend.position = "bottom") +
    scale_colour_manual(values = c("Concomitant" = "#D55E00", "Separate" = "#009E73"),
                        labels = c("One report of AEFI", "At least one report of AEFI"))
  
  summ_f <- function(x) paste0(round(median(x), 2), " (", round(quantile(x, 0.025), 2), ",", round(quantile(x, 0.975), 2), ")")
  dat_vis[, summ_f(p_marg), by = .(Outcome, s, t)]

  dat_vis_cont <- dcast.data.table(dat_vis[Outcome == "P(at least one event)"], sim + t ~ s, value.var = "p_marg")
  dat_vis_cont[, contrast := Concomitant - Separate] 
  
  ggplot(dat_vis_cont, aes(x = contrast)) +
    facet_grid(~forcats::fct_relevel(t, "2 months", "4 months", "6 months", "12 months")~.) +
    geom_density() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab("Change in Probability") +
    ylab("Density")
  
  dat_vis_cont[, summ_f(contrast), by = t]

########################################################################


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



