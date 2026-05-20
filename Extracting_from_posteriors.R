###Visualising the posteriors - marginalised estimates
library(data.table)
library(ggplot2)
library(ggdist)
library(distributional)
library(tidyverse)
library(stringr)
library(posterior)

#load in the relevant posterior
draws <- readRDS(file ="C:/Users/ETay/Documents/postr_concom_impact.rds")
#draws <- readRDS(file ="C:/Users/ETay/Documents/postr_concom_ma_SA.rds")# SA with N(0),10) priors

#load in the STAN dat
dat <- readRDS("C:/Users/ETay/Documents/dat_C_impact.rds")
#marginalise - using option 2.	The event probability (e.g., one/two MAs) for an average person 
#from the population of survey responders receiving their X month schedule.

ndraws <- nrow(draws)
ninds <- 1000
p_post <- array(data = NA, dim = c(ndraws, 4, 2, 4))
for(t in 1:4){
  inds <- sample(which(dat$t == t), size = ninds, replace = TRUE)
  mu1 <- draws[1:ndraws,paste0("mu[", t, ",1]")]
  mu2 <- draws[1:ndraws,paste0("mu[", t, ",2]")]
  alpha <- draws[1:ndraws,paste0("alpha[", t, "]")]
  epsilon1 <- sapply(1:ninds, function(r) rnorm(ndraws, sd = draws[1:ndraws,"sigma_epsilon[1]"]))
  epsilon2 <- sapply(1:ninds, function(r) rnorm(ndraws, sd = draws[1:ndraws,"sigma_epsilon[2]"]))
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
  for(s in 1:2){
    p_tmp <- array(NA, dim = c(ndraws, 4, ninds))
    for(r in 1:ninds){
      eta <- array(NA, dim = c(ndraws, 3))
      eta[,1] <- 0
      eta[,2] <- mu1 + (s == 2)*alpha + epsilon1[,r] + dat$w[inds[r]]*beta1 + dat$x[inds[r]]*gamma1 + as.vector(dat$q[inds[r],]%*%t(rho1)) + as.vector(dat$c[inds[r],]%*%t(tau1)) + dat$z[inds[r]]*delta1
      if(s == 1){
        eta[,3] <- mu2 + epsilon2[,r] + dat$w[inds[r]]*beta2 + dat$x[inds[r]]*gamma2 + as.vector(dat$q[inds[r],]%*%t(rho2)) + as.vector(dat$c[inds[r],]%*%t(tau2)) + dat$z[inds[r]]*delta2
      } else {
        eta[,3] <- - Inf
      }
      p_tmp[,1:3,r] <- exp(eta - apply(eta, 1, function(x) max(x) + log(sum(exp(x - max(x))))))
      p_tmp[,4,r] <- p_tmp[,2,r] + p_tmp[,3,r]
    }
    p_post[,,s,t] <- apply(p_tmp, c(1,2), mean)
  }
}
  
dat_vis <- data.table(sim = rep(1:ndraws, 4*2*4),
                      p_marg   = as.vector(p_post),
                      Outcome = factor(rep(rep(c("P(0 events)", "P(1 event)", "P(2 events)", "P(At least 1 event)"), each = ndraws), 2*4), levels = c("P(0 events)", "P(1 event)", "P(2 events)", "P(At least 1 event)")),
                      s   = factor(rep(rep(c("Separate", "Concomitant"), each = ndraws*4), 4), levels = c("Separate", "Concomitant")),
                      t   = factor(rep(c("2 months", "4 months", "6 months", "12 months"), each = ndraws*4*2), levels = c("2 months", "4 months", "6 months", "12 months")))
  
#Example code for the Probability distribution figures in the Supplementary Information
ggplot(dat_vis[Outcome == "P(At least 1 event)"], aes(x = p_marg, colour = s)) +
  facet_grid(s~forcats::fct_relevel(t, "2 months", "4 months", "6 months", "12 months")) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("Concomitant" = "#D55E00", "Separate" = "#009E73"),
                      labels = c("One report of AEFI", "At least one report of AEFI"))

#code to derive the data files (dat_vis_cred) for Figure 3
summ_f <- function(x) paste0(round(mean(x), 3), " (", round(quantile(x, 0.025), 3), ",", round(quantile(x, 0.975), 3), ")")
dat_vis_cred <- dat_vis[, summ_f(p_marg), by = .(Outcome, s, t)]


#contrast values
dat_vis_cont <- dcast.data.table(dat_vis[Outcome == "P(At least 1 event)"], sim + t ~ s, value.var = "p_marg")
dat_vis_cont[, contrast := Concomitant - Separate] 

#for the probability distributions of the "contrast" figures
  
ggplot(dat_vis_cont, aes(x = contrast)) +
    facet_grid(~forcats::fct_relevel(t, "2 months", "4 months", "6 months", "12 months")~.) +
    geom_density() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab("Change in Probability") +
    ylab("Density")

ggplot2::ggsave("AEFI_impact.png", dpi = 400, width = 8, height = 4, units = "in")
  
dat_vis_cont[, summ_f(contrast), by = t]

#######################################################################
