## posterior predictive distributions

library(ggplot2)
library(ggdist)
library(distributional)

## draws is the posterior distribution
## dat is the stan data set you used (it should be a list with N_R, N_I, N_sched, etc.)
start_time <- Sys.time()
posterior_predictive <- function(draws, dat){
  ndraws <- nrow(draws)
  y_pred <- array(data = NA_integer_, dim = c(ndraws, dat$N_R))
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
    y_pred[,r] <- apply(p, 1, function(p_sim) sample.int(3, size = 1, prob = p_sim) - 1)
  }

  sim_df <- data.table(sim = rep(1:ndraws, dat$N_R),#can add eg sex
                       y   = as.vector(y_pred),
                       s   = factor(rep(dat$s, each = ndraws), labels = c("Separate", "Concomitant")),
                       t   = factor(rep(dat$t, each = ndraws), labels = c("2 months", "4 months", "6 months", "12 months")))
  sim_summ <- sim_df[, .(`0` = mean(y == 0), `1` = mean(y == 1), `2` = mean(y == 2)), by = .(sim, s, t)]##avge over the covariates
  sim_summ <- melt.data.table(sim_summ, id.vars = c("sim", "s", "t"), measure.vars = c("0", "1", "2"), variable.name = "Outcome")
  sim_summ <- sim_summ[!(s == "Concomitant" & Outcome == "2"), .(s,t,Outcome, value)]
  sim_summ <- sim_summ[, .(value = list(value)), by = .(s, t, Outcome)]
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
  dat_vis <- merge(sim_summ,dat_raw,by = c("s","t","Outcome"))
  p <- ggplot(dat_vis, aes(y = Outcome, xdist = dist_sample(value))) +
    facet_grid(t ~ s) +
    stat_histinterval(aes(fill = after_stat(level)), .width = c(.50, .95, 1)) +
    stat_pointinterval(.width = c(.50, .95, 1)) +
    scale_fill_brewer(breaks = c(1, .95, .50),name = "Proportion") +
    geom_point(aes(x = y_mean), colour = "red", shape = 17) +
    xlab("Probability")
  return(p)
}

posterior_predictive(draws, dat)
time <- Sys.time() - start_time
## SNAP code below if needed

ppc_dat <- as_tibble(dat_stan)
ppc_dat <- bind_cols(ppc_dat, y_rep)
ppc <- as_tibble(ppc_dat) |>
    group_by(s, t) |>
    summarise(y = mean(y), y_mid = median(y), y_rep = rvar_mean(value), n = n())
ggplot(ppc, aes(y = 1:nrow(ppc), xdist = y_rep)) +
  stat_histinterval(aes(fill = after_stat(level)), .width = c(.50, .95, 1)) +
  stat_pointinterval(.width = c(.50, .95, 1)) +
  scale_fill_brewer(breaks = c(1, .95, .50)) +
  geom_point(aes(x = y, y = 1:nrow(ppc)), colour = "red", shape = 23)

ggsave("PPD_NIPMenB.png", dpi = 400, width = 10, height = 8, units = "in")
