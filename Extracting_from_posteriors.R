###Visualising and tabulating the posteriors

#load in the relevant posterior

draws_full <- 
#marginalise - using option 2.	The event probability (e.g., one/two MAs) for an average person from the population of survey responders receiving their X month schedule.

draws_marg <- array(NA, dim = c(nrow(draws_full), 2, 4, 2), 
                    dimnames = list(draw = 1:nrow(draws_full), 
                                    strategy = c("concom", "separate"), 
                                    schedule = c("2 months", "4 months", "6 months", "12 months"), 
                                    parameter = c("one", "two")))
sex_weights <- dat[, .(mn = mean(sex)), by = schedule]
indig_weights <- dat[, .(mn = mean(indig)), by = schedule]
comorb_weights <- dat[, .(mn = mean(pmh)), by = schedule]

for(strat in 1:2){
  for(sched in 1:4){
    for(par in 1:2){
      sex_par <- draws_full[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
      indig_par <- draws_full[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
      comorb_par <- draws_full[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
      mu_par <- draws_full[,paste0("mu[", strat, ",", sched, ",", par, "]")]
      draws_marg[, strat, sched, par] <- as.vector(mu_par + sex_par + indig_par + comorb_par)
    }
  }
}

dat_vis <- data.table(x = plogis(as.vector(draws_marg)),
                      samp = rep(1:8000,2*4*2),
                      strategy = rep(rep(c("Concomitant", "Separate"), each = 8000), 4*2),
                      schedule = rep(rep(c("2 months", "4 months", "6 months", "12 months"), each = 8000*2), 2),
                      Parameter = rep(c("P(k = 1)", "P(k = 2)"), each = 8000*2*4))
dat_vis <- dat_vis[!(strategy == "Concomitant" & Parameter == "P(k = 2)")]

dat_vis <- dcast.data.table(dat_vis, samp + strategy + schedule ~ Parameter, value.var = 'x')
#dat_vis[, 'P(k >= 1)' := 'P(k = 1)' + 'P(k = 2)']
dat_vis$"P(k >= 1)" <- dat_vis$'P(k = 1)' + dat_vis$'P(k = 2)'
dat_vis <- melt.data.table(dat_vis, id.vars = c('samp', 'strategy', 'schedule'), measure.vars = c('P(k = 1)', 'P(k = 2)', 'P(k >= 1)'), variable.name = 'Parameter', value.name = 'x')
dat_vis <- dat_vis[!(strategy == 'Concomitant' & Parameter %in% c('P(k = 2)', 'P(k >= 1)'))]
dat_vis <- dat_vis[!(strategy == 'Separate' & Parameter %in% c('P(k = 1)', 'P(k = 2)'))]

#AEFI
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(strategy~forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k >= 1)" = "#009E73"),
                      labels = c("One report of AEFI", "At least one report of AEFI"))

ggsave("AEFI_3.png", dpi = 400, width = 6, height = 5, units = "in")

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
