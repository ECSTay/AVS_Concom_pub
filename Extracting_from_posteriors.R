###Visualising and tabulating the posteriors - marginalised estimates
library(data.table)
library(ggplot2)
library(tidyverse)
library(stringr)
library(posterior)
     
#load in the relevant posterior
draws_full <- readRDS(file = paste0("C:/Users/ETay/Documents/postr_concom_AEFI.rds"))



#load in dat
dat <- fread(file = "C:/Users/ETay/Documents/NIP_MenB_dat.csv")
dat$res <- as.integer(factor(dat$uid_person, levels = unique(dat$uid_person)))
dat$vax_sequence[dat$vax_sequence == "2"] <- 0
dat$vax_sequence <- as.integer(dat$vax_sequence)

# N_state  = 7                    
# N_clinic = 2
# 
# q <- cbind(as.integer(dat$clinic_state == "ACT"), #col1
#            as.integer(dat$clinic_state == "NT"),  #col2
#            as.integer(dat$clinic_state == "QLD"), #col3
#            as.integer(dat$clinic_state == "SA"),  #col4
#            as.integer(dat$clinic_state == "TAS"), #col5
#            as.integer(dat$clinic_state == "VIC"), #col6
#            as.integer(dat$clinic_state == "WA"))  #col7
# 
# c <- cbind(as.integer(dat$clinic_type == "Aboriginal Health Service"),
#            as.integer(dat$clinic_type == "State Health"))


#marginalise - using option 2.	The event probability (e.g., one/two MAs) for an average person from the population of survey responders receiving their X month schedule.

draws_marg <- array(NA, dim = c(nrow(draws_full), 2, 4, 2), 
                    dimnames = list(draw = 1:nrow(draws_full), 
                                    strategy = c("concom", "separate"), 
                                    schedule = c("2 months", "4 months", "6 months", "12 months"), 
                                    parameter = c("one", "two")))
sex_weights <- dat[, .(mn = mean(sex)), by = schedule]
indig_weights <- dat[, .(mn = mean(indig)), by = schedule]

#state_weights
a <- table(dat$schedule,dat$clinic_state)
a <- a[, -2]#excludeNSW
b <- prop.table(a, margin = 2)
colnames(b)[c(1:7)] <- c("1", "2", "3","4","5","6","7")#ACT,NT,QLD,SA,TAS,VIC,WA
state_weights <- as.data.table(b)
setnames(state_weights, old = "V1", new = "schedule")
setnames(state_weights, old = "V2", new = "state")
state_weights$schedule <- as.integer(state_weights$schedule)
state_weights$state <- as.integer(state_weights$state)
setnames(state_weights, old = "N", new = "mn")

#clinic_weights
c <- table(dat$schedule,dat$clinic_type)
c <- c[, -2]#excludeGP
d <- prop.table(c, margin = 2)
colnames(d)[c(1:2)] <- c("1", "2")#AHS,State
clinic_weights <- as.data.table(d)
setnames(clinic_weights, old = "V1", new = "schedule")
setnames(clinic_weights, old = "V2", new = "clinic")
clinic_weights$schedule <- as.integer(clinic_weights$schedule)
clinic_weights$clinic <- as.integer(clinic_weights$clinic)
setnames(clinic_weights, old = "N", new = "mn")

comorb_weights <- dat[, .(mn = mean(pmh)), by = schedule]

as.vector(state_weights[schedule==sched]$mn%*%t(draws_full[,paste0("rho[", state, ",", strat, "]")]))

for(strat in 1:2){
  for(sched in 1:4){
    for(state in 1:7){
      for(clinic in 1:2){
        for(par in 1:2){
          sex_par <- draws_full[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
          indig_par <- draws_full[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
          state_par <- as.vector(state_weights[schedule==sched]$mn%*%t(draws_full[,paste0("rho[", state, ",", strat, "]")]))
          clinic_par <- as.vector(clinic_weights[schedule==sched]$mn%*%t(draws_full[,paste0("tau[", clinic, ",", strat, "]")]))
          comorb_par <- draws_full[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
          mu_par <- draws_full[,paste0("mu[", strat, ",", sched, ",", par, "]")]
          draws_marg[, strat, sched, par] <- as.vector(mu_par + sex_par + indig_par + state_par + clinic_par + comorb_par)####
        } 
      }
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

#contrasts - separate strategy 2 events postr - concomitant 1 event postr then summarise as another column in Tabble 3

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






