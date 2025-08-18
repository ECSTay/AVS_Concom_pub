##############Model C analysis script - prob of reporting at least one AEFI/Impact following (1) concom (2) NIP first (3) Men B first
library(cmdstanr)
library(posterior)
library(ggplot2)
library(forcats)


#load in data
dat <- read.csv(file = "C:/Users/ETay/Documents/Work documents/AVS work/Thuy_concom/AVS_Concom_pub/dat_modelC.csv")



N_C = nrow(dat)                   ## number of responders
N_strat_C  = 2                    ## number of strategies
N_sched_C  = 4                    ## number of schedules

s_C <- dat$vax_sequence           ## vaccine strategy, 1 = "Concomitant vaccination", 2 = "Separate"
#s_C <- dat$group                 ## vaccine strategy 1 = "Concomitant vaccination", 2 = "Separate"
t_C <- dat$sched                  ## schedule - 1 = 2 months, 2 = 4 months, 3 = 6 months, 4 = 12 months
w_C <- dat$sex                    ## sex - 0 = "Male", 1 = "Female"
x_C <- dat$indig                  ## Indigenous status - 0 = Non-Indig, 1 = Aboriginal and Torres Strait Islander
z_C <- dat$pmh                    ## comorbidity - 0 - None, 1 - at least one

#y_C <- dat$any_event + 1           ## outcome - 0,1,2 p(1) = 0.49, p(2) = 0.03
y_C <- dat$impact + 1              ## outcome - 0,1,2 p(1) = 0.04, p(2) = 0.0005
#y_C <- dat$medical_attention + 1                  ## outcome - 0,1,2 p(1) = 0.02
#y_C <- dat$local + 1               ## outcome - 0,1,2 p(1) = 0.27, p(2) = 0.02
#y_C <- dat$fever + 1               ## outcome - 0,1,2 p(1) = 0.27, p(2) = 0.01

#a = -3, b = 1  for P(impact), P(MA)
#a = -1, b = 1  for P(AEFI), P(local), P(fever)

a <- -3                       ## prior distribution mean for one event - depends on the schedule and the vaccine strategy
b <- 1                         ## prior distribution standard deviation for two events
c <- -1                        ## prior distribution mean for two events
d <- 1                         ## prior distribution for standard deviation for two events

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


##SAP code
model_C <- cmdstan_model("C:/Users/etay/Documents/Work documents/AVS work/Thuy_concom/AVS_Concom_pub/model_C.stan")

fit_C <- model_C$sample(dat_C, chains = 8, parallel_chains = 8) 
draws_full <- as_draws_matrix(fit_C$draws(c("mu", "beta", "gamma", "delta")))
#saveRDS(draws_full, "draws_full_any_event_2025-08-12.rds") # any AE
saveRDS(draws_full, "draws_full_impact_2025-08-12.rds") # any days of impact
#saveRDS(draws_full, "draws_full_ma_2025-08-12.rds") # MA
#saveRDS(draws_full, "draws_full_local_2025-08-12.rds") # local reaction
#saveRDS(draws_full, "draws_full_fever_2025-08-12.rds") # fever

#load in the relevant draws_full
## marginalise - using option b) - 2.	The event probability (e.g., one/two MAs) for an average person from the population of survey responders receiving their X month schedule.

draws_marg <- array(NA, dim = c(nrow(draws_full), 2, 4, 2), 
                    dimnames = list(draw = 1:nrow(draws_full), 
                                    strategy = c("concom", "separate"), 
                                    schedule = c("2 months", "4 months", "6 months", "12 months"), 
                                    parameter = c("one", "two")))
sex_weight <- mean(dat$sex == 1)
indig_weight <- mean(dat$indig == 1)
comorb_weight <- mean(dat$pmh == 1)


for(strat in 1:2){
  for(sched in 1:4){
    for(par in 1:2){
      sex_par <- draws_full[,paste0("beta[", par, "]")]*sex_weight
      indig_par <- draws_full[,paste0("gamma[", par, "]")]*indig_weight
      comorb_par <- draws_full[,paste0("delta[", par, "]")]*comorb_weight
      mu_par <- draws_full[,paste0("mu[", strat, ",", sched, ",", par, "]")]
      draws_marg[, strat, sched, par] <- as.vector(mu_par + sex_par + indig_par + comorb_par)
    }
  }
}

dat_vis <- data.frame(x = plogis(as.vector(draws_marg)),
                      strategy = rep(rep(c("Concomitant", "Separate"), each = 8000), 4*2),
                      schedule = rep(rep(c("2 months", "4 months", "6 months", "12 months"), each = 8000*2), 2),
                      Parameter = rep(c("P(k = 1)", "P(k = 2)"), each = 8000*2*4))
dat_vis <- dat_vis[!(dat_vis$strategy == "Concomitant" & dat_vis$Parameter == "P(k = 2)"),]

#AEFI
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")~strategy) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k = 2)" = "#009E73"),
                      labels = c("One AEFI", "Two AEFI"))

ggsave("AEFI.png", dpi = 400, width = 6, height = 5, units = "in")

#local reaction

ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")~strategy) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k = 2)" = "#009E73"),
                      labels = c("One local reaction", "Two local reactions"))

ggsave("Local.png", dpi = 400, width = 6, height = 5, units = "in")
#fever
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")~strategy) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k = 2)" = "#009E73"),
                      labels = c("One report of fever", "Two reports of fever"))

ggsave("Fever.png", dpi = 400, width = 6, height = 5, units = "in")
#MA
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")~strategy) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k = 2)" = "#009E73"),
                      labels = c("One report of MA", "Two reports of MA"))

ggsave("MA.png", dpi = 400, width = 6, height = 5, units = "in")
#impact
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")~strategy) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k = 2)" = "#009E73"),
                      labels = c("One report of carer impact", "Two reports of carer impact"))

ggsave("Impact.png", dpi = 400, width = 6, height = 5, units = "in")

results <- c()

for (sch in 1:schedules)

summary_fun_C <- function(strat, k, dat){
  mn <- format(round(mean(dat[dat$strategy == strat & dat$Parameter == paste0("P(k = ", k, ")"),]$x), 2), nsmall = 2)
  cr <- format(round(quantile(dat[dat$strategy == strat & dat$Parameter == paste0("P(k = ", k, ")"),]$x, c(0.025, 0.975)), 2), nsmall = 2)
  paste0(strat, " strategy P(k = ", k, "): ", mn, " (", cr[1], ", ", cr[2], ")")
}


results <- c()
schedules <- c("2 months","4 months", "6 months", "12 months")

for (sch in 1:schedules) {
  results <- summary_fun_C("Concomitant", 1, dat_vis)
  results
}

summary_fun_C("Concomitant", 1, dat_vis)
summary_fun_C("Separate", 1, dat_vis)
summary_fun_C("Separate", 2, dat_vis)

#####################################

hist(plogis(rnorm(10000, -2, 1)))
hist(plogis(rnorm(10000, -3, 1)))
summary(plogis(rnorm(10000, -3, 1)))
quantile(plogis(rnorm(10000, -3, 1)), probs = c(0.025, 0.975))
quantile(plogis(rnorm(10000, -3.5, 1)), probs = c(0.025, 0.975))
quantile(plogis(rnorm(10000, -4, 1)), probs = c(0.025, 0.975))
hist(plogis(rnorm(10000, -4, 1)))
hist(plogis(rnorm(10000, -1, 1)))
hist(plogis(rnorm(10000, -1, 0.5)))

