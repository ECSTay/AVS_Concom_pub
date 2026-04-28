##############Model C analysis script - prob of reporting at least one AEFI/Impact following (1) concom (2) NIP first (3) Men B first
library(cmdstanr)
library(posterior)
library(ggplot2)
library(forcats)
library(data.table)
library(bayesplot)

#load in data
dat <- fread(file = "C:/Users/ETay/Documents/NIP_MenB_dat.csv")
dat$res <- as.integer(factor(dat$uid_person, levels = unique(dat$uid_person)))
dat$vax_sequence[dat$vax_sequence == "2"] <- 0
dat$vax_sequence <- as.integer(dat$vax_sequence)
N_R <- dat[,.N]                                            ## number of responses
N_I <- dat[,length(unique(dat$uid_person))]                ## number of unique respondents

N_strat  = 2                    ## number of strategies
N_sched  = 4                    ## number of schedules
N_state  = 7                    
N_clinic = 2

infant <- dat$res

s <- dat$vax_sequence           ## vaccine strategy, 1 = "Concomitant vaccination" --->1, 2 = "Separate" --->0
t <- dat$sched                  ## schedule - 1 = 2 months, 2 = 4 months, 3 = 6 months, 4 = 12 months
w <- dat$sex                    ## sex - 0 = "Male", 1 = "Female"
x <- dat$indig                  ## Indigenous status - 0 = Non-Indig, 1 = Aboriginal and Torres Strait Islander

q <- cbind(as.integer(dat$clinic_state == "ACT"),
as.integer(dat$clinic_state == "NT"),
as.integer(dat$clinic_state == "QLD"),
as.integer(dat$clinic_state == "SA"),
as.integer(dat$clinic_state == "TAS"),
as.integer(dat$clinic_state == "VIC"),
as.integer(dat$clinic_state == "WA"))

c <- cbind(as.integer(dat$clinic_type == "Aboriginal Health Service"),
           as.integer(dat$clinic_type == "State Health"))

z <- dat$pmh                    ## comorbidity - 0 - None, 1 - at least one

y <- dat$any_event + 1           ## outcome - 0,1,2 p(1) = 0.49, p(2) = 0.03
#y_C <- dat$impact + 1              ## outcome - 0,1,2 p(1) = 0.04, p(2) = 0.0005
#y_C <- dat$ma + 1                  ## outcome - 0,1,2 p(1) = 0.02
#y_C <- dat$local + 1               ## outcome - 0,1,2 p(1) = 0.27, p(2) = 0.02
#y_C <- dat$fever + 1               ## outcome - 0,1,2 p(1) = 0.27, p(2) = 0.01

#a = -3, b = 1  for P(impact), P(MA)
#a = -1, b = 1  for P(AEFI), P(local), P(fever)

# a <- -3                       ## prior distribution mean for one event - depends on the schedule and the vaccine strategy
# b <- 1                         ## prior distribution standard deviation for two events
# c <- -1                        ## prior distribution mean for two events
# d <- 1                         ## prior distribution for standard deviation for two events

dat_C <- list(N_R = N_R,
              N_I = N_I,
              N_strat = N_strat,
              N_sched = N_sched,
              N_state = N_state,
              N_clinic = N_clinic,
              infant = infant,
              s = s,
              t = t,
              w = w,
              x = x,
              q = q,
              c = c,
              z = z,
              y = y)
                # a = a,
                # b = b,
                # c = c,
                # d = d)


##SAP code
model_C <- cmdstan_model("C:/Users/ETay/Documents/Work documents/NIP_MenB revisions/Concom_new_model.stan")

time_start <- Sys.time()
#fit_C <- model_C$sample(dat_C, chains = 2, parallel_chains = 2) #43.3 minutes
fit_C <- model_C$sample(dat_C, chains = 8, parallel_chains = 8) 
time_finish <- Sys.time() - time_start

postr <-fit_C$draws()

saveRDS(fit_C, file = "fit_concom_AEFI_test.rds")
saveRDS(postr, file = "postr_concom_AEFI_test.rds")


saveRDS(postr, "postr_concom_AEFI.rds") # any AE
#saveRDS(postr, "postr_concom_impact.rds") # any days of impact
#saveRDS(postr, "postr_concom_ma.rds") # MA
#saveRDS(postr, "postr_concom_local.rds") # local reaction
#saveRDS(postr, "postr_concom_fever.rds") # fever
###########################################################################################
#load in the relevant postr
## marginalise - using option b) The event probability (e.g., one/two MAs) for an average person from the population of survey responders receiving their X month schedule.

postr <- posterior::as_draws_matrix(postr)
postr <- as_draws_matrix(fit_C$draws(c("mu", "beta", "gamma", "rho", "tau","delta")))

draws_marg <- array(NA, dim = c(nrow(postr), 2, 4, 2), 
                    dimnames = list(draw = 1:nrow(postr), 
                                    strategy = c("separate", "concom"), 
                                    schedule = c("2 months", "4 months", "6 months", "12 months"), 
                                    parameter = c("one", "two")))
sex_weights <- dat[, .(mn = mean(sex)), by = schedule]
indig_weights <- dat[, .(mn = mean(indig)), by = schedule]
comorb_weights <- dat[, .(mn = mean(pmh)), by = schedule]

for(strat in 1:2){
  for(sched in 1:4){
    for(par in 1:2){
      sex_par <- postr[,paste0("beta[", par, "]")]*sex_weights[schedule == sched]$mn
      indig_par <- postr[,paste0("gamma[", par, "]")]*indig_weights[schedule == sched]$mn
      comorb_par <- postr[,paste0("delta[", par, "]")]*comorb_weights[schedule == sched]$mn
      mu_par <- postr[,paste0("mu[", strat, ",", sched, ",", par, "]")]
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

#local reaction

ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(strategy~forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k >= 1)" = "#009E73"),
                      labels = c("One report of local reaction", "At least one report of local reaction"))


ggsave("Local_3.png", dpi = 400, width = 6, height = 5, units = "in")
#fever
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(strategy~forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k >= 1)" = "#009E73"),
                      labels = c("One report of fever", "At least one report of fever"))

ggsave("Fever_3.png", dpi = 400, width = 6, height = 5, units = "in")
#MA
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(strategy~forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k >= 1)" = "#009E73"),
                      labels = c("One report of MA", "At least one report of MA"))

ggsave("MA_3.png", dpi = 400, width = 6, height = 5, units = "in")
#impact
ggplot(dat_vis, aes(x = x, colour = Parameter)) +
  facet_grid(strategy~forcats::fct_relevel(schedule, "2 months", "4 months", "6 months", "12 months")) +
  geom_density() +
  xlab("Probability") +
  ylab("Density") +
  labs(color = NULL) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("P(k = 1)" = "#D55E00", "P(k >= 1)" = "#009E73"),
                      labels = c("One report of carer impact", "At least one report of carer impact"))

ggsave("Impact_3.png", dpi = 400, width = 6, height = 5, units = "in")

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

twomonths <- dat_vis[strategy == "Separate" & schedule == "2 months"]
one <- twomonths[Parameter == "P(k = 1)"]$x
two <- twomonths[Parameter == "P(k = 2)"]$x
atleastone <- one + two

fourmonths <- dat_vis[strategy == "Separate" & schedule == "4 months"]
one_4 <- fourmonths[Parameter == "P(k = 1)"]$x
two_4 <- fourmonths[Parameter == "P(k = 2)"]$x
atleastone_4 <- one_4 + two_4

sixmonths <- dat_vis[strategy == "Separate" & schedule == "6 months"]
one_6 <- sixmonths[Parameter == "P(k = 1)"]$x
two_6 <- sixmonths[Parameter == "P(k = 2)"]$x
atleastone_6 <- one_6 + two_6

twelvemonths <- dat_vis[strategy == "Separate" & schedule == "12 months"]
one_12 <- twelvemonths[Parameter == "P(k = 1)"]$x
two_12 <- twelvemonths[Parameter == "P(k = 2)"]$x
atleastone_12 <- one_12 + two_12

atleastone <- c(atleastone, atleastone_4, atleastone_6, atleastone_12)
strategy <- rep("Separate",32000)
schedule <- c(rep("2 months",8000), rep("4 months", 8000), rep("6 months", 8000), rep ("12 months", 8000))
Parameter <- rep("P(>= 1)")
dat_vis_sub <- data.frame(x = atleastone, )

######Diagnostic - to be included in the code to the thinkstation


postr <- fit_C$draws()
postr <- posterior::as_draws_matrix(fit_C$draws())
pars <- colnames(postr)[str_detect(colnames(postr), "mu|alpha|beta|delta|rho|tau|gamma|sigma")]

bayesplot::mcmc_neff_hist(neff_ratio(fit, pars = pars)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # plot area
    plot.background  = element_rect(fill = "white", color = NA)       # outer background
  )
ggplot2::ggsave("Arthralgia_flu_6m_neff_ratio.png", dpi = 400, width = 8, height = 9, units = "in")
bayesplot::mcmc_rhat_hist(rhat(fit, pars = pars)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # plot area
    plot.background  = element_rect(fill = "white", color = NA)       # outer background
  )
ggplot2::ggsave("Arthralgia_flu_6m_rhat.png", dpi = 400, width = 8, height = 9, units = "in")
bayesplot::mcmc_trace(postr, pars = sample(pars, size = 20)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # plot area
    plot.background  = element_rect(fill = "white", color = NA)       # outer background
  )
ggplot2::ggsave("Arthralgia_flu_6m_trace.png", dpi = 400, width = 8, height = 9, units = "in")

postr <- postr[,str_detect(colnames(postr), "mu|alpha|beta|delta|rho|tau|gamma|sigma")]

saveRDS(postr, file = "postr_6m_flu_arthralgia.rds")



