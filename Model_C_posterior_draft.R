##############Model C analysis script - prob of reporting at least one AEFI/Impact following (1) concom (2) NIP first (3) Men B first
library(cmdstanr)
library(posterior)
library(ggplot2)


#load in data
#dat <- read.csv(file = "C:/Users/ETay/Documents/Work documents/AVS work/Thuy_concom/dat_modelC.csv")

N_C = nrow(dat)                    ## number of responders
N_strat_C  = 3                     ## number of strategies
N_sched_C  = 4                     ## number of schedules

s_C <- dat$vax_sequence            ## vaccine strategy, 1 = "Concomitant vaccination", 2 = "NIP first", 3 = "Men B first"
t_C <- dat$sched                   ## schedule - 1 = 2 months, 2 = 4 months, 3 = 6 months, 4 = 12 months
w_C <- dat$sex                     ## sex - 0 = "Male", 1 = "Female"
x_C <- dat$indig                   ## Indigenous status -0 = Non-indig, 1 = Aboriginal and Torres Strait Islander
z_C <- dat$PMH                     ## comorbidity

y <- dat$any_event        ## outcome - 0,1,2
y <- dat$impact           ## outcome - 0,1,2
y <- dat$ma               ## outcome - 0,1,2

a <- qlogis(0.3)           ## prior distribution mean for one event - depends on the schedule and the vaccine strategy?
b <-                       ## prior distribution standard deviation for two events
c <- qlogis(0.05)          ## prior distribution mean for two events
d <-                       ## prior distribution for standard deviation for two events

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



model_C <- cmdstan_model(write_stan_file(readLines("C:/Users/etay/Documents/Work documents/AVS work/Thuy_concom/model_C.stan")))###

fit_C <- model_C$sample(data = dat_C, 
                  chains = 8, parallel_chains = 8)###

postr <- posterior::as_draws_matrix(fit_C$draws())

##SAP code
model_C <- cmdstan_model("C:/Users/etay/Documents/Work documents/AVS work/Thuy_concom/model_C.stan")

drp <- capture.output({fit_C <- model_C$sample(dat_C, chains = 8, parallel_chains = 8)})
draws_C <- as_draws_matrix(fit_C$draws(c("mu")))
dat_vis_C <- data.frame(x = plogis(as.vector(draws_C)),
                        strategy = rep(rep(c("Concomitant", "MenB first", "NIP first"), each = 8000), 3),
                        Parameter = rep(c("P(k = 1)", "P(k = 2)"), each = 8000*3))
dat_vis_C <- dat_vis_C[!(dat_vis_C$strategy == "Concomitant" & dat_vis_C$Parameter == "P(k = 2)"),]

ggplot(dat_vis_C, aes(x = x, colour = Parameter)) +
  facet_wrap(~strategy) +
  geom_density() +
  scale_colour_manual(values = colours)

summary_fun_C <- function(strat, k, dat){
  mn <- format(round(mean(dat[dat$strategy == strat & dat$Parameter == paste0("P(k = ", k, ")"),]$x), 2), nsmall = 2)
  cr <- format(round(quantile(dat[dat$strategy == strat & dat$Parameter == paste0("P(k = ", k, ")"),]$x, c(0.025, 0.975)), 2), nsmall = 2)
  paste0(strat, " strategy P(k = ", k, "): ", mn, " (", cr[1], ", ", cr[2], ")")
}

summary_fun_C("Concomitant", 1, dat_vis_C)
summary_fun_C("NIP first", 1, dat_vis_C)
summary_fun_C("NIP first", 2, dat_vis_C)
summary_fun_C("MenB first", 1, dat_vis_C)
summary_fun_C("MenB first", 2, dat_vis_C)
