###Concom - ORs/effects

##Fixed effects for categorical variables 

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

############################################################

fnames <- c("AEFI", "MA", "impact", "fever", "local")
tab <- lapply(fnames, function(fn){
  postr <- readRDS(file = paste0("C:/Users/ETay/Documents/postr_concom", fn, ".rds"))
  postr <- posterior::as_draws_matrix(postr)
  postr <- postr[,str_detect(colnames(postr), "mu|alpha|beta|delta|rho|tau|gamma|sigma")]
  postr_fixed <- as.data.frame(postr[,str_detect(colnames(postr), "mu|alpha|beta|delta|rho|tau|gamma|sigma")])
  postr_fixed_summary <- sapply(exp(postr_fixed), function(par) 
    paste0(format(round(median(par), 1), nsmall = 1), " (", format(round(quantile(par, 0.025), 1), nsmall = 1), ", ", format(round(quantile(par, 0.975), 1), nsmall = 1), ")"))
  tmp <- data.table(parameter = names(postr_fixed_summary), postr_fixed_summary)
  names(tmp)[2] <- paste0("aOR_", fn)
  tmp
})
tab <- setmergelist(tab, on = "parameter")
write.csv(tab, "CORs.csv")