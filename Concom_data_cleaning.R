
##########################################

#Model C data sim from SAP
N_C <- 1000
N_strat_C <- 2
N_sched_C <- 1
s_C <- sample(1:N_strat_C, size = N_C, replace = TRUE)
t_C <- rep(1:N_sched_C, N_C/N_sched_C)
w_C <- sample(c(0,1), size = N_C, replace = TRUE)
x_C <- sample(c(0,1), size = N_C, replace = TRUE)
z_C <- sample(c(0,1), size = N_C, replace = TRUE)
mu_C <- array(NA, dim = c(N_strat_C, N_sched_C, 2))
mu_C[,,1] <- qlogis(0.3)
mu_C[2,,2] <- qlogis(0.05)
beta1_C <- beta2_C <- 0
gamma1_C <- gamma2_C <- 0
delta1_C <- delta2_C <- 0
p_C <- eta_C <- matrix(0, nrow = N_C, ncol = 3)
for(i in 1:N_C){
  eta_C[i, 2] <- mu_C[s_C[i], t_C[i], 1] +
    w_C[i]*beta1_C +
    x_C[i]*gamma1_C +
    z_C[i]*delta1_C
  if(s_C[i] == 1){
    eta_C[i, 3] <- -Inf
  } else {
    eta_C[i, 3] <- mu_C[s_C[i], t_C[i], 2] +
      w_C[i]*beta2_C +
      x_C[i]*gamma2_C +
      z_C[i]*delta2_C
  }
  p_C[i,] <- exp(eta_C[i,])/sum(exp(eta_C[i,]))
}
y_C <- sapply(1:N_C, function(i) sample(1:3, size = 1, prob = p_C[i,]))
a <- qlogis(0.3)
b <- d <- 2
c <- qlogis(0.05)
dat_C <- list(N = N_C,
              N_strat = N_strat_C,
              N_sched = N_sched_C,
              N_clinic = N_clinic,
              s = s_C, #strategy/vax seqeunce - int
              t = t_C, #schedule - int
              w = w_C, #sex
              x = x_C, #indig
              z = z_C, #comorbidity
              y = y_C, #outcome - int
              a = a,
              b = b,
              c = c,
              d = d)

str(dat_C)
##########################################################################################################
##Model C prob of reporting at least one AEFI/IMPACT  following (1)Concom, (2) NIP first or (3) Men B first

library(tidyverse)
library(stringr)
library(data.table)
library(here)

load("Z:/Analyses/Concomitant vaccination/Infant_NIP_MenB_concom_modelling_MA added.rda")#connect to USyd RDS)
str(infant)

dat <- as.data.table(infant)
colnames(dat) <- tolower(colnames(dat))
#11101
dat <- dat[!(is.na(any_event)) &
             !(is.na(schedule)) &
             !(is.na(atsi)) &
             !(is.na(sex)) & 
             !(is.na(group)) &
             !(is.na(clinic_state)) &
             !(is.na(clinic_type)) &
             #!is.na(fever)) & 4 records missing a value for fever
             !(is.na(pmh)), ]#11088
dat$vax_time_diff[is.na(dat$vax_time_diff)] <- 0

table(dat$clinic_state)
# ACT  NSW   NT  QLD   SA  TAS  VIC   WA 
# 702 3212   24 3163 1785  219 1000  996

table(dat$clinic_type)
# Aboriginal Health Service          General Practice              State Health 
#              214                     10419                       468 

dat <- dat[dat$vax_time_diff != "367",] #11087

# dat$group <- str_replace_all(dat$group, c("Concomitant vaccination" = "1", "Seperate vaccination" = "2"))
# dat$group <- as.integer(dat$group)
# 
# table(dat$group,dat$schedule)

#     1    2    3    4
# 1 2118 3113 1156 2351
# 2  653  715  250  731

dat$vax_sequence[is.na(dat$vax_sequence)] <- "Concomitant vaccination"
dat$vax_sequence <- str_replace_all(dat$vax_sequence, c("Concomitant vaccination" = "1", "NIP first" = "2", "MenB first" = "3"))
dat$vax_sequence <- as.integer(dat$vax_sequence)
table(dat$vax_sequence)

#1    2    3 
#8738 2208  141 

#exclude Men B first
dat <- dat %>%
  filter(!(vax_sequence=="3")) #10946
#s <- dat$vax_sequence 

dat$any_event <- as.integer(dat$any_event)
table(dat$any_event)
#  0    1    2 
# 6305 4301  340 

dat$impact <- as.integer(dat$impact)
table(dat$impact)
#   0     1     2 
# 10496   444     6 

dat$ma <- as.integer(dat$ma)
table(dat$ma)

# 0     1     2 
# 10541   401     4 

dat <- dat %>%
  mutate(ma = case_when(ma == "2" ~ 1, .default = ma))

table(dat$ma)

# 0     1 
# 10541   405 

dat$local <- as.integer(dat$local)
dat$local[is.na(dat$local)] <- 0
table(dat$local)

#  0    1    2 
# 7773 3000  173 

dat$fever <- as.integer(dat$fever)
dat$fever[is.na(dat$fever)] <- 0
table(dat$fever)

# 0    1    2 
# 7921 2910  115 

dat$schedule <- str_replace_all(dat$schedule, c("2 months" = "1", "4 months" = "2",
                                                "6 months" = "3","11" = "4"))
dat$schedule <- as.integer(dat$schedule)
table(dat$vax_sequence, dat$schedule)


#    1    2    3    4
# 1 2118 3113 1156 2351
# 2  653  630  242  683

dat$sex <- str_replace_all(dat$sex, c("Female" = "1", "Male" = "0"))
dat$sex <- as.numeric(dat$sex)
setnames(dat, "atsi", "indig")

##clinic state 
# ACT  NSW   NT  QLD   SA  TAS  VIC   WA 
# 688 3162   24 3129 1783  216  977  967 

##clinic type 
# Aboriginal Health Service          General Practice              State Health 
# 214                     10265                       467 


###########################################################################
write.csv(dat, file = "C:/Users/ETay/Documents/NIP_MenB_dat.csv", row.names = FALSE)
###########################################################################
#looking for dupes in the data
table(table(infant$UID_PERSON))

1    2    3    4 
6602 1754  329    1 
# include these numbers - one line in the results who contributed data to 1,2,3,4 schedule points..
# 1 person responded 4 times...so ~ 8500 unique people

table(table(infant[infant$SCHEDULE == "2 months",]$UID_PERSON))


##################################
unique_values <- infant %>%
     group_by(UID_PERSON) %>%
     filter(n() == 1) 

duplicate_values <- infant %>%
  group_by(UID_PERSON) %>%
  filter(n() > 1) #4499

duplicate_values_2 <- duplicate_values %>%
  group_by(UID_PERSON) %>%
  filter(n() < 3) #3508

duplicate_values_3or4 <- duplicate_values %>%
  group_by(UID_PERSON) %>%
  filter(n() > 2) #991

duplicate_values_3 <- duplicate_values_3or4[duplicate_values_3or4$UID_PERSON != "HHS-24064",]

subset_df <- duplicate_values_3[duplicate_values_3$ANY_EVENT == 1, ] #425
subset_df_3 <- subset_df %>%
  group_by(UID_PERSON) %>%
  filter(n() > 2) #156

subset_df_1 <- subset_df %>%
  group_by(UID_PERSON) %>%
  filter(n() <2) #95

subset_df_2 <- subset_df %>%
  group_by(UID_PERSON) %>%
  filter(n() == 2)

duplicate_values_4 <- infant %>%
      group_by(UID_PERSON) %>%
     filter(n() > 3) #4
#####################
# q[1:10,]
# [,1] [,2] [,3] [,4] [,5] [,6] [,7]
# [1,]    0    0    1    0    0    0    0
# [2,]    0    0    1    0    0    0    0
# [3,]    0    0    1    0    0    0    0
# [4,]    0    0    1    0    0    0    0
# [5,]    0    0    1    0    0    0    0
# [6,]    0    0    1    0    0    0    0
# [7,]    0    0    1    0    0    0    0
# [8,]    0    0    1    0    0    0    0
# [9,]    0    0    1    0    0    0    0
# [10,]    0    0    1    0    0    0    0
# colSums(q)
# [1]  688   24 3129 1783  216  977  967
# table(rowSums(q))
# 
# 0    1 
# 3162 7784 