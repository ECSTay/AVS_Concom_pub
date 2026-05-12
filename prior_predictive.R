
# prior predictive distributions

rm(list = ls())

library(data.table)
library(ggplot2)

set.seed(5426)

## define fixed variables

N_sim <- 10000
N_R <- 1000
N_sched <- 4
t <- rep(1:N_sched, N_R/N_sched)
infant <- c(1:600, rep(601:750, each = 2), rep(751:770, each = 3), rep(771:780, each = 4))
N_I <- length(unique(infant))
s <- sample(0:1, size = N_R, replace = TRUE, prob = c(0.2, 0.8))
w <- sample(0:1, size = N_R, replace = TRUE, prob = c(0.5, 0.5))
x <- sample(0:1, size = N_R, replace = TRUE, prob = c(0.85, 0.15))
q_tmp <- sample(1:8, size = N_R, replace = TRUE, prob = c(0.05, 0.25, 0.23, 0.1, 0.2, 0.05, 0.02, 0.1)) ## include NSW
q <- matrix(0, nrow = N_R, ncol = 8)
q[cbind(1:N_R, q_tmp)] <- 1
q <- q[,-2]
c_tmp <- sample(1:3, size = N_R, replace = TRUE, prob = c(0.05, 0.8, 0.15)) ## include GP
c <- matrix(0, nrow = N_R, ncol = 3)
c[cbind(1:N_R, c_tmp)] <- 1
c <- c[,-2]
z <- sample(0:1, size = N_R, replace = TRUE, prob = c(0.95, 0.05))

## simulate parameters from their prior distributions

mu <- array(rnorm(N_sim*N_sched*2, 0, 2), dim = c(N_sim, N_sched, 2))
alpha_star <- rnorm(N_sim)
sigma_alpha <- 1/rgamma(N_sim, 3, 1)
alpha <- sapply(1:N_sched, function(t) rnorm(N_sim, alpha_star, sigma_alpha))
sigma_epsilon <- array(1/rgamma(N_sim*2, 3, 1), dim = c(N_sim, 2))
epsilon <- aperm(array(unlist(lapply(1:N_I, function(i)
                                 sapply(1:2, function(out)
                                    rnorm(N_sim, 0, sigma_epsilon[,out])))), dim = c(N_sim, 2, N_I)), c(1, 3, 2))
beta <- array(rnorm(N_sim*2), dim = c(N_sim, 2))
gamma <- array(rnorm(N_sim*2), dim = c(N_sim, 2))
rho <- array(rnorm(N_sim*7*2), dim = c(N_sim, 7, 2))
tau <- array(rnorm(N_sim*2*2), dim = c(N_sim, 2, 2))
delta <- array(rnorm(N_sim*2), dim = c(N_sim, 2))

## construct linear predictors for each strategy, schedule and outcome, and transform to probabilities

p <- array(NA, dim = c(N_sim, N_R, 3))
for(r in 1:N_R){
  eta <- array(NA, dim = c(N_sim, 3))
  eta[,1] <- 0
  eta[,2] <- mu[,t[r],1] + s[r]*alpha[,t[r]] + epsilon[,infant[r],1] + w[r]*beta[,1] + x[r]*gamma[,1] + as.vector(q[r,]%*%t(rho[,,1])) + as.vector(c[r,]%*%t(tau[,,1])) + z[r]*delta[,1]
  if(s[r] == 0){
    eta[,3] <- mu[,t[r],2] + epsilon[,infant[r],2] + w[r]*beta[,2] + x[r]*gamma[,2] + as.vector(q[r,]%*%t(rho[,,2])) + as.vector(c[r,]%*%t(tau[,,2])) + z[r]*delta[,2]
  } else {
    eta[,3] <- - Inf
  }
  p[,r,] <- exp(eta - apply(eta, 1, function(x) max(x) + log(sum(exp(x - max(x))))))
}

## simulate data

y <- array(NA_integer_, dim = c(N_sim, N_R))
for(r in 1:N_R) y[,r] <- apply(p[,r,], 1, function(p_sim) sample.int(3, size = 1, prob = p_sim) - 1)

## extract by strategy, schedule and outcome, and visualise

sim_df <- data.table(sim = rep(1:N_sim, N_R),
                     y   = as.vector(y),
                     s   = rep(s, each = N_sim),
                     t   = rep(t, each = N_sim))
sim_summ <- sim_df[, .(`0` = mean(y == 0), `1` = mean(y == 1), `2` = mean(y == 2)), by = .(sim, s, t)]
sim_summ[s == 1, `2` := NA]
sim_summ <- melt.data.table(sim_summ, id.vars = c("sim", "s", "t"), measure.vars = c("0", "1", "2"), variable.name = "Outcome")
sim_summ[, `:=`(
  s = factor(s, labels = c("Separate", "Concomitant")),
  t = factor(t, labels = c("2 months", "4 months", "6 months", "12 months"))
)]

ggplot(sim_summ, aes(x = value, colour = Outcome)) +
  scale_color_discrete(labels = c("No events","One event", "Two events")) +
  facet_grid(t ~ s) +
  xlab("Probability")+
  ylab("Density") +
  geom_density()

ggsave("Prior_NIPMenB.png", dpi = 400, width = 8, height = 8, units = "in")
