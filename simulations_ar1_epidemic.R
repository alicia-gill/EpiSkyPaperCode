devtools::load_all()

stop_time <- 30
death_rate <- 0.1

set.seed(1)
epi <- epi_con(birth_rate = 0.3, death_rate = death_rate, stop_time = stop_time, x0 = 1)
prev <- prevalence(epidemic = epi, stop_time = stop_time)

ptree <- phylo_tree(epi, stop_time)
set.seed(3)
sample05 <- sample_phylo(ptree = ptree, ptree_lag = 0, pi0 = 0, pi1 = 0.05)

# contant rho=0.15
set.seed(2)
noisy_prev0 <- sample_prevalence(prev, 0.15)
pobs_true0 <- sum(noisy_prev0[-1,2])/sum(prev[-1,2])

sum <- 0
for (i in 1:1000) {
  sum <- 0.1 + sum/3
}

p_x <- 0.1
p_y <- 1/3

# AR(1) + binomial mixture
noisy_prev <- prev
set.seed(2)
noisy_prev[2,2] <- rbinom(n=1, size=prev[2,2], p=p_x)
for (i in 3:nrow(noisy_prev)) {
  noisy_prev[i,2] <- rbinom(n=1, size=prev[i,2], p=p_x) + rbinom(n=1, size=noisy_prev[i-1,2], p=p_y)
}
plot(noisy_prev, type="l")
pobs_true <- sum(noisy_prev[-1,2])/sum(prev[-1,2])

width <- 6
height <- width*0.75
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sim_ar1_epi_prev.pdf", height=height, width=width)
par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
plot(noisy_prev0, type="l", ylim=c(0,250), xlab="Day", ylab="Observed prevalence")
lines(noisy_prev, col="red")
legend("topleft", legend=c("correctly specified","misspecified"), col=c("black", "red"), lty=1)
dev.off()

iter <- 10000
max_time <- 10*60
sigma0 <- 0.05
reporting_prob0 <- 0.1
x0 <- 1
# death_rate <- 0.1
n_particles <- 1000
print <- T

library(doParallel)
registerDoParallel(cores=4)

chains <- foreach(i=0:3) %dopar% {
  if (i != 3) {
    set.seed(i+4)
    pmmh(iter = iter,
         max_time = max_time,
         sigma0 = sigma0,
         reporting_prob0 = reporting_prob0,
         x0 = x0,
         death_rate = death_rate,
         ptree = sample05$ptree,
         ptree_lag = round(sample05$ptree_lag,0),
         sample_prevalence = noisy_prev,
         n_particles = n_particles,
         print = print)
  } else {
    set.seed(5)
    pmmh(iter = iter,
         max_time = max_time,
         sigma0 = sigma0,
         reporting_prob0 = reporting_prob0,
         x0 = x0,
         death_rate = death_rate,
         ptree = sample05$ptree,
         ptree_lag = round(sample05$ptree_lag,0),
         sample_prevalence = noisy_prev0,
         n_particles = n_particles,
         print = print)
  }
}

chain1 <- chains[[1]]
chain2 <- chains[[2]]
chain3 <- chains[[3]]
chain50 <- chains[[4]]

par(mfrow=c(3,3))
par(mar=c(4,4,3,1))
plot(chain1$reporting_prob, type="l", xlab="iteration", ylab="reporting prob", main="seed 4")
plot(chain2$reporting_prob, type="l", xlab="iteration", ylab="reporting prob", main="seed 5")
plot(chain3$reporting_prob, type="l", xlab="iteration", ylab="reporting prob", main="seed 6")
plot(chain1$sigma, type="l", xlab="iteration", ylab="smoothness")
plot(chain2$sigma, type="l", xlab="iteration", ylab="smoothness")
plot(chain3$sigma, type="l", xlab="iteration", ylab="smoothness")
plot(chain1$x0, type="l", xlab="iteration", ylab="day 0 prev")
plot(chain2$x0, type="l", xlab="iteration", ylab="day 0 prev")
plot(chain3$x0, type="l", xlab="iteration", ylab="day 0 prev")

burn_in <- 5000
br_mean0 <- apply(chain50$birth_rate[-(1:burn_in),],2,mean)
br_lcl0 <- matrixStats::colQuantiles(chain50$birth_rate[-(1:burn_in),], probs=0.025)
br_ucl0 <- matrixStats::colQuantiles(chain50$birth_rate[-(1:burn_in),], probs=0.975)
bt_mean <- apply(chain1$birth_rate[-(1:burn_in),],2,mean)
bt_lcl <- matrixStats::colQuantiles(chain1$birth_rate[-(1:burn_in),], probs=0.025)
bt_ucl <- matrixStats::colQuantiles(chain1$birth_rate[-(1:burn_in),], probs=0.975)

width <- 8
height <- width/3
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sim_ar1_epi_post.pdf", height=height, width=width)
par(mfrow=c(1,3))
par(mar=c(4,4,1,1))
plot(density(chain50$reporting_prob[-(1:burn_in)]), ylim=c(0,15), xlab="Reporting probability", main="", col="grey30")
lines(density(chain1$reporting_prob[-(1:burn_in)]), col="grey70")
plot(density(chain50$sigma[-(1:burn_in)]), ylim=c(0,60), xlab="Smoothness", main="", col="grey30")
lines(density(chain1$sigma[-(1:burn_in)]), col="grey70")
hist(chain50$x0[-(1:burn_in)], ylim=c(0,0.4), xlab="Day 0 prevalence", main="", freq=F, xlim=c(0,15), col="grey30")
hist(chain1$x0[-(1:burn_in)], col=rgb(0.7,0.7,0.7,0.5), add=T, freq=F)
dev.off()

width <- 8
height <- width/2
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sim_ar1_epi_birth.pdf", height=height, width=width)
par(mfrow=c(1,2))
par(mar=c(4,4,2,1))
plot(br_mean0, type="l", xlab="Day", ylab="Birth rate", main="Correctly specified", ylim=c(0.1, 0.5))
lines(br_lcl0, lty=2)
lines(br_ucl0, lty=2)
abline(h=0.3, col="red")
plot(bt_mean, type="l", xlab="Day", ylab="Birth rate", main="Misspecified", ylim=c(0.1, 0.5))
lines(bt_lcl, lty=2)
lines(bt_ucl, lty=2)
abline(h=0.3, col="red")
dev.off()
