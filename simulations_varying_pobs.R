devtools::load_all()

stop_time <- 30
death_rate <- 0.1

# peaked <- function(t) {
#   stop_time <- 40
#   peak_time <- 20
#   min <- 0.2
#   max <- 0.2
#
#   if (t < peak_time) {
#     slope <- (max-min)/(peak_time)
#     return(min + slope*t)
#   } else {
#     slope <- (min-max)/(stop_time-peak_time)
#     return(max + slope*(t-peak_time))
#   }
# }

set.seed(1)
epi <- epi_con(birth_rate = 0.3, death_rate = death_rate, stop_time = stop_time, x0 = 1)
prev <- prevalence(epidemic = epi, stop_time = stop_time)

ptree <- phylo_tree(epi, stop_time)
set.seed(3)
sample05 <- sample_phylo(ptree = ptree, ptree_lag = 0, pi0 = 0, pi1 = 0.05)

par(mfrow=c(1,1))
plot(prev, type="l")

# contant rho=0.15
set.seed(2)
noisy_prev0 <- sample_prevalence(prev, 0.15)
lines(noisy_prev0, col="red")
pobs_true0 <- sum(noisy_prev0[-1,2])/sum(prev[-1,2])

# rho=0.05 to rho=0.25
set.seed(2)
noisy_prev1 <- prev
noisy_prev1[2:(1+stop_time/2),2] <- rbinom(stop_time/2, prev[2:(1+stop_time/2),2], 0.05)
noisy_prev1[(2+stop_time/2):(1+stop_time),2] <- rbinom(stop_time/2, prev[(2+stop_time/2):(1+stop_time),2], 0.25)
lines(noisy_prev1, col="blue")
pobs_true1 <- sum(noisy_prev1[-1,2])/sum(prev[-1,2])

# rho=0.25 to rho=0.05
set.seed(2)
noisy_prev2 <- prev
noisy_prev2[2:(1+stop_time/2),2] <- rbinom(stop_time/2, prev[2:(1+stop_time/2),2], 0.25)
noisy_prev2[(2+stop_time/2):(1+stop_time),2] <- rbinom(stop_time/2, prev[(2+stop_time/2):(1+stop_time),2], 0.05)
lines(noisy_prev2, col="darkgreen")
pobs_true2 <- sum(noisy_prev2[-1,2])/sum(prev[-1,2])

iter <- 10000
max_time <- 10*60
sigma0 <- 0.05
reporting_prob0 <- 0.15
x0 <- 1
# death_rate <- 0.1
n_particles <- 1000
print <- T

library(doParallel)
registerDoParallel(cores=9)

chains <- foreach(i=0:8) %dopar% {
  if (i %in% 0:2) {
    seed <- 4
  }
  if (i %in% 3:5) {
    seed <- 5
  }
  if (i %in% 6:8) {
    seed <- 6
  }
  set.seed(seed)
  pmmh(iter = iter,
       max_time = max_time,
       sigma0 = sigma0,
       reporting_prob0 = reporting_prob0,
       x0 = x0,
       death_rate = death_rate,
       ptree = sample05$ptree,
       ptree_lag = round(sample05$ptree_lag,0),
       sample_prevalence = eval(parse(text=paste0("noisy_prev",i %% 3))),
       n_particles = n_particles,
       print = print)
}

chain40 <- chains[[1]]
chain41 <- chains[[2]]
chain42 <- chains[[3]]
chain50 <- chains[[4]]
chain51 <- chains[[5]]
chain52 <- chains[[6]]
chain60 <- chains[[7]]
chain61 <- chains[[8]]
chain62 <- chains[[9]]

beepr::beep(3)

library(scales)
width <- 8
height <- width
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sim_varying_pobs_trace.pdf", height=height, width=width)
par(mfrow=c(3,3))
par(mar=c(4,4,2,1))
plot(chain50$reporting_prob, type="l", ylim=c(0,0.6), xlab="Iteration", ylab="Reporting probability", main="rho=0.15", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain50$reporting_prob, na.rm=T), col="red")
# abline(h=pobs_true0, col="blue")
plot(chain51$reporting_prob, type="l", ylim=c(0,0.6), xlab="Iteration", ylab="Reporting probability", main="rho=0.05 then rho=0.25", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain51$reporting_prob, na.rm=T), col="red")
# abline(h=pobs_true1, col="blue")
plot(chain52$reporting_prob, type="l", ylim=c(0,0.6), xlab="Iteration", ylab="Reporting probability", main="rho=0.25 then rho=0.05", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain52$reporting_prob, na.rm=T), col="red")
# abline(h=pobs_true2, col="blue")
plot(chain50$sigma, type="l", ylim=c(0,0.4), xlab="Iteration", ylab="Smoothness", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain50$sigma, na.rm=T), col="red")
plot(chain51$sigma, type="l", ylim=c(0,0.4), xlab="Iteration", ylab="Smoothness", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain51$sigma, na.rm=T), col="red")
plot(chain52$sigma, type="l", ylim=c(0,0.4), xlab="Iteration", ylab="Smoothness", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain52$sigma, na.rm=T), col="red")
plot(chain50$x0, type="l", ylim=c(0,20), xlab="Iteration", ylab="Day 0 prevalence", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain50$x0, na.rm=T), col="red")
plot(chain51$x0, type="l", ylim=c(0,20), xlab="Iteration", ylab="Day 0 prevalence", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain51$x0, na.rm=T), col="red")
plot(chain52$x0, type="l", ylim=c(0,20), xlab="Iteration", ylab="Day 0 prevalence", col=alpha(colour="black",alpha=0.75))
# abline(h=mean(chain52$x0, na.rm=T), col="red")
dev.off()

burn_in <- 5000
par(mfrow=c(1,3))
plot(density(chain50$reporting_prob[-(1:burn_in)]), xlim=c(0,0.5), ylim=c(0,25), xlab="Reporting probability", main="", col="grey10")
lines(density(chain51$reporting_prob[-(1:burn_in)]), col="grey40")
lines(density(chain52$reporting_prob[-(1:burn_in)]), col="grey70")
plot(density(chain50$sigma[-(1:burn_in)]), xlim=c(0,0.4), ylim=c(0,60), xlab="Smoothness", main="", col="grey10")
lines(density(chain51$sigma[-(1:burn_in)]), col="grey40")
lines(density(chain52$sigma[-(1:burn_in)]), col="grey70")
hist(chain50$x0[-(1:burn_in)], freq=F, xlim=c(0,20), ylim=c(0,0.4), xlab="Day 0 prevalence", col=rgb(0.1,0.1,0.1,0.5))
hist(chain51$x0[-(1:burn_in)], freq=F, col=rgb(0.4,0.4,0.4,0.5), add=T)
hist(chain52$x0[-(1:burn_in)], freq=F, col=rgb(0.7,0.7,0.7,0.5), add=T)


width <- 8
height <- width/3
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sim_varying_pobs_br.pdf", height=height, width=width)
br_mean0 <- apply(chain50$birth_rate[-(1:burn_in),], 2, mean, na.rm=T)
br_lcl0 <- matrixStats::colQuantiles(chain50$birth_rate[-(1:burn_in),], probs=0.025, na.rm=T)
br_ucl0 <- matrixStats::colQuantiles(chain50$birth_rate[-(1:burn_in),], probs=0.975, na.rm=T)
br_mean1 <- apply(chain51$birth_rate[-(1:burn_in),], 2, mean, na.rm=T)
br_lcl1 <- matrixStats::colQuantiles(chain51$birth_rate[-(1:burn_in),], probs=0.025, na.rm=T)
br_ucl1 <- matrixStats::colQuantiles(chain51$birth_rate[-(1:burn_in),], probs=0.975, na.rm=T)
br_mean2 <- apply(chain52$birth_rate[-(1:burn_in),], 2, mean, na.rm=T)
br_lcl2 <- matrixStats::colQuantiles(chain52$birth_rate[-(1:burn_in),], probs=0.025, na.rm=T)
br_ucl2 <- matrixStats::colQuantiles(chain52$birth_rate[-(1:burn_in),], probs=0.975, na.rm=T)
par(mfrow=c(1,3))
par(mar=c(4,4,2,1))
plot(br_mean0, type="l", xlab="Day", ylab="Birth rate", ylim=c(0,0.8), main="rho=0.15")
lines(br_lcl0, lty=2)
lines(br_ucl0, lty=2)
abline(h=0.3, col="red")
# lines(1:30, c(seq(0.1, 0.3, length.out=16)[-1],seq(0.3,0.1,length.out=16)[-1]), col="red")
plot(br_mean1, type="l", xlab="Day", ylab="Birth rate", ylim=c(0,0.8), main="rho=0.05 then rho=0.25")
lines(br_lcl1, lty=2)
lines(br_ucl1, lty=2)
abline(h=0.3, col="red")
# lines(1:30, c(seq(0.1, 0.3, length.out=16)[-1],seq(0.3,0.1,length.out=16)[-1]), col="red")
plot(br_mean2, type="l", xlab="Day", ylab="Birth rate", ylim=c(0,0.8), main="rho=0.25 then rho=0.05")
lines(br_lcl2, lty=2)
lines(br_ucl2, lty=2)
abline(h=0.3, col="red")
# lines(1:30, c(seq(0.1, 0.3, length.out=16)[-1],seq(0.3,0.1,length.out=16)[-1]), col="red")
dev.off()

which(br_lcl1 > 0.3)
which(br_ucl2 < 0.3)

xt_mean0 <- apply(chain50$prevalence[-(1:burn_in),], 2, mean, na.rm=T)
xt_lcl0 <- matrixStats::colQuantiles(chain50$prevalence[-(1:burn_in),], probs=0.025, na.rm=T)
xt_ucl0 <- matrixStats::colQuantiles(chain50$prevalence[-(1:burn_in),], probs=0.975, na.rm=T)
xt_mean1 <- apply(chain51$prevalence[-(1:burn_in),], 2, mean, na.rm=T)
xt_lcl1 <- matrixStats::colQuantiles(chain51$prevalence[-(1:burn_in),], probs=0.025, na.rm=T)
xt_ucl1 <- matrixStats::colQuantiles(chain51$prevalence[-(1:burn_in),], probs=0.975, na.rm=T)
xt_mean2 <- apply(chain52$prevalence[-(1:burn_in),], 2, mean, na.rm=T)
xt_lcl2 <- matrixStats::colQuantiles(chain52$prevalence[-(1:burn_in),], probs=0.025, na.rm=T)
xt_ucl2 <- matrixStats::colQuantiles(chain52$prevalence[-(1:burn_in),], probs=0.975, na.rm=T)
par(mfrow=c(1,3))
plot(0:stop_time,xt_mean0, type="l", ylab="Prevalence", ylim=c(0,2500), main="0.15")
lines(0:stop_time,xt_lcl0, lty=2)
lines(0:stop_time,xt_ucl0, lty=2)
lines(prev, col="red")
plot(0:stop_time,xt_mean1, type="l", ylab="Prevalence", ylim=c(0,8000), main="0.05 then 0.25")
lines(0:stop_time,xt_lcl1, lty=2)
lines(0:stop_time,xt_ucl1, lty=2)
lines(prev, col="red")
plot(0:stop_time,xt_mean2, type="l", ylab="Prevalence", ylim=c(0,2000), main="0.25 then 0.05")
lines(0:stop_time,xt_lcl2, lty=2)
lines(0:stop_time,xt_ucl2, lty=2)
lines(prev, col="red")

posteriors <- array(NA, dim=c(3,3,3), dimnames=list(c("15%", "5% then 25%", "25% then 5%"), c("pobs", "sigma", "x0"), c("mean", "lcl", "ucl")))
for (i in 1:3) {
  current_chain5 <- eval(parse(text=paste0("chain5",i-1)))
  for (j in 1:3) {
    current_variable <- ifelse(j==1, eval(parse(text="current_chain5$reporting_prob")), ifelse(j==2, eval(parse(text="current_chain5$sigma")), eval(parse(text="current_chain5$x0"))))
    posteriors[i,j,1] <- mean(current_variable[-(1:burn_in)])
    posteriors[i,j,2:3] <- quantile(current_variable[-(1:burn_in)], c(0.025, 0.975))
  }
}
posteriors
