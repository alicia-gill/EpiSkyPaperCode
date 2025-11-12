files <- list.files(path="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/RData/", pattern="(hiv)")
for (i in 1:length(files)) {
  print(i)
  load(paste0("~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/RData/", files[i]))
}
gc()
devtools::load_all()

accept_g0e1 <- c(hiv_g0e1_1$acceptance_rate, hiv_g0e1_2$acceptance_rate, hiv_g0e1_3$acceptance_rate, hiv_g0e1_4$acceptance_rate, hiv_g0e1_5$acceptance_rate)*100
accept_g1e0 <- c(hiv_g1e0_1$acceptance_rate, hiv_g1e0_2$acceptance_rate, hiv_g1e0_3$acceptance_rate, hiv_g1e0_4$acceptance_rate, hiv_g1e0_5$acceptance_rate)*100
accept_g1e1_short <- c(hiv_g1e1_short_1$acceptance_rate, hiv_g1e1_short_2$acceptance_rate, hiv_g1e1_short_3$acceptance_rate, hiv_g1e1_short_4$acceptance_rate, hiv_g1e1_short_5$acceptance_rate)*100
accept_g1e1_long <- c(hiv_g1e1_long_1$acceptance_rate, hiv_g1e1_long_2$acceptance_rate, hiv_g1e1_long_3$acceptance_rate, hiv_g1e1_long_4$acceptance_rate, hiv_g1e1_long_5$acceptance_rate)*100

burn_in <- 10000
smcvar_g0e1 <- c(var(hiv_g0e1_1$smc_llik[-(1:burn_in)]),
                 var(hiv_g0e1_2$smc_llik[-(1:burn_in)]),
                 var(hiv_g0e1_3$smc_llik[-(1:burn_in)]),
                 var(hiv_g0e1_4$smc_llik[-(1:burn_in)]),
                 var(hiv_g0e1_5$smc_llik[-(1:burn_in)]))
smcvar_g1e0 <- c(var(hiv_g1e0_1$smc_llik[-(1:burn_in)]),
                 var(hiv_g1e0_2$smc_llik[-(1:burn_in)]),
                 var(hiv_g1e0_3$smc_llik[-(1:burn_in)]),
                 var(hiv_g1e0_4$smc_llik[-(1:burn_in)]),
                 var(hiv_g1e0_5$smc_llik[-(1:burn_in)]))
smcvar_g1e1_long <- c(var(hiv_g1e1_long_1$smc_llik[-(1:burn_in)],na.rm=T),
                      var(hiv_g1e1_long_2$smc_llik[-(1:burn_in)],na.rm=T),
                      var(hiv_g1e1_long_3$smc_llik[-(1:burn_in)],na.rm=T),
                      var(hiv_g1e1_long_4$smc_llik[-(1:burn_in)],na.rm=T),
                      var(hiv_g1e1_long_5$smc_llik[-(1:burn_in)],na.rm=T))

par(mfrow=c(5,3))
par(mar=c(4,4,1,1))
plot(hiv_g0e1_1$sigma, type="l")
plot(hiv_g0e1_1$reporting_prob, type="l")
plot(hiv_g0e1_1$x0, type="l")
plot(hiv_g0e1_2$sigma, type="l")
plot(hiv_g0e1_2$reporting_prob, type="l")
plot(hiv_g0e1_3$x0, type="l")
plot(hiv_g0e1_3$sigma, type="l")
plot(hiv_g0e1_3$reporting_prob, type="l")
plot(hiv_g0e1_3$x0, type="l")
plot(hiv_g0e1_4$sigma, type="l")
plot(hiv_g0e1_4$reporting_prob, type="l")
plot(hiv_g0e1_4$x0, type="l")
plot(hiv_g0e1_5$sigma, type="l")
plot(hiv_g0e1_5$reporting_prob, type="l")
plot(hiv_g0e1_5$x0, type="l")
accept_g0e1
smcvar_g0e1

par(mfrow=c(5,2))
par(mar=c(4,4,1,1))
plot(hiv_g1e0_1$sigma, type="l")
plot(hiv_g1e0_1$x0, type="l")
plot(hiv_g1e0_2$sigma, type="l")
plot(hiv_g1e0_3$x0, type="l")
plot(hiv_g1e0_3$sigma, type="l")
plot(hiv_g1e0_3$x0, type="l")
plot(hiv_g1e0_4$sigma, type="l")
plot(hiv_g1e0_4$x0, type="l")
plot(hiv_g1e0_5$sigma, type="l")
plot(hiv_g1e0_5$x0, type="l")
accept_g1e0
smcvar_g1e0

par(mfrow=c(5,3))
par(mar=c(4,4,1,1))
plot(hiv_g1e1_long_1$sigma, type="l")
plot(hiv_g1e1_long_1$reporting_prob, type="l")
plot(hiv_g1e1_long_1$x0, type="l")
plot(hiv_g1e1_long_2$sigma, type="l")
plot(hiv_g1e1_long_2$reporting_prob, type="l")
plot(hiv_g1e1_long_3$x0, type="l")
plot(hiv_g1e1_long_3$sigma, type="l")
plot(hiv_g1e1_long_3$reporting_prob, type="l")
plot(hiv_g1e1_long_3$x0, type="l")
plot(hiv_g1e1_long_4$sigma, type="l")
plot(hiv_g1e1_long_4$reporting_prob, type="l")
plot(hiv_g1e1_long_4$x0, type="l")
plot(hiv_g1e1_long_5$sigma, type="l")
plot(hiv_g1e1_long_5$reporting_prob, type="l")
plot(hiv_g1e1_long_5$x0, type="l")
accept_g1e1_long
smcvar_g1e1_long

hiv_g0e1 <- hiv_g0e1_2
hiv_g1e0 <- hiv_g1e0_3

burnin_g1e1_long <- 6777
hiv_g1e1_long <- list()
for (i in c(1,3,4,5)) {
  current <- eval(parse(text=paste0("hiv_g1e1_long_",i)))
  iter <- which.max(is.na(current$sigma))-1

  hiv_g1e1_long$birth_rate <- rbind(hiv_g1e1_long$birth_rate, current$birth_rate[burnin_g1e1_long:iter,])
  hiv_g1e1_long$prevalence <- rbind(hiv_g1e1_long$prevalence, current$prevalence[burnin_g1e1_long:iter,])
  hiv_g1e1_long$sigma <- c(hiv_g1e1_long$sigma, current$sigma[burnin_g1e1_long:iter])
  hiv_g1e1_long$reporting_prob <- c(hiv_g1e1_long$reporting_prob, current$reporting_prob[burnin_g1e1_long:iter])
  hiv_g1e1_long$x0 <- c(hiv_g1e1_long$x0, current$x0[burnin_g1e1_long:iter])
  hiv_g1e1_long$accept <- c(hiv_g1e1_long$accept, current$accept[burnin_g1e1_long:iter])
  hiv_g1e1_long$run_time <- c(hiv_g1e1_long$run_time, current$run_time)
  hiv_g1e1_long$smc_llik <- c(hiv_g1e1_long$smc_llik, current$smc_llik[burnin_g1e1_long:iter])
  hiv_g1e1_long$chain_no <- c(hiv_g1e1_long$chain_no, rep(i, iter-burnin_g1e1_long+1))
  hiv_g1e1_long$n_particles <- c(hiv_g1e1_long$n_particles, current$n_particles)
}

burnin_g1e1_long_nobi <- 0
hiv_g1e1_long_nobi <- list()
for (i in c(1,3,4,5)) {
  current <- eval(parse(text=paste0("hiv_g1e1_long_",i)))
  iter <- which.max(is.na(current$sigma))-1

  hiv_g1e1_long_nobi$birth_rate <- rbind(hiv_g1e1_long_nobi$birth_rate, current$birth_rate[burnin_g1e1_long_nobi:iter,])
  hiv_g1e1_long_nobi$prevalence <- rbind(hiv_g1e1_long_nobi$prevalence, current$prevalence[burnin_g1e1_long_nobi:iter,])
  hiv_g1e1_long_nobi$sigma <- c(hiv_g1e1_long_nobi$sigma, current$sigma[burnin_g1e1_long_nobi:iter])
  hiv_g1e1_long_nobi$reporting_prob <- c(hiv_g1e1_long_nobi$reporting_prob, current$reporting_prob[burnin_g1e1_long_nobi:iter])
  hiv_g1e1_long_nobi$x0 <- c(hiv_g1e1_long_nobi$x0, current$x0[burnin_g1e1_long_nobi:iter])
  hiv_g1e1_long_nobi$accept <- c(hiv_g1e1_long_nobi$accept, current$accept[burnin_g1e1_long_nobi:iter])
  hiv_g1e1_long_nobi$run_time <- c(hiv_g1e1_long_nobi$run_time, current$run_time)
  hiv_g1e1_long_nobi$smc_llik <- c(hiv_g1e1_long_nobi$smc_llik, current$smc_llik[burnin_g1e1_long_nobi:iter])
  hiv_g1e1_long_nobi$chain_no <- c(hiv_g1e1_long_nobi$chain_no, rep(i, iter-burnin_g1e1_long_nobi))
  hiv_g1e1_long_nobi$n_particles <- c(hiv_g1e1_long_nobi$n_particles, current$n_particles)
}


#Trace plots
# height <- 540 #best height for jpeg/png
library(scales)
iter <- 100000
height <- 3*2.5
width <- height
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_trace.pdf", width=width, height=height)
par(mfrow=c(3,3))
par(mar=c(4,4,2,1))
plot(hiv_g0e1$reporting_prob, type="l", ylab="Reporting probability", xlab="Iteration (10^3)", xaxt="n", yaxt="n", main="Epidemic only", ylim=c(0.2, 1), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,iter/1000,by=25), scientific=F))
axis(2, at=seq(0.2,1,by=0.2))
plot(hiv_g1e0$reporting_prob, type="l", ylab="Reporting probability", xlab="Iteration (10^3)", xaxt="n", main="Genetic only", col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,iter/1000,by=25), scientific=F))
plot(hiv_g1e1_long_nobi$reporting_prob, type="l", ylab="Reporting probability", xlab="Iteration (10^3)", xaxt="n", main="Combination", ylim=c(0.92,1), xlim=c(0,iter+25000), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter+25000,by=25000), labels=format(seq(0,iter/1000+25,by=25), scientific=F))
abline(v=1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==1))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==3))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==4))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==5))+1, lty=2)
plot(hiv_g0e1$sigma, type="l", ylab="Sigma", xlab="Iteration (10^3)", xaxt="n", main="Epidemic only", ylim=c(0, 0.05), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,iter/1000,by=25), scientific=F))
plot(hiv_g1e0$sigma, type="l", ylab="Sigma", xlab="Iteration (10^3)", xaxt="n", main="Genetic only", ylim=c(0, 0.08), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,iter/1000,by=25), scientific=F))
plot(hiv_g1e1_long_nobi$sigma, type="l", ylab="Sigma", xlab="Iteration (10^3)", xaxt="n", main="Combination", ylim=c(0.02, 0.2), xlim=c(0,iter+25000), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter+25000,by=25000), labels=format(seq(0,iter/1000+25,by=25), scientific=F))
abline(v=1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==1))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==3))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==4))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==5))+1, lty=2)
plot(hiv_g0e1$x0, type="l", ylab="Day 0 prevalence", xlab="Iteration (10^3)", xaxt="n", yaxt="n", main="Epidemic only", ylim=c(0, 150000), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,iter/1000,by=25), scientific=F))
axis(2, at=seq(0,150000,by=50000))
plot(hiv_g1e0$x0, type="l", ylab="Day 0 prevalence", xlab="Iteration (10^3)", xaxt="n", main="Genetic only", ylim=c(2500, 4000), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,iter/1000,by=25), scientific=F))
plot(hiv_g1e1_long_nobi$x0, type="l", ylab="Day 0 prevalence", xlab="Iteration (10^3)", xaxt="n", main="Combination", ylim=c(3000, 4500), xlim=c(0,iter+25000), col=alpha(colour="black",alpha=0.75))
axis(1, at=seq(0,iter+25000,by=25000), labels=format(seq(0,iter/1000+25,by=25), scientific=F))
abline(v=1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==1))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==3))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==4))+1, lty=2)
abline(v=max(which(hiv_g1e1_long_nobi$chain_no==5))+1, lty=2)
dev.off()

width <- 10.2
height <- width/2
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_data.pdf", width=width, height=height)
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,1.5,1.5))
library(ape)
plot(prev_paper[-1,], type="l", ylab="Observed prevalence", xlim=c(2008,2020), ylim=c(0,40000))
par(mar=c(4.5,1.5,1.5,1.5))
ptree_raw$root.time <- 2019 - max(distToRoot(ptree_raw))
plot(ptree_raw, show.tip.label = F)
axisPhylo(backward=F)
dev.off()

burn_in <- iter/5

width <- 7.5
height <- width*2/3
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_inference.pdf", width=width, height=height)
layout.matrix <- matrix(c(1,1,2,2,3,3,4,4,4,5,5,5),nrow=2,ncol=6,byrow=T)
layout(mat=layout.matrix, heights=c(2,2),widths=c(1,1,1,1,1,1))
par(mar=c(4.5,4.5,1.5,1.5))
plot(density(hiv_g1e1_long$reporting_prob), type="l", xlab="Reporting probability", main="", xlim=c(0.94, 1), ylim=c(0,200))
abline(v=mean(hiv_g1e1_long$reporting_prob), lty=2)
# lines(seq(0,1,length.out=1000), rep(1,1000), col="red")
# abline(v=0.5, lty=2, col="red")
title("A", adj=0)
plot(density(hiv_g1e1_long$sigma), type="l", xlab="Sigma", main="", xlim=c(0,0.11), ylim=c(0,60))
abline(v=mean(hiv_g1e1_long$sigma), lty=2)
# lines(seq(0,0.3,length.out=1000), dexp(seq(0,0.3,length.out=1000), rate=10), col="red")
# abline(v=0.1, lty=2, col="red")
title("B", adj=0)
hist(hiv_g1e1_long$x0, breaks=seq(3000,5000,by=50), freq=FALSE, xlab="Day 0 prevalence", main="", ylim=c(0,0.004), yaxt="n")
axis(side=2, at=seq(0,0.004,by=0.001))
abline(v=mean(hiv_g1e1_long$x0), lty=2)
# lines(seq(0,25,length.out=1000), rep(0,1000), col="red")
title("C", adj=0)
par(mar=c(4.5,4.5,1.5,1.5))
plot(1981:2019, apply(hiv_g1e1_long$birth_rate/death_rate, 2, mean), type="l", ylim=c(0,5), xlab="Year", ylab="R(t)", yaxt="n")
axis(side=2, at=seq(0,5,by=1))
lines(1981:2019, matrixStats::colQuantiles(hiv_g1e1_long$birth_rate, probs=0.025)/death_rate, lty=2)
lines(1981:2019, matrixStats::colQuantiles(hiv_g1e1_long$birth_rate, probs=0.975)/death_rate, lty=2)
# abline(h=1, col="red")
title("D", adj=0)
par(mar=c(4.5,4.5,1.5,1.5))
plot(1980:2019, matrixStats::colQuantiles(hiv_g1e1_long$prevalence, probs=0.975), lty=2, type="l", yaxt="n", xlab="Year", ylab="Prevalence", ylim=c(0,40000))
axis(2, at=seq(0,40000,by=10000), labels=format(seq(0,40000,by=10000), scientific=F))
lines(1980:2019, matrixStats::colQuantiles(hiv_g1e1_long$prevalence, probs=0.025), lty=2)
lines(1980:2019, apply(hiv_g1e1_long$prevalence, 2, mean))
title("E", adj=0)
dev.off()

100*mean(hiv_g1e1_long$reporting_prob)
mean(hiv_g1e1_long$x0)
mean(hiv_g1e1_long$sigma)
2*mean(hiv_g1e1_long$sigma)/death_rate

quantile(hiv_g1e1_long$x0, c(0.025,0.975))

lcl <- matrixStats::colQuantiles(hiv_g1e1_long$birth_rate, probs=0.025)/death_rate
ucl <- matrixStats::colQuantiles(hiv_g1e1_long$birth_rate, probs=0.975)/death_rate
df <- data.frame("Year"=1981:2019, "LCL"=lcl, "UCL"=ucl)
mean(df[1:29,3]-df[1:29,2])
mean(df[30:39,3]-df[30:39,2])

### epiestim comparison ###

library(EpiEstim)
prev_epiestim <- prev_paper
names(prev_epiestim) <- c("dates", "I")
prev_epiestim$dates <- as.numeric(prev_epiestim$dates)
#prev_epiestim$dates <- lubridate::ymd(prev_epiestim$dates, truncated = 2L)
T <- nrow(prev_epiestim)
t_start <- seq(2, T-0) #1-year window?
t_end <- t_start + 0
si_sample <- data.frame("Time"=seq(0,10,by=1),"Frequency"=c(0,0.34,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.26))
#si_sample <- t(si_sample[,2])
config <- make_config(list(t_start=t_start, t_end=t_end, n2=1))
epiestimR <- estimate_R(prev_epiestim, method="si_from_sample", si_sample=matrix(rep(si_sample$Freq,2), ncol=2), config=config)
plot(epiestimR)

par(mfrow=c(1,3))
par(mar=c(4.5, 4.5, 3, 1))
plot(2010:2019, epiestimR$R$`Mean(R)`, type="l", ylim=c(0,4), xlab="Year", ylab="R(t)", main="EpiEstim")
lines(2010:2019, epiestimR$R$`Quantile.0.025(R)`, lty=2)
lines(2010:2019, epiestimR$R$`Quantile.0.975(R)`, lty=2)


epiestim_df <- data.frame("Year"=c(2010:2019, rep(2010:2019,2)),
                          "Mean"=c(epiestimR$R$`Mean(R)`, apply(hiv_g0e1$birth_rate[-(1:burn_in),], 2, mean)/death_rate, apply(hiv_g1e1_long$birth_rate[-(1:burn_in), 30:39], 2, mean)/death_rate),
                          "LCL"=c(epiestimR$R$`Quantile.0.025(R)`, matrixStats::colQuantiles(hiv_g0e1$birth_rate[-(1:burn_in),], probs=0.025)/death_rate, matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),30:39], probs=0.025)/death_rate),
                          "UCL"=c(epiestimR$R$`Quantile.0.975(R)`, matrixStats::colQuantiles(hiv_g0e1$birth_rate[-(1:burn_in),], probs=0.975)/death_rate, matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),30:39], probs=0.975)/death_rate),
                          "Method"=c(rep("EpiEstim",10), rep("EpiSky - Epidemic only",10), rep("EpiSky - Combination",10)))

# pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_epiestim.pdf", width=8, height=4)
# library(ggplot2)
# ggplot() +
#   geom_line(data=epiestim_df, aes(x=Year, y=Mean, color=Method)) +
#   geom_ribbon(data=epiestim_df, aes(x=Year, ymin=LCL, ymax=UCL, fill=Method), alpha=0.2) +
#   labs(y="R(t)") +
#   scale_x_continuous(breaks=seq(1997,2019,by=2)) +
#   scale_fill_grey(start=0.7, end=0.2) +
#   scale_color_grey(start=0.7, end=0.2) +
#   theme_light()
# dev.off()

height <- 2.5
width <- height*3
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_epiestim.pdf", width=width, height=height)
par(mfrow=c(1,3))
par(mar=c(4.5,4.5,1.5,1.5))
plot(2010:2019, epiestim_df$Mean[epiestim_df$Method=="EpiEstim"], ylim=c(1,3.5), xlim=c(2010,2020), type="l", main="EpiEstim", xlab="Year", ylab="R(t)")
lines(2010:2019, epiestim_df$LCL[epiestim_df$Method=="EpiEstim"], lty=2)
lines(2010:2019, epiestim_df$UCL[epiestim_df$Method=="EpiEstim"], lty=2)
plot(2010:2019, epiestim_df$Mean[epiestim_df$Method=="EpiSky - Epidemic only"], ylim=c(1,3.5), xlim=c(2010,2020), type="l", main="EpiSky - Epidemic only", xlab="Year", ylab="R(t)")
lines(2010:2019, epiestim_df$LCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
lines(2010:2019, epiestim_df$UCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
plot(2010:2019, epiestim_df$Mean[epiestim_df$Method=="EpiSky - Combination"], ylim=c(1,3.5), xlim=c(2010,2020), type="l", main="EpiSky - Combination", xlab="Year", ylab="R(t)")
lines(2010:2019, epiestim_df$LCL[epiestim_df$Method=="EpiSky - Combination"], lty=2)
lines(2010:2019, epiestim_df$UCL[epiestim_df$Method=="EpiSky - Combination"], lty=2)
dev.off()

### skygrowth comparison ###

library(skygrowth)
set.seed(1)
r <- skygrowth.mcmc(ptree_raw, mhsteps=100000, res=72)
r1 <- computeR(r, gamma=death_rate)

skygrowth_rt_df3 <- data.frame("Year"=c(r1$time+2019, 1981:2019, 1981:2019),
                               "Mean"=c(apply(r1$R,2,mean), apply(hiv_g1e0$birth_rate[-(1:burn_in),], 2, mean)/death_rate, apply(hiv_g1e1_long$birth_rate[-(1:burn_in),], 2, mean)/death_rate),
                               "LCL"=c(r1$R_ci[,1], matrixStats::colQuantiles(hiv_g1e0$birth_rate[-(1:burn_in),], probs=0.025)/death_rate, matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),], probs=0.025)/death_rate),
                               "UCL"=c(r1$R_ci[,3], matrixStats::colQuantiles(hiv_g1e0$birth_rate[-(1:burn_in),], probs=0.975)/death_rate, matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),], probs=0.975)/death_rate),
                               "Method"=c(rep("SkyGrowth", nrow(r1$R_ci)), rep("EpiSky - Genetic only", length(1981:2019)), rep("EpiSky - Combination", length(1981:2019))))

# pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_skygrowth.pdf", width=8, height=4)
# library(ggplot2)
# ggplot(skygrowth_rt_df3) + geom_line(aes(x=Year, y=Mean, colour=Method)) +
#   geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Year, fill=Method), alpha=0.3) +
#   ylab("R(t)")
# dev.off()

height <- 2.5
width <- height*3
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_skygrowth.pdf", width=width, height=height)
par(mfrow=c(1,3))
par(mar=c(4.5,4.5,1.5,1.5))
plot(r1$time+2019, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="SkyGrowth"], ylim=c(0,6.3), xlim=c(1980,2010), type="l", main="SkyGrowth", xlab="Year", ylab="R(t)")
lines(r1$time+2019, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
lines(r1$time+2019, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
plot(1981:2019, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], ylim=c(0,6.3), xlim=c(1980,2010), type="l", main="EpiSky - Genetic only", xlab="Year", ylab="R(t)")
lines(1981:2019, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
lines(1981:2019, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
plot(1981:2019, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Combination"], ylim=c(0,6.3), xlim=c(1980,2010), type="l", main="EpiSky - Combination", xlab="Year", ylab="R(t)")
lines(1981:2019, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
lines(1981:2019, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
dev.off()

comparison_df <- data.frame("Year"=c(r1$time+2019, rep(1981:2019,2), 2010:2019),
                          "Mean"=c(apply(r1$R,2,mean),
                                   apply(hiv_g1e0$birth_rate[-(1:burn_in),], 2, mean)/death_rate,
                                   apply(hiv_g1e1_long$birth_rate[-(1:burn_in),], 2, mean)/death_rate,
                                   apply(hiv_g0e1$birth_rate[-(1:burn_in),], 2, mean)/death_rate),
                          "LCL"=c(r1$R_ci[,1],
                                  matrixStats::colQuantiles(hiv_g1e0$birth_rate[-(1:burn_in),], probs=0.025)/death_rate,
                                  matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),], probs=0.025)/death_rate,
                                  matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),30:39], probs=0.025)/death_rate),
                          "UCL"=c(r1$R_ci[,3],
                                  matrixStats::colQuantiles(hiv_g1e0$birth_rate[-(1:burn_in),], probs=0.975)/death_rate,
                                  matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),], probs=0.975)/death_rate,
                                  matrixStats::colQuantiles(hiv_g1e1_long$birth_rate[-(1:burn_in),30:39], probs=0.975)/death_rate),
                          "Method"=c(rep("SkyGrowth", nrow(r1$R_ci)),
                                     rep("EpiSky - Genetic only", length(1981:2019)),
                                     rep("EpiSky - Combination", length(1981:2019)),
                                     rep("EpiSky - Epidemic only", length(2010:2019))))

# pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_comparison.pdf", width=8, height=4)
# library(ggplot2)
# ggplot(comparison_df) + geom_line(aes(x=Year, y=Mean, colour=Method)) +
#   geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Year, fill=Method), alpha=0.3) +
#   ylab("R(t)") +
#   scale_color_hue(labels=c("SkyGrowth", "EpiEstim", "EpiSky - Epidemic only", "EpiSky - Genetic only", "EpiSky - Combination")) +
#   scale_fill_hue(labels=c("SkyGrowth", "EpiEstim", "EpiSky - Epidemic only", "EpiSky - Genetic only", "EpiSky - Combination")) +
#   scale_x_continuous(breaks=seq(1981,2010,by=5), minor_breaks=seq(1981,2010,by=5))
# dev.off()

height <- 2.5
width <- height*5
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/hiv_comparison.pdf", width=width, height=height)
par(mfrow=c(1,5))
par(mar=c(4.5,4.5,1.5,1.5))
plot(r1$time+2019, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="SkyGrowth"], ylim=c(0,6.3), xlim=c(1980,2020), type="l", main="SkyGrowth", xlab="Year", ylab="R(t)")
lines(r1$time+2019, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
lines(r1$time+2019, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
plot(2010:2019, epiestim_df$Mean[epiestim_df$Method=="EpiEstim"], ylim=c(0,6.3), xlim=c(1980,2020), type="l", main="EpiEstim", xlab="Year", ylab="R(t)")
lines(2010:2019, epiestim_df$LCL[epiestim_df$Method=="EpiEstim"], lty=2)
lines(2010:2019, epiestim_df$UCL[epiestim_df$Method=="EpiEstim"], lty=2)
plot(2010:2019, epiestim_df$Mean[epiestim_df$Method=="EpiSky - Epidemic only"], ylim=c(0,6.3), xlim=c(1980,2020), type="l", main="EpiSky - Epidemic only", xlab="Year", ylab="R(t)")
lines(2010:2019, epiestim_df$LCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
lines(2010:2019, epiestim_df$UCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
plot(1981:2019, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], ylim=c(0,6.3), xlim=c(1980,2020), type="l", main="EpiSky - Genetic only", xlab="Year", ylab="R(t)")
lines(1981:2019, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
lines(1981:2019, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
plot(1981:2019, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Combination"], ylim=c(0,6.3), xlim=c(1980,2020), type="l", main="EpiSky - Combination", xlab="Year", ylab="R(t)")
lines(1981:2019, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
lines(1981:2019, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
dev.off()


library(coda)
mcmc_g0e1 <- mcmc(matrix(data=c(hiv_g0e1$sigma, hiv_g0e1$reporting_prob, hiv_g0e1$x0, hiv_g0e1$birth_rate, hiv_g0e1$prevalence[,-1]),
                         ncol=23,
                         dimnames=list(1:100000,c("sigma","pobs","x0",paste0("br",1:10),paste0("prev",1:10)))))
mcmc_g1e0 <- mcmc(matrix(data=c(hiv_g1e0$sigma, hiv_g1e0$reporting_prob, hiv_g1e0$x0, hiv_g1e0$birth_rate, hiv_g1e0$prevalence[,-1]),
                         ncol=81,
                         dimnames=list(1:100000,c("sigma","pobs","x0",paste0("br",1:39),paste0("prev",1:39)))))
mcmc_g1e1 <- mcmc(matrix(data=c(hiv_g1e1_long$sigma, hiv_g1e1_long$reporting_prob, hiv_g1e1_long$x0, hiv_g1e1_long$birth_rate, hiv_g1e1_long$prevalence[,-1]),
                         ncol=81,
                         dimnames=list(1:80000,c("sigma","pobs","x0",paste0("br",1:39),paste0("prev",1:39)))))

n1 <- which.max(is.na(hiv_g1e1_long_1$sigma))-1
mcmc_g1e1_1 <- mcmc(matrix(data=c(hiv_g1e1_long_1$sigma[1:n1], hiv_g1e1_long_1$reporting_prob[1:n1], hiv_g1e1_long_1$x0[1:n1], hiv_g1e1_long_1$birth_rate[1:n1,], hiv_g1e1_long_1$prevalence[1:n1,-1]),
                           ncol=81,
                           dimnames=list(1:n1,c("sigma","pobs","x0",paste0("br",1:39),paste0("prev",1:39)))))
n2 <- which.max(is.na(hiv_g1e1_long_2$sigma))-1
mcmc_g1e1_2 <- mcmc(matrix(data=c(hiv_g1e1_long_2$sigma[1:n2], hiv_g1e1_long_2$reporting_prob[1:n2], hiv_g1e1_long_2$x0[1:n2], hiv_g1e1_long_2$birth_rate[1:n2,], hiv_g1e1_long_2$prevalence[1:n2,-1]),
                           ncol=81,
                           dimnames=list(1:n2,c("sigma","pobs","x0",paste0("br",1:39),paste0("prev",1:39)))))
n3 <- which.max(is.na(hiv_g1e1_long_3$sigma))-1
mcmc_g1e1_3 <- mcmc(matrix(data=c(hiv_g1e1_long_3$sigma[1:n3], hiv_g1e1_long_3$reporting_prob[1:n3], hiv_g1e1_long_3$x0[1:n3], hiv_g1e1_long_3$birth_rate[1:n3,], hiv_g1e1_long_3$prevalence[1:n3,-1]),
                           ncol=81,
                           dimnames=list(1:n3,c("sigma","pobs","x0",paste0("br",1:39),paste0("prev",1:39)))))
n4 <- which.max(is.na(hiv_g1e1_long_4$sigma))-1
mcmc_g1e1_4 <- mcmc(matrix(data=c(hiv_g1e1_long_4$sigma[1:n4], hiv_g1e1_long_4$reporting_prob[1:n4], hiv_g1e1_long_4$x0[1:n4], hiv_g1e1_long_4$birth_rate[1:n4,], hiv_g1e1_long_4$prevalence[1:n4,-1]),
                           ncol=81,
                           dimnames=list(1:n4,c("sigma","pobs","x0",paste0("br",1:39),paste0("prev",1:39)))))
n5 <- which.max(is.na(hiv_g1e1_long_5$sigma))-1
mcmc_g1e1_5 <- mcmc(matrix(data=c(hiv_g1e1_long_5$sigma[1:n5], hiv_g1e1_long_5$reporting_prob[1:n5], hiv_g1e1_long_5$x0[1:n5], hiv_g1e1_long_5$birth_rate[1:n5,], hiv_g1e1_long_5$prevalence[1:n5,-1]),
                           ncol=81,
                           dimnames=list(1:n5,c("sigma","pobs","x0",paste0("br",1:39),paste0("prev",1:39)))))

mcmc_g1e1_trim_1 <- mcmc(mcmc_g1e1_1[1:26000,])
mcmc_g1e1_trim_2 <- mcmc(mcmc_g1e1_2[1:26000,])
mcmc_g1e1_trim_3 <- mcmc(mcmc_g1e1_3[1:26000,])
mcmc_g1e1_trim_4 <- mcmc(mcmc_g1e1_4[1:26000,])
mcmc_g1e1_trim_5 <- mcmc(mcmc_g1e1_5[1:26000,])

mcmc_g1e1_trim2_1 <- mcmc(mcmc_g1e1_1[(n1-19999):n1,])
mcmc_g1e1_trim2_2 <- mcmc(mcmc_g1e1_2[(n2-19999):n2,])
mcmc_g1e1_trim2_3 <- mcmc(mcmc_g1e1_3[(n3-19999):n3,])
mcmc_g1e1_trim2_4 <- mcmc(mcmc_g1e1_4[(n4-19999):n4,])
mcmc_g1e1_trim2_5 <- mcmc(mcmc_g1e1_5[(n5-19999):n5,])


mcmc_g1e1_trim_list <- mcmc.list(mcmc_g1e1_trim_1, mcmc_g1e1_trim_3, mcmc_g1e1_trim_4, mcmc_g1e1_trim_5)
mcmc_g1e1 <- mcmc.list(mcmc_g1e1_trim2_1, mcmc_g1e1_trim2_3, mcmc_g1e1_trim2_4, mcmc_g1e1_trim2_5)

plot(mcmc_g1e1)
gelman.diag(mcmc_g1e1)

gelman.diag(mcmc_g1e1)
summary(mcmc_g1e1)
1 - mean(rejectionRate(mcmc_g1e1))
effectiveSize(mcmc_g1e1)
autocorr.plot(mcmc_g1e1)
plotEssBurn(mcmc_g1e1)

library(rstan)
rhat_g0e1 <- rep(NA,3+2*ncol(hiv_g0e1$birth_rate))
rhat_g1e0 <- rep(NA,2+2*ncol(hiv_g1e0$birth_rate))
rhat_g1e1 <- rep(NA,3+2*ncol(hiv_g1e1_long$birth_rate))
rhat_g1e1_all <- rep(NA,3+2*ncol(hiv_g1e1_long_1$birth_rate))

rhat_g0e1[1] <- Rhat(hiv_g0e1$sigma)
rhat_g0e1[2] <- Rhat(hiv_g0e1$reporting_prob)
rhat_g0e1[3] <- Rhat(hiv_g0e1$x0)
for (i in 1:((length(rhat_g0e1)-3)/2)) {
  rhat_g0e1[3+i] <- Rhat(hiv_g0e1$birth_rate[,i])
  rhat_g0e1[3+((length(rhat_g0e1)-3)/2)+i] <- Rhat(hiv_g0e1$prevalence[,i+1])
}
rhat_g1e0[1] <- Rhat(hiv_g1e0$sigma)
rhat_g1e0[2] <- Rhat(hiv_g1e0$x0)
for (i in 1:((length(rhat_g1e0)-2)/2)) {
  rhat_g1e0[2+i] <- Rhat(hiv_g1e0$birth_rate[,i])
  rhat_g1e0[2+((length(rhat_g1e0)-2)/2)+i] <- Rhat(hiv_g1e0$prevalence[,i+1])
}
rhat_g1e1[1] <- Rhat(hiv_g1e1_long$sigma)
rhat_g1e1[2] <- Rhat(hiv_g1e1_long$reporting_prob)
rhat_g1e1[3] <- Rhat(hiv_g1e1_long$x0)
for (i in 1:((length(rhat_g1e1)-3)/2)) {
  rhat_g1e1[3+i] <- Rhat(hiv_g1e1_long$birth_rate[,i])
  rhat_g1e1[3+((length(rhat_g1e1)-3)/2)+i] <- Rhat(hiv_g1e1_long$prevalence[,i+1])
}

max(rhat_g0e1)
max(rhat_g1e0)
max(rhat_g1e1)

par(mfrow=c(1,3))
plot(rhat_g0e1, ylim=c(0.99999, 1.04))
abline(h=1.01, lty=2)
plot(rhat_g1e0, ylim=c(0.99999, 1.04))
abline(h=1.01, lty=2)
plot(rhat_g1e1, ylim=c(0.99999, 1.04))
abline(h=1.01, lty=2)

sum(rhat_g1e1<1.01)
sum(rhat_g1e1>=1.01)

rhat_g1e1_all <- c()
rhat_g1e1_all <- c(rhat_g1e1_all,
                   Rhat(matrix(c(hiv_g1e1_long_1$sigma,
                                 hiv_g1e1_long_3$sigma,
                                 hiv_g1e1_long_4$sigma,
                                 hiv_g1e1_long_5$sigma), ncol=4)[1:26000,]))
rhat_g1e1_all <- c(rhat_g1e1_all,
                   Rhat(matrix(c(hiv_g1e1_long_1$reporting_prob,
                                 hiv_g1e1_long_3$reporting_prob,
                                 hiv_g1e1_long_4$reporting_prob,
                                 hiv_g1e1_long_5$reporting_prob), ncol=4)[1:26000,]))
rhat_g1e1_all <- c(rhat_g1e1_all,
                   Rhat(matrix(c(hiv_g1e1_long_1$x0,
                                 hiv_g1e1_long_3$x0,
                                 hiv_g1e1_long_4$x0,
                                 hiv_g1e1_long_5$x0), ncol=4)[1:26000,]))
for (i in 1:39) {
  rhat_g1e1_all <- c(rhat_g1e1_all,
                     Rhat(matrix(c(hiv_g1e1_long_1$birth_rate[,i],
                                   hiv_g1e1_long_3$birth_rate[,i],
                                   hiv_g1e1_long_4$birth_rate[,i],
                                   hiv_g1e1_long_5$birth_rate[,i]), ncol=4)[1:26000,])
  )
}
for (i in 1:39) {
  rhat_g1e1_all <- c(rhat_g1e1_all,
                     Rhat(matrix(c(hiv_g1e1_long_1$prevalence[,i+1],
                                   hiv_g1e1_long_3$prevalence[,i+1],
                                   hiv_g1e1_long_4$prevalence[,i+1],
                                   hiv_g1e1_long_5$prevalence[,i+1]), ncol=4)[1:26000,])
  )
}
plot(rhat_g1e1_all)

ess_bulk(matrix(c(hiv_g1e1_long_1$sigma,
                  hiv_g1e1_long_3$sigma,
                  hiv_g1e1_long_4$sigma,
                  hiv_g1e1_long_5$sigma), ncol=4)[1:26000,])
ess_bulk(matrix(c(hiv_g1e1_long_1$reporting_prob,
                  hiv_g1e1_long_3$reporting_prob,
                  hiv_g1e1_long_4$reporting_prob,
                  hiv_g1e1_long_5$reporting_prob), ncol=4)[1:26000,])
ess_bulk(matrix(c(hiv_g1e1_long_1$x0,
                  hiv_g1e1_long_3$x0,
                  hiv_g1e1_long_4$x0,
                  hiv_g1e1_long_5$x0), ncol=4)[1:26000,])
