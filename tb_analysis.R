files <- list.files(path="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/RData/", pattern="(tb)")
for (i in 1:length(files)) {
  print(i)
  load(paste0("~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/RData/", files[i]))
}
gc()

devtools::load_all()

accept_g0e1 <- c(tb_g0e1_1$acceptance_rate, tb_g0e1_2$acceptance_rate, tb_g0e1_3$acceptance_rate, tb_g0e1_4$acceptance_rate, tb_g0e1_5$acceptance_rate)*100
accept_g1e0 <- c(tb_g1e0_1$acceptance_rate, tb_g1e0_2$acceptance_rate, tb_g1e0_3$acceptance_rate, tb_g1e0_4$acceptance_rate, tb_g1e0_5$acceptance_rate)*100
accept_g1e1_short <- c(tb_g1e1_short_1$acceptance_rate, tb_g1e1_short_2$acceptance_rate, tb_g1e1_short_3$acceptance_rate, tb_g1e1_short_4$acceptance_rate, tb_g1e1_short_5$acceptance_rate)*100
accept_g1e1_long <- c(tb_g1e1_long_1$acceptance_rate, tb_g1e1_long_2$acceptance_rate, tb_g1e1_long_3$acceptance_rate, tb_g1e1_long_4$acceptance_rate, tb_g1e1_long_5$acceptance_rate)*100

burn_in <- 10000
smcvar_g0e1 <- c(var(tb_g0e1_1$smc_llik[-(1:burn_in)]),
                 var(tb_g0e1_2$smc_llik[-(1:burn_in)]),
                 var(tb_g0e1_3$smc_llik[-(1:burn_in)]),
                 var(tb_g0e1_4$smc_llik[-(1:burn_in)]),
                 var(tb_g0e1_5$smc_llik[-(1:burn_in)]))
smcvar_g1e0 <- c(var(tb_g1e0_1$smc_llik[-(1:burn_in)]),
                 var(tb_g1e0_2$smc_llik[-(1:burn_in)]),
                 var(tb_g1e0_3$smc_llik[-(1:burn_in)]),
                 var(tb_g1e0_4$smc_llik[-(1:burn_in)]),
                 var(tb_g1e0_5$smc_llik[-(1:burn_in)]))
smcvar_g1e1_long <- c(var(tb_g1e1_long_1$smc_llik[-(1:burn_in)]),
                      var(tb_g1e1_long_2$smc_llik[-(1:burn_in)]),
                      var(tb_g1e1_long_3$smc_llik[-(1:burn_in)]),
                      var(tb_g1e1_long_4$smc_llik[-(1:burn_in)]),
                      var(tb_g1e1_long_5$smc_llik[-(1:burn_in)]))

par(mfrow=c(5,3))
par(mar=c(4,4,1,1))
plot(tb_g0e1_1$sigma, type="l")
plot(tb_g0e1_1$reporting_prob, type="l")
plot(tb_g0e1_1$x0, type="l")
plot(tb_g0e1_2$sigma, type="l")
plot(tb_g0e1_2$reporting_prob, type="l")
plot(tb_g0e1_3$x0, type="l")
plot(tb_g0e1_3$sigma, type="l")
plot(tb_g0e1_3$reporting_prob, type="l")
plot(tb_g0e1_3$x0, type="l")
plot(tb_g0e1_4$sigma, type="l")
plot(tb_g0e1_4$reporting_prob, type="l")
plot(tb_g0e1_4$x0, type="l")
plot(tb_g0e1_5$sigma, type="l")
plot(tb_g0e1_5$reporting_prob, type="l")
plot(tb_g0e1_5$x0, type="l")
accept_g0e1
smcvar_g0e1

par(mfrow=c(5,2))
par(mar=c(4,4,1,1))
plot(tb_g1e0_1$sigma, type="l")
plot(tb_g1e0_1$x0, type="l")
plot(tb_g1e0_2$sigma, type="l")
plot(tb_g1e0_3$x0, type="l")
plot(tb_g1e0_3$sigma, type="l")
plot(tb_g1e0_3$x0, type="l")
plot(tb_g1e0_4$sigma, type="l")
plot(tb_g1e0_4$x0, type="l")
plot(tb_g1e0_5$sigma, type="l")
plot(tb_g1e0_5$x0, type="l")
accept_g1e0
smcvar_g1e0

par(mfrow=c(5,3))
par(mar=c(4,4,1,1))
plot(tb_g1e1_long_1$sigma, type="l")
plot(tb_g1e1_long_1$reporting_prob, type="l")
plot(tb_g1e1_long_1$x0, type="l")
plot(tb_g1e1_long_2$sigma, type="l")
plot(tb_g1e1_long_2$reporting_prob, type="l")
plot(tb_g1e1_long_3$x0, type="l")
plot(tb_g1e1_long_3$sigma, type="l")
plot(tb_g1e1_long_3$reporting_prob, type="l")
plot(tb_g1e1_long_3$x0, type="l")
plot(tb_g1e1_long_4$sigma, type="l")
plot(tb_g1e1_long_4$reporting_prob, type="l")
plot(tb_g1e1_long_4$x0, type="l")
plot(tb_g1e1_long_5$sigma, type="l")
plot(tb_g1e1_long_5$reporting_prob, type="l")
plot(tb_g1e1_long_5$x0, type="l")
accept_g1e1_long
smcvar_g1e1_long

tb_g0e1 <- tb_g0e1_2
tb_g1e0 <- tb_g1e0_2
tb_g1e1 <- tb_g1e1_long_4

tb_g1e1$n_particles
tb_g1e1$run_time/60/60
tb_g1e1$acceptance_rate

width <- 10.2
height <- width/2
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_data.pdf", width=width, height=height)
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,1.5,1.5))
library(ape)
plot(prev_raw, type="l", ylab="Observed prevalence", xlim=c(1995,2010), ylim=c(0,40))
par(mar=c(4.5,1.5,1.5,1.5))
ptree_raw$root.time <- 2009 - max(distToRoot(ptree_raw))
plot(ptree_raw, show.tip.label = F)
axisPhylo(backward=F)
dev.off()

#Trace plots
# height <- 540 #best height for jpeg/png
library(scales)
height <- 7.5
width <- height
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_trace.pdf", width=width, height=height)
par(mfrow=c(3,3))
par(mar=c(4,4,2,1))
plot(tb_g0e1$reporting_prob, type="l", ylab="Reporting probability", xlab="Iteration (10^3)", xaxt="n", yaxt="n", main="Epidemic only", ylim=c(0, 0.25), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
axis(2, at=seq(0,0.25,by=0.05))
plot(tb_g1e0$reporting_prob, type="l", ylab="Reporting probability", xlab="Iteration (10^3)", xaxt="n", main="Genetic only", col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
plot(tb_g1e1$reporting_prob, type="l", ylab="Reporting probability", xlab="Iteration (10^3)", xaxt="n", main="Combination", ylim=c(0,0.12), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
plot(tb_g0e1$sigma, type="l", ylab="Sigma", xlab="Iteration (10^3)", xaxt="n", main="Epidemic only", ylim=c(0, 1.25), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
plot(tb_g1e0$sigma, type="l", ylab="Sigma", xlab="Iteration (10^3)", xaxt="n", main="Genetic only", ylim=c(0, 0.2), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
plot(tb_g1e1$sigma, type="l", ylab="Sigma", xlab="Iteration (10^3)", xaxt="n", main="Combination", ylim=c(0, 0.3), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
plot(tb_g0e1$x0, type="l", ylab="Day 0 prevalence", xlab="Iteration (10^3)", xaxt="n", yaxt="n", main="Epidemic only", ylim=c(0, 4000), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
axis(2, at=seq(0,4000,by=1000))
plot(tb_g1e0$x0, type="l", ylab="Day 0 prevalence", xlab="Iteration (10^3)", xaxt="n", main="Genetic only", ylim=c(0, 30), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
plot(tb_g1e1$x0, type="l", ylab="Day 0 prevalence", xlab="Iteration (10^3)", xaxt="n", main="Combination", ylim=c(0, 50), col=alpha("black", 0.75))
axis(1, at=seq(0,iter,by=25000), labels=format(seq(0,100,by=25), scientific=F))
dev.off()

width <- 10.2
height <- width/2
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_data.pdf", width=width, height=height)
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,1.5,1.5))
library(ape)
plot(prev_raw, type="l", ylab="Observed prevalence", xlim=c(1996,2010))
par(mar=c(4.5,1.5,1.5,1.5))
ptree_raw$root.time <- 1970.661
plot(ptree_raw, show.tip.label = F)
axisPhylo(backward=F)
dev.off()

burn_in <- iter/5

width <- 9
height <- width*2/3
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_inference.pdf", width=width, height=height)
layout.matrix <- matrix(c(1,1,2,2,3,3,4,4,4,5,5,5),nrow=2,ncol=6,byrow=T)
layout(mat=layout.matrix, heights=c(2,2),widths=c(1,1,1,1,1,1))
par(mar=c(4.5,4.5,1.5,1.5))
plot(density(tb_g1e1$reporting_prob[-(1:burn_in)]), type="l", xlab="Reporting probability", main="", xlim=c(0.005, 0.025), ylim=c(0,200))
abline(v=mean(tb_g1e1$reporting_prob[-(1:burn_in)]), lty=2)
# lines(seq(0,1,length.out=1000), rep(1,1000), col="red")
# abline(v=0.5, lty=2, col="red")
title("A", adj=0)
plot(density(tb_g1e1$sigma[-(1:burn_in)]), type="l", xlab="Sigma", main="", xlim=c(0,0.3), ylim=c(0,20))
abline(v=mean(tb_g1e1$sigma[-(1:burn_in)]), lty=2)
# lines(seq(0,0.3,length.out=1000), dexp(seq(0,0.3,length.out=1000), rate=10), col="red")
# abline(v=0.1, lty=2, col="red")
title("B", adj=0)
hist(tb_g1e1$x0[-(1:burn_in)], breaks=seq(0,50,by=2), freq=FALSE, xlab="Day 0 prevalence", main="", ylim=c(0,0.2), yaxt="n")
axis(side=2, at=seq(0,0.2,by=0.05))
abline(v=mean(tb_g1e1$x0[-(1:burn_in)]), lty=2)
# lines(seq(0,25,length.out=1000), rep(0,1000), col="red")
title("C", adj=0)
par(mar=c(4.5,4.5,1.5,1.5))
plot(1970:2009, apply(tb_g1e1$birth_rate[-(1:burn_in),]/death_rate, 2, mean), type="l", ylim=c(0,3), xlab="Year", ylab="R(t)", yaxt="n")
axis(side=2, at=seq(0,3,by=1))
lines(1970:2009, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], probs=0.025)/death_rate, lty=2)
lines(1970:2009, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], probs=0.975)/death_rate, lty=2)
# abline(h=1, col="red")
title("D", adj=0)
par(mar=c(4.5,4.5,1.5,1.5))
plot(matrixStats::colQuantiles(tb_g1e1$prevalence[-(1:burn_in),], probs=0.975), lty=2, type="l", xaxt="n", yaxt="n", xlab="Year", ylab="Prevalence", ylim=c(0,4000))
axis(1, at=seq(1,42,by=10), labels=seq(1969, 2009, by=10))
axis(2, at=seq(0,4000,by=1000), labels=format(seq(0,4000,by=1000), scientific=F))
lines(1:41, matrixStats::colQuantiles(tb_g1e1$prevalence[-(1:burn_in),], probs=0.025), lty=2)
lines(1:41, apply(tb_g1e1$prevalence[-(1:burn_in),], 2, mean))
title("E", adj=0)
dev.off()

100*mean(tb_g1e1$reporting_prob[-(1:burn_in)])
mean(tb_g1e1$x0[-(1:burn_in)])
mean(tb_g1e1$sigma[-(1:burn_in)])
2*mean(tb_g1e1$sigma[-(1:burn_in)])/death_rate

mean <- apply(tb_g1e1$birth_rate[-(1:burn_in),],2,mean)/death_rate
lcl <- matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], prob=0.025)/death_rate
ucl <- matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], prob=0.975)/death_rate
df <- data.frame("year"=1970:2009, mean, lcl, ucl)
df

### epiestim comparison ###

library(EpiEstim)
prev_epiestim <- prev_raw
names(prev_epiestim) <- c("dates", "I")
prev_epiestim$dates <- as.numeric(prev_epiestim$dates)
#prev_epiestim$dates <- lubridate::ymd(prev_epiestim$dates, truncated = 2L)
T <- nrow(prev_epiestim)
t_start <- seq(2, T-0) #1-year window?
t_end <- t_start + 0
epiestimR <- estimate_R(prev_epiestim, method="parametric_si", config = make_config(list(mean_si = 2, std_si = 2, t_start=t_start, t_end=t_end)))

epiestim_df <- data.frame("Year"=c(1997:2009, rep(1996:2009,2)),
                          "Mean"=c(epiestimR$R$`Mean(R)`, apply(tb_g0e1$birth_rate[-(1:burn_in),], 2, mean)/death_rate, apply(tb_g1e1$birth_rate[-(1:burn_in), 27:40], 2, mean)/death_rate),
                          "LCL"=c(epiestimR$R$`Quantile.0.025(R)`, matrixStats::colQuantiles(tb_g0e1$birth_rate[-(1:burn_in),], probs=0.025)/death_rate, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),27:40], probs=0.025)/death_rate),
                          "UCL"=c(epiestimR$R$`Quantile.0.975(R)`, matrixStats::colQuantiles(tb_g0e1$birth_rate[-(1:burn_in),], probs=0.975)/death_rate, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),27:40], probs=0.975)/death_rate),
                          "Method"=c(rep("EpiEstim",13), rep("EpiSky - Epidemic only",14), rep("EpiSky - Combination",14)))

# pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_epiestim.pdf", width=8, height=4)
# library(ggplot2)
# ggplot() +
#   geom_line(data=epiestim_df, aes(x=Year, y=Mean, color=Method)) +
#   geom_ribbon(data=epiestim_df, aes(x=Year, ymin=LCL, ymax=UCL, fill=Method), alpha=0.2) +
#   labs(y="R(t)") +
#   scale_x_continuous(breaks=seq(1997,2009,by=2)) +
#   scale_fill_grey(start=0.7, end=0.2) +
#   scale_color_grey(start=0.7, end=0.2) +
#   theme_light()
# dev.off()

height <- 2.5
width <- height*3
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_epiestim.pdf", width=width, height=height)
par(mfrow=c(1,3))
par(mar=c(4.5,4.5,1.5,1.5))
plot(1997:2009, epiestim_df$Mean[epiestim_df$Method=="EpiEstim"], ylim=c(0,30), xlim=c(1996,2010), type="l", main="EpiEstim", xlab="Year", ylab="R(t)")
lines(1997:2009, epiestim_df$LCL[epiestim_df$Method=="EpiEstim"], lty=2)
lines(1997:2009, epiestim_df$UCL[epiestim_df$Method=="EpiEstim"], lty=2)
plot(1996:2009, epiestim_df$Mean[epiestim_df$Method=="EpiSky - Epidemic only"], ylim=c(0,30), xlim=c(1996,2010), type="l", main="EpiSky - Epidemic only", xlab="Year", ylab="R(t)")
lines(1996:2009, epiestim_df$LCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
lines(1996:2009, epiestim_df$UCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
plot(1996:2009, epiestim_df$Mean[epiestim_df$Method=="EpiSky - Combination"], ylim=c(0,30), xlim=c(1996,2010), type="l", main="EpiSky - Combination", xlab="Year", ylab="R(t)")
lines(1996:2009, epiestim_df$LCL[epiestim_df$Method=="EpiSky - Combination"], lty=2)
lines(1996:2009, epiestim_df$UCL[epiestim_df$Method=="EpiSky - Combination"], lty=2)
dev.off()

### skygrowth comparison ###

library(skygrowth)
set.seed(1)
r <- skygrowth.mcmc(ptree, mhsteps=100000, res=40)
# skygrowth_ne_df2 <- data.frame("Year"=c(r$time+2009-40, 1970:2009),
#                                "Mean"=c(r$ne, apply(tb_g1e0$prevalence[-(1:50000),-1], 2, mean)),
#                                "LCL"=c(r$ne_ci[,1], matrixStats::colQuantiles(tb_g1e0$prevalence[-(1:50000),-1], probs=0.025)),
#                                "UCL"=c(r$ne_ci[,3], matrixStats::colQuantiles(tb_g1e0$prevalence[-(1:50000),-1], probs=0.975)),
#                                "Group"=c(rep("skygrowth", length(r$ne)), rep("gen only", length(1970:2009))))
# skygrowth_ne_df3 <- data.frame("Year"=c(r$time+2009-40, 1970:2009, 1970:2009),
#                                "Mean"=c(r$ne, apply(tb_g1e0$prevalence[-(1:50000),-1], 2, mean), apply(tb_g1e1$prevalence[-(1:50000),-1], 2, mean)),
#                                "LCL"=c(r$ne_ci[,1], matrixStats::colQuantiles(tb_g1e0$prevalence[-(1:50000),-1], probs=0.025), matrixStats::colQuantiles(tb_g1e1$prevalence[-(1:50000),-1], probs=0.025)),
#                                "UCL"=c(r$ne_ci[,3], matrixStats::colQuantiles(tb_g1e0$prevalence[-(1:50000),-1], probs=0.975), matrixStats::colQuantiles(tb_g1e1$prevalence[-(1:50000),-1], probs=0.975)),
#                                "Group"=c(rep("skygrowth", length(r$ne)), rep("gen only", length(1970:2009)), rep("epi+gen", length(1970:2009))))
# library(ggplot2)
# ggplot(skygrowth_ne_df2) + geom_line(aes(x=Year, y=Mean, colour=Group)) + geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Year, fill=Group), alpha=0.3) + ylab("Effective population size")
# library(ggplot2)
# ggplot(skygrowth_ne_df3) + geom_line(aes(x=Year, y=Mean, colour=Group)) + geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Year, fill=Group), alpha=0.3) + ylab("Effective population size")

r1 <- computeR(r, gamma=death_rate)

skygrowth_rt_df3 <- data.frame("Year"=c(r1$time+2009, 1970:2009, 1970:2009),
                               "Mean"=c(apply(r1$R,2,mean), apply(tb_g1e0$birth_rate[-(1:burn_in),], 2, mean)/death_rate, apply(tb_g1e1$birth_rate[-(1:burn_in),], 2, mean)/death_rate),
                               "LCL"=c(r1$R_ci[,1], matrixStats::colQuantiles(tb_g1e0$birth_rate[-(1:burn_in),], probs=0.025)/death_rate, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], probs=0.025)/death_rate),
                               "UCL"=c(r1$R_ci[,3], matrixStats::colQuantiles(tb_g1e0$birth_rate[-(1:burn_in),], probs=0.975)/death_rate, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], probs=0.975)/death_rate),
                               "Method"=c(rep("SkyGrowth", nrow(r1$R_ci)), rep("EpiSky - Genetic only", length(1970:2009)), rep("EpiSky - Combination", length(1970:2009))))

# pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_skygrowth.pdf", width=8, height=4)
# library(ggplot2)
# ggplot(skygrowth_rt_df3) + geom_line(aes(x=Year, y=Mean, colour=Method)) +
#   geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Year, fill=Method), alpha=0.3) +
#   ylab("R(t)")
# dev.off()

height <- 2.5
width <- height*3
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_skygrowth.pdf", width=width, height=height)
par(mfrow=c(1,3))
par(mar=c(4.5,4.5,1.5,1.5))
plot(r1$time+2009, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="SkyGrowth"], ylim=c(0,5), xlim=c(1969,2010), type="l", main="SkyGrowth", xlab="Year", ylab="R(t)")
lines(r1$time+2009, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
lines(r1$time+2009, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
plot(1970:2009, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], ylim=c(0,5), xlim=c(1969,2010), type="l", main="EpiSky - Genetic only", xlab="Year", ylab="R(t)")
lines(1970:2009, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
lines(1970:2009, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
plot(1970:2009, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Combination"], ylim=c(0,5), xlim=c(1969,2010), type="l", main="EpiSky - Combination", xlab="Year", ylab="R(t)")
lines(1970:2009, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
lines(1970:2009, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
dev.off()

comparison_df <- data.frame("Year"=c(r1$time+2009, 1997:2009, 1996:2009, 1970:2009, 1970:2009),
                            "Mean"=c(apply(r1$R,2,mean), epiestimR$R$`Mean(R)`, apply(tb_g0e1$birth_rate[-(1:burn_in),], 2, mean)/death_rate, apply(tb_g1e0$birth_rate[-(1:burn_in),], 2, mean)/death_rate, apply(tb_g1e1$birth_rate[-(1:burn_in),], 2, mean)/death_rate),
                            "LCL"=c(r1$R_ci[,1], epiestimR$R$`Quantile.0.025(R)`, matrixStats::colQuantiles(tb_g0e1$birth_rate[-(1:burn_in),], probs=0.025)/death_rate, matrixStats::colQuantiles(tb_g1e0$birth_rate[-(1:burn_in),], probs=0.025)/death_rate, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], probs=0.025)/death_rate),
                            "UCL"=c(r1$R_ci[,3], epiestimR$R$`Quantile.0.975(R)`, matrixStats::colQuantiles(tb_g0e1$birth_rate[-(1:burn_in),], probs=0.975)/death_rate, matrixStats::colQuantiles(tb_g1e0$birth_rate[-(1:burn_in),], probs=0.975)/death_rate, matrixStats::colQuantiles(tb_g1e1$birth_rate[-(1:burn_in),], probs=0.975)/death_rate),
                            "Method"=c(rep("SkyGrowth", nrow(r1$R_ci)), rep("EpiEstim", length(1997:2009)), rep("EpiSky - ", length(1996:2009)), rep("d", length(1970:2009)), rep("e", length(1970:2009))))

# pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_comparison.pdf", width=8, height=4)
# library(ggplot2)
# ggplot(comparison_df) + geom_line(aes(x=Year, y=Mean, colour=Method)) +
#   geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Year, fill=Method), alpha=0.3) +
#   ylab("R(t)") +
#   scale_color_hue(labels=c("SkyGrowth", "EpiEstim", "EpiSky - Epidemic only", "EpiSky - Genetic only", "EpiSky - Combination")) +
#   scale_fill_hue(labels=c("SkyGrowth", "EpiEstim", "EpiSky - Epidemic only", "EpiSky - Genetic only", "EpiSky - Combination")) +
#   scale_x_continuous(breaks=seq(1970,2010,by=5), minor_breaks=seq(1970,2010,by=5))
# dev.off()

height <- 2.5
width <- height*5
pdf(file = "~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/tb_comparison.pdf", width=width, height=height)
par(mfrow=c(1,5))
par(mar=c(4.5,4.5,1.5,1.5))
plot(r1$time+2009, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="SkyGrowth"], ylim=c(0,30), xlim=c(1969,2010), type="l", main="SkyGrowth", xlab="Year", ylab="R(t)")
lines(r1$time+2009, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
lines(r1$time+2009, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="SkyGrowth"], lty=2)
plot(1997:2009, epiestim_df$Mean[epiestim_df$Method=="EpiEstim"], ylim=c(0,30), xlim=c(1969,2010), type="l", main="EpiEstim", xlab="Year", ylab="R(t)")
lines(1997:2009, epiestim_df$LCL[epiestim_df$Method=="EpiEstim"], lty=2)
lines(1997:2009, epiestim_df$UCL[epiestim_df$Method=="EpiEstim"], lty=2)
plot(1996:2009, epiestim_df$Mean[epiestim_df$Method=="EpiSky - Epidemic only"], ylim=c(0,30), xlim=c(1969,2010), type="l", main="EpiSky - Epidemic only", xlab="Year", ylab="R(t)")
lines(1996:2009, epiestim_df$LCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
lines(1996:2009, epiestim_df$UCL[epiestim_df$Method=="EpiSky - Epidemic only"], lty=2)
plot(1970:2009, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], ylim=c(0,30), xlim=c(1969,2010), type="l", main="EpiSky - Genetic only", xlab="Year", ylab="R(t)")
lines(1970:2009, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
lines(1970:2009, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Genetic only"], lty=2)
plot(1970:2009, skygrowth_rt_df3$Mean[skygrowth_rt_df3$Method=="EpiSky - Combination"], ylim=c(0,30), xlim=c(1969,2010), type="l", main="EpiSky - Combination", xlab="Year", ylab="R(t)")
lines(1970:2009, skygrowth_rt_df3$LCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
lines(1970:2009, skygrowth_rt_df3$UCL[skygrowth_rt_df3$Method=="EpiSky - Combination"], lty=2)
dev.off()

