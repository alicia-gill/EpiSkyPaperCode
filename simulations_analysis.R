#Load data
files <- list.files(path="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/RData/", pattern="(sim).+(RData)")
for (i in 1:length(files)) {
  print(i)
  load(paste0("~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/RData/", files[i]))
}
gc()
devtools::load_all()
iter <- 100000

#Plot simulated data
library(ape)
width <- 10.5
height <- width*2/5
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sim_data.pdf", height=height, width=width)
par(mfcol=c(2,5))
par(mar=c(4.5,4.5,2.5,1.5))
for(i in 1:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  noisy_prev <- eval(parse(text=paste0("noisy",I)))[,2]
  plot(0:40, noisy_prev, type="l", xlab="Day", ylab="Observed prevalence", main=paste0(i,"% cases observed"), ylim=c(0,25))
  #  axis(side=2, at=seq(0,1000,by=200),labels=seq(0,1000,by=200))
  sample <- eval(parse(text=paste0("sample",I)))
  sample_tree <- sample$ptree
  sample_day <- sample$ptree_lag
  sample_tree$root_time <- stop_time - max(distToRoot(sample_tree)) - sample_day
  plot(sample_tree, main=paste0(i,"% sample - ", length(sample_tree$tip.label), " leaves"), show.tip.label=F)
  axisPhylo(backward=F, root.time=sample_tree$root_time)
}
#dev.off()

#Matrices of acceptance rates, number of particles, variance of log-likelihood and run time
burn_in <- iter/5
accept <- matrix(nrow=6,ncol=6,dimnames=list(paste0("g",perc),paste0("e",perc)))
nparticles <- matrix(nrow=6,ncol=6,dimnames=list(paste0("g",perc),paste0("e",perc)))
smcvar <- matrix(nrow=6,ncol=6,dimnames=list(paste0("g",perc),paste0("e",perc)))
runtime <- matrix(nrow=6,ncol=6,dimnames=list(paste0("g",perc),paste0("e",perc)))
npart_rtime <- data.frame("nparticles"=rep(NA,36), "runtime"=rep(NA,36), "epi"=rep(NA,36), "gen"=rep(NA,36))
for (i in 1:6) {
  I <- formatC(i-1, width=2, format="d", flag="0")
  for (j in 1:6) {
    J <- formatC(j-1, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    accept[i,j] <- round(100*current$acceptance_rate,3)
    nparticles[i,j] <- current$n_particles
    smcvar[i,j] <- round(var(current$smc_llik[-(1:burn_in)]),3)
    runtime[i,j] <- round(current$run_time/60,0)
    k <- (i-1)*6 + j
    npart_rtime[k,] <- c(current$n_particles, round(current$run_time/60,0), as.numeric(I), as.numeric(J))
  }
}
max(runtime)/60
accept
nparticles
smcvar
runtime
summary(lm(nparticles ~ epi + gen, npart_rtime))
lm <- lm(runtime ~ nparticles + epi + gen, npart_rtime)
# lm <- lm(runtime ~ nparticles, npart_rtime)
summary(lm)
summary(lm(runtime ~ nparticles + epi + gen, npart_rtime))

plot(NULL, xlim=c(0,600), ylim=c(0,200), xlab="Number of particles", ylab="Run time (mins)")
for (i in 1:6) {
  points(nparticles[,i], runtime[,i], col=i)
}
abline(lm(runtime~nparticles, npart_rtime))

#Trace plots of reporting probability
width <- 6*2
height <- width
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/pobs_trace.pdf", height=height, width=width)
par(mfcol=c(6,6))
par(mar=c(4,4,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(current$reporting_prob, type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xaxt="n", xlab="Iteration (10^3)", ylab="Reporting probability")
    axis(1, at=seq(0,iter,by=iter/4), labels=format(seq(0,iter/1000,by=iter/4000), scientific=F))
  }
}
#dev.off()

library(scales)
# min_pobs <- rep(0,6)
max_pobs <- rep(-Inf,6)
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    # min <- min(current$reporting_prob)
    max <- max(current$reporting_prob)
    # if (min < min_pobs[i+1]) {
    #   min_pobs[i+1] <- min
    # }
    if (max > max_pobs[i+1]) {
      max_pobs[i+1] <- max
    }
  }
}
max_pobs <- c(1,1,1,0.5,0.2,0.2)
height <- width <- 11.5
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/pobs_trace_test.pdf", height=height, width=width)
layout.matrix <- matrix(rbind(cbind(37:42,matrix(1:36,nrow=6,byrow=T)),43:49), nrow=7)
layout(mat=layout.matrix, widths=c(0.6,rep(2,6)), heights=c(rep(2,6),0.6))
par(mar=c(1,1,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(current$reporting_prob, xlim=c(0,100000), ylim=c(0,max_pobs[i+1]), type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xaxt="n", yaxt="n", xlab="", ylab="", col=alpha(colour="black",alpha=0.75))
    # plot(current$reporting_prob, xlim=c(0,100000), ylim=c(0,1), type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xlab="", ylab="")
  }
}
par(mar=c(1,4,2,0))
for (i in 1:6) {
  plot(NULL, xlim=c(0,0), ylim=c(0,max_pobs[i]), xaxt="n",frame=F,xlab="",ylab="Reporting probability")
}
plot(NULL,xaxt="n",yaxt="n",xlab="",ylab="",main="",frame=F,xlim=c(0,0),ylim=c(0,0))
par(mar=c(4,1,0,1))
for (i in 1:6) {
  plot(NULL, xlim=c(0,100000), ylim=c(0,0), xaxt="n",yaxt="n",frame=F,xlab="Iteration (10^3)",ylab="")
  axis(1, at=seq(0,iter,by=iter/4), labels=format(seq(0,iter/1000,by=iter/4000), scientific=F))
}
dev.off()

#Trace plots of smoothness
width <- 6*2
height <- width
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sigma_trace.pdf", height=height, width=width)
par(mfcol=c(6,6))
par(mar=c(4,4,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(current$sigma, type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xaxt="n", xlab="Iteration (10^3)", ylab="Sigma")
    axis(1, at=seq(0,iter,by=iter/4), labels=format(seq(0,iter/1000,by=iter/4000), scientific=F))
  }
}
#dev.off()

min_sigma <- rep(Inf,6)
max_sigma <- rep(-Inf,6)
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    min <- min(current$sigma)
    max <- max(current$sigma)
    if (min < min_sigma[i+1]) {
      min_sigma[i+1] <- min
    }
    if (max > max_sigma[i+1]) {
      max_sigma[i+1] <- max
    }
  }
}
max_sigma <- c(1.2,1,0.8,0.8,0.6,0.5)
height <- width <- 11.5
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sigma_trace_test.pdf", height=height, width=width)
layout.matrix <- matrix(rbind(cbind(37:42,matrix(1:36,nrow=6,byrow=T)),43:49), nrow=7)
layout(mat=layout.matrix, widths=c(0.6,rep(2,6)), heights=c(rep(2,6),0.6))
par(mar=c(1,1,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(current$sigma, xlim=c(0,100000), ylim=c(0,max_sigma[i+1]), type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xaxt="n", yaxt="n", xlab="", ylab="", col=alpha(colour="black",alpha=0.75))
    # plot(current$sigma, xlim=c(0,100000), ylim=c(0,1), type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xlab="", ylab="")
  }
}
par(mar=c(1,4,2,0))
for (i in 1:6) {
  plot(NULL, xlim=c(0,0), ylim=c(0,max_sigma[i]), xaxt="n",frame=F,xlab="",ylab="Sigma")
}
plot(NULL,xaxt="n",yaxt="n",xlab="",ylab="",main="",frame=F,xlim=c(0,0),ylim=c(0,0))
par(mar=c(4,1,0,1))
for (i in 1:6) {
  plot(NULL, xlim=c(0,100000), ylim=c(0,0), xaxt="n",yaxt="n",frame=F,xlab="Iteration (10^3)",ylab="")
  axis(1, at=seq(0,iter,by=iter/4), labels=format(seq(0,iter/1000,by=iter/4000), scientific=F))
}
dev.off()

#Trace plots of day 0 prevalence
width <- 6*2
height <- width
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/x0_trace.pdf", height=height, width=width)
par(mfcol=c(6,6))
par(mar=c(4,4,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(current$x0, type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xaxt="n", xlab="Iteration (10^3)", ylab="Day 0 prevalence")
    axis(1, at=seq(0,iter,by=iter/4), labels=format(seq(0,iter/1000,by=iter/4000), scientific=F))
  }
}
#dev.off()

min_x0 <- rep(Inf,6)
max_x0 <- rep(-Inf,6)
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    min <- min(current$x0)
    max <- max(current$x0)
    if (min < min_x0[i+1]) {
      min_x0[i+1] <- min
    }
    if (max > max_x0[i+1]) {
      max_x0[i+1] <- max
    }
  }
}
max_x0 <- c(100,60,30,40,30,30)
height <- width <- 11.5
pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/x0_trace_test.pdf", height=height, width=width)
layout.matrix <- matrix(rbind(cbind(37:42,matrix(1:36,nrow=6,byrow=T)),43:49), nrow=7)
layout(mat=layout.matrix, widths=c(0.6,rep(2,6)), heights=c(rep(2,6),0.6))
par(mar=c(1,1,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(current$x0, xlim=c(0,100000), ylim=c(0,max_x0[i+1]), type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xaxt="n", yaxt="n", xlab="", ylab="", col=alpha(colour="black",alpha=0.75))
    # plot(current$x0, xlim=c(0,100000), ylim=c(0,1), type="l", main=paste0("Gen ", i, "%, Epi ", j, "%"), xlab="", ylab="")
  }
}
par(mar=c(1,4,2,0))
for (i in 1:6) {
  plot(NULL, xlim=c(0,0), ylim=c(0,max_x0[i]), xaxt="n",frame=F,xlab="",ylab="Day 0 prevalence")
}
plot(NULL,xaxt="n",yaxt="n",xlab="",ylab="",main="",frame=F,xlim=c(0,0),ylim=c(0,0))
par(mar=c(4,1,0,1))
for (i in 1:6) {
  plot(NULL, xlim=c(0,100000), ylim=c(0,0), xaxt="n",yaxt="n",frame=F,xlab="Iteration (10^3)",ylab="")
  axis(1, at=seq(0,iter,by=iter/4), labels=format(seq(0,iter/1000,by=iter/4000), scientific=F))
}
dev.off()

#Reporting probability density
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  rho <- c(eval(parse(text=paste0("sim_g00e",I,"_2")))$reporting_prob[-(1:burn_in)],
           eval(parse(text=paste0("sim_g01e",I,"_2")))$reporting_prob[-(1:burn_in)],
           eval(parse(text=paste0("sim_g02e",I,"_2")))$reporting_prob[-(1:burn_in)],
           eval(parse(text=paste0("sim_g03e",I,"_2")))$reporting_prob[-(1:burn_in)],
           eval(parse(text=paste0("sim_g04e",I,"_2")))$reporting_prob[-(1:burn_in)],
           eval(parse(text=paste0("sim_g05e",I,"_2")))$reporting_prob[-(1:burn_in)])
  obsgen <- c(rep("0%",iter-burn_in),
              rep("1%",iter-burn_in),
              rep("2%",iter-burn_in),
              rep("3%",iter-burn_in),
              rep("4%",iter-burn_in),
              rep("5%",iter-burn_in))
  df <- data.frame("Reporting probability"=rho, "Observed genetic data"=obsgen)
  assign(paste0("pobs_e",I,"_df"), df)
  rm(df)
}
width <- 15
height <- width/5
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/pobs_density.pdf", height=height, width=width)
library(ggplot2)
library(ggridges)
library(gridExtra)
grid.arrange(
  ggplot(data=pobs_e01_df, aes(x=Reporting.probability, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.1)) +
    geom_vline(xintercept=0.01, linetype=2) +
    labs(x="Reporting probability", y="Proportion of genetic data", title="1% epidemic data") +
    scale_fill_grey(),
  ggplot(data=pobs_e02_df, aes(x=Reporting.probability, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.1)) +
    geom_vline(xintercept=0.02, linetype=2) +
    labs(x="Reporting probability", y="Proportion of genetic data", title="2% epidemic data") +
    scale_fill_grey(),
  ggplot(data=pobs_e03_df, aes(x=Reporting.probability, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.1)) +
    geom_vline(xintercept=0.03, linetype=2) +
    labs(x="Reporting probability", y="Proportion of genetic data", title="3% epidemic data") +
    scale_fill_grey(),
  ggplot(data=pobs_e04_df, aes(x=Reporting.probability, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.1)) +
    geom_vline(xintercept=0.04, linetype=2) +
    labs(x="Reporting probability", y="Proportion of genetic data", title="4% epidemic data") +
    scale_fill_grey(),
  ggplot(data=pobs_e05_df, aes(x=Reporting.probability, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.1)) +
    geom_vline(xintercept=0.05, linetype=2) +
    labs(x="Reporting probability", y="Proportion of genetic data", title="5% epidemic data") +
    scale_fill_grey(),
  nrow=1
)
#dev.off()
pobs_mean <- array(dim=c(6,6), dimnames=list(paste0("g",perc), paste0("e",perc)))
for (i in 1:6) {
  I <- formatC(i-1, width=2, format="d", flag="0")
  for (j in 1:6) {
    J <- formatC(j-1, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    pobs_mean[i,j] <- mean(current$reporting_prob[-(1:burn_in)])
  }
}
pobs_mean
apply(pobs_mean[-1,],2,mean)

#Smoothness density plots
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  sigma <- c(eval(parse(text=paste0("sim_g00e",I,"_2")))$sigma[-(1:burn_in)],
             eval(parse(text=paste0("sim_g01e",I,"_2")))$sigma[-(1:burn_in)],
             eval(parse(text=paste0("sim_g02e",I,"_2")))$sigma[-(1:burn_in)],
             eval(parse(text=paste0("sim_g03e",I,"_2")))$sigma[-(1:burn_in)],
             eval(parse(text=paste0("sim_g04e",I,"_2")))$sigma[-(1:burn_in)],
             eval(parse(text=paste0("sim_g05e",I,"_2")))$sigma[-(1:burn_in)])
  obsgen <- c(rep("0%",iter-burn_in),
              rep("1%",iter-burn_in),
              rep("2%",iter-burn_in),
              rep("3%",iter-burn_in),
              rep("4%",iter-burn_in),
              rep("5%",iter-burn_in))
  df <- data.frame("sigma"=sigma, "Observed genetic data"=obsgen)
  assign(paste0("sigma_e",I,"_df"), df)
  rm(df)
}
width <- 15
height <- width/5
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/sigma_density.pdf", height=height, width=width)
library(ggridges)
library(gridExtra)
grid.arrange(
  ggplot(data=sigma_e01_df, aes(x=sigma, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.3)) +
    #    geom_vline(xintercept=0.0577, linetype=2) +
    labs(x="Sigma", y="Proportion of genetic data", title="1% epidemic data") +
    scale_fill_grey(),
  ggplot(data=sigma_e02_df, aes(x=sigma, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.3)) +
    #    geom_vline(xintercept=0.0577, linetype=2) +
    labs(x="Sigma", y="Proportion of genetic data", title="2% epidemic data") +
    scale_fill_grey(),
  ggplot(data=sigma_e03_df, aes(x=sigma, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.3)) +
    #    geom_vline(xintercept=0.0577, linetype=2) +
    labs(x="Sigma", y="Proportion of genetic data", title="3% epidemic data") +
    scale_fill_grey(),
  ggplot(data=sigma_e04_df, aes(x=sigma, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.3)) +
    #    geom_vline(xintercept=0.0577, linetype=2) +
    labs(x="Sigma", y="Proportion of genetic data", title="4% epidemic data") +
    scale_fill_grey(),
  ggplot(data=sigma_e05_df, aes(x=sigma, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,0.3)) +
    #    geom_vline(xintercept=0.0577, linetype=2) +
    labs(x="Sigma", y="Proportion of genetic data", title="5% epidemic data") +
    scale_fill_grey(),
  nrow=1
)
#dev.off()

#Day 0 prevalence density plots
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  x0 <- c(eval(parse(text=paste0("sim_g00e",I,"_2")))$x0[-(1:burn_in)],
          eval(parse(text=paste0("sim_g01e",I,"_2")))$x0[-(1:burn_in)],
          eval(parse(text=paste0("sim_g02e",I,"_2")))$x0[-(1:burn_in)],
          eval(parse(text=paste0("sim_g03e",I,"_2")))$x0[-(1:burn_in)],
          eval(parse(text=paste0("sim_g04e",I,"_2")))$x0[-(1:burn_in)],
          eval(parse(text=paste0("sim_g05e",I,"_2")))$x0[-(1:burn_in)])
  obsgen <- c(rep("0%",iter-burn_in),
              rep("1%",iter-burn_in),
              rep("2%",iter-burn_in),
              rep("3%",iter-burn_in),
              rep("4%",iter-burn_in),
              rep("5%",iter-burn_in))
  df <- data.frame("x0"=x0, "Observed genetic data"=obsgen)
  assign(paste0("x0_e",I,"_df"), df)
  rm(df)
}
width <- 15
height <- width/5
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/x0_density.pdf", height=height, width=width)
library(ggridges)
library(gridExtra)
grid.arrange(
  ggplot(data=x0_e01_df, aes(x=x0, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(stat="binline", alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,25)) +
    geom_vline(xintercept=1, linetype=2) +
    labs(x="Day 0 prevalence", y="Proportion of genetic data", title="1% epidemic data") +
    scale_fill_grey(),
  ggplot(data=x0_e02_df, aes(x=x0, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(stat="binline", alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,25)) +
    geom_vline(xintercept=1, linetype=2) +
    labs(x="Day 0 prevalence", y="Proportion of genetic data", title="2% epidemic data") +
    scale_fill_grey(),
  ggplot(data=x0_e03_df, aes(x=x0, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(stat="binline", alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,25)) +
    geom_vline(xintercept=1, linetype=2) +
    labs(x="Day 0 prevalence", y="Proportion of genetic data", title="3% epidemic data") +
    scale_fill_grey(),
  ggplot(data=x0_e04_df, aes(x=x0, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(stat="binline", alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,25)) +
    geom_vline(xintercept=1, linetype=2) +
    labs(x="Day 0 prevalence", y="Proportion of genetic data", title="4% epidemic data") +
    scale_fill_grey(),
  ggplot(data=x0_e05_df, aes(x=x0, y=Observed.genetic.data, fill=Observed.genetic.data)) +
    geom_density_ridges(stat="binline", alpha=0.6) +
    theme_ridges() +
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(0,25)) +
    geom_vline(xintercept=1, linetype=2) +
    labs(x="Day 0 prevalence", y="Proportion of genetic data", title="5% epidemic data") +
    scale_fill_grey(),
  nrow=1
)
#dev.off()

for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  day <- rep(1:stop_time, 6)
  mean <- c(apply(eval(parse(text=paste0("sim_g00e",I,"_2")))$birth_rate[-(1:burn_in),], 2, mean),
            apply(eval(parse(text=paste0("sim_g01e",I,"_2")))$birth_rate[-(1:burn_in),], 2, mean),
            apply(eval(parse(text=paste0("sim_g02e",I,"_2")))$birth_rate[-(1:burn_in),], 2, mean),
            apply(eval(parse(text=paste0("sim_g03e",I,"_2")))$birth_rate[-(1:burn_in),], 2, mean),
            apply(eval(parse(text=paste0("sim_g04e",I,"_2")))$birth_rate[-(1:burn_in),], 2, mean),
            apply(eval(parse(text=paste0("sim_g05e",I,"_2")))$birth_rate[-(1:burn_in),], 2, mean))
  lcl <- c(matrixStats::colQuantiles(eval(parse(text=paste0("sim_g00e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.025),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g01e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.025),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g02e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.025),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g03e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.025),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g04e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.025),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g05e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.025))
  ucl <- c(matrixStats::colQuantiles(eval(parse(text=paste0("sim_g00e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.975),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g01e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.975),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g02e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.975),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g03e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.975),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g04e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.975),
           matrixStats::colQuantiles(eval(parse(text=paste0("sim_g05e",I,"_2")))$birth_rate[-(1:burn_in),], probs=0.975))
  obsgen <- c(rep("0%",stop_time),
              rep("1%",stop_time),
              rep("2%",stop_time),
              rep("3%",stop_time),
              rep("4%",stop_time),
              rep("5%",stop_time))
  df <- data.frame("Day"=day, "mean"=mean, "lcl"=lcl, "ucl"=ucl, "Observed genetic data"=obsgen)
  assign(paste0("bt_e",I,"_df"), df)
  rm(df)
}
gc()

p0 <- ggplot() +
  geom_ribbon(data=bt_e00_df, aes(x=Day, ymin=lcl, ymax=ucl, color=Observed.genetic.data, fill=Observed.genetic.data), alpha=0.4, linetype=2) +
  geom_line(data=bt_e00_df, aes(x=Day, y=mean, color=Observed.genetic.data), linewidth=1) +
  labs(y="Birth rate", title="0% epidemic data", fill="Observed genetic data", color="Observed genetic data") +
  geom_line(data=data.frame("Day"=1:40, "Truth"=c(seq(0.1,0.3,length.out=21)[-1],seq(0.3,0.1,length.out=21)[-1])), aes(x=Day, y=Truth), linetype=4, linewidth=1) +
  scale_fill_grey(start=0.8, end=0.2) +
  scale_color_grey(start=0.8, end=0.2)
p1 <- ggplot() +
  geom_ribbon(data=bt_e01_df, aes(x=Day, ymin=lcl, ymax=ucl, color=Observed.genetic.data, fill=Observed.genetic.data), alpha=0.4, linetype=2) +
  geom_line(data=bt_e01_df, aes(x=Day, y=mean, color=Observed.genetic.data), linewidth=1) +
  labs(y="Birth rate", title="1% epidemic data", fill="Observed genetic data", color="Observed genetic data") +
  geom_line(data=data.frame("Day"=1:40, "Truth"=c(seq(0.1,0.3,length.out=21)[-1],seq(0.3,0.1,length.out=21)[-1])), aes(x=Day, y=Truth), linetype=4, linewidth=1) +
  scale_fill_grey(start=0.8, end=0.2) +
  scale_color_grey(start=0.8, end=0.2) +
  ylim(c(0,0.8))
p2 <- ggplot() +
  geom_ribbon(data=bt_e02_df, aes(x=Day, ymin=lcl, ymax=ucl, color=Observed.genetic.data, fill=Observed.genetic.data), alpha=0.4, linetype=2) +
  geom_line(data=bt_e02_df, aes(x=Day, y=mean, color=Observed.genetic.data), linewidth=1) +
  labs(y="Birth rate", title="2% epidemic data", fill="Observed genetic data", color="Observed genetic data") +
  geom_line(data=data.frame("Day"=1:40, "Truth"=c(seq(0.1,0.3,length.out=21)[-1],seq(0.3,0.1,length.out=21)[-1])), aes(x=Day, y=Truth), linetype=4, linewidth=1) +
  scale_fill_grey(start=0.8, end=0.2) +
  scale_color_grey(start=0.8, end=0.2) +
  ylim(c(0,0.8))
p3 <- ggplot() +
  geom_ribbon(data=bt_e03_df, aes(x=Day, ymin=lcl, ymax=ucl, color=Observed.genetic.data, fill=Observed.genetic.data), alpha=0.4, linetype=2) +
  geom_line(data=bt_e03_df, aes(x=Day, y=mean, color=Observed.genetic.data), linewidth=1) +
  labs(y="Birth rate", title="3% epidemic data", fill="Observed genetic data", color="Observed genetic data") +
  geom_line(data=data.frame("Day"=1:40, "Truth"=c(seq(0.1,0.3,length.out=21)[-1],seq(0.3,0.1,length.out=21)[-1])), aes(x=Day, y=Truth), linetype=4, linewidth=1) +
  scale_fill_grey(start=0.8, end=0.2) +
  scale_color_grey(start=0.8, end=0.2) +
  ylim(c(0,0.8))
p4 <- ggplot() +
  geom_ribbon(data=bt_e04_df, aes(x=Day, ymin=lcl, ymax=ucl, color=Observed.genetic.data, fill=Observed.genetic.data), alpha=0.4, linetype=2) +
  geom_line(data=bt_e04_df, aes(x=Day, y=mean, color=Observed.genetic.data), linewidth=1) +
  labs(y="Birth rate", title="4% epidemic data", fill="Observed genetic data", color="Observed genetic data") +
  geom_line(data=data.frame("Day"=1:40, "Truth"=c(seq(0.1,0.3,length.out=21)[-1],seq(0.3,0.1,length.out=21)[-1])), aes(x=Day, y=Truth), linetype=4, linewidth=1) +
  scale_fill_grey(start=0.8, end=0.2) +
  scale_color_grey(start=0.8, end=0.2) +
  ylim(c(0,0.8))
p5 <- ggplot() +
  geom_ribbon(data=bt_e05_df, aes(x=Day, ymin=lcl, ymax=ucl, color=Observed.genetic.data, fill=Observed.genetic.data), alpha=0.4, linetype=2) +
  geom_line(data=bt_e05_df, aes(x=Day, y=mean, color=Observed.genetic.data), linewidth=1) +
  labs(y="Birth rate", title="5% epidemic data", fill="Observed genetic data", color="Observed genetic data") +
  geom_line(data=data.frame("Day"=1:40, "Truth"=c(seq(0.1,0.3,length.out=21)[-1],seq(0.3,0.1,length.out=21)[-1])), aes(x=Day, y=Truth), linetype=4, linewidth=1) +
  scale_fill_grey(start=0.8, end=0.2) +
  scale_color_grey(start=0.8, end=0.2) +
  ylim(c(0,0.8))
height <- 3
width <- height*6
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/bt.pdf", height=height, width=width)
library(gridExtra)
library(lemon)
grid_arrange_shared_legend(p0,p1,p2,p3,p4,p5,nrow=1,position="right")
#dev.off()

bt_array <- array(dim=c(6,6,40,3), dimnames=list(paste0("g",perc),paste0("e",perc),paste0("day",1:40),c("mean","lcl","ucl")))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    bt_array[i+1,j+1,,1] <- apply(current$birth_rate, 2, mean)
    bt_array[i+1,j+1,,2] <- matrixStats::colQuantiles(current$birth_rate, probs=0.025)
    bt_array[i+1,j+1,,3] <- matrixStats::colQuantiles(current$birth_rate, probs=0.975)
  }
}
truth <- c(seq(0.1,0.3,length.out=21)[-1],seq(0.3,0.1,length.out=21)[-1])
height <- 6*2
width <- height
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/bt_grid.pdf", height=height, width=width)
par(mfrow=c(6,6))
par(mar=c(4,4,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(i, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(1:40, bt_array[i+1,j+1,,1], type="l", ylim=c(0,0.8), xlab="Day", ylab="Birth rate", main=paste0("Gen ", i, "%, Epi ", j, "%"))
    lines(1:40, bt_array[i+1,j+1,,2], lty=2)
    lines(1:40, bt_array[i+1,j+1,,3], lty=2)
    lines(1:40, truth, col="red")
  }
}
#dev.off()

height <- width <- 11.5
#pdf(file="~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/Plots/bt_grid_test.pdf", height=height, width=width)
layout.matrix <- matrix(rbind(cbind(37:42,matrix(1:36,nrow=6,byrow=T)),43:49), nrow=7)
layout(mat=layout.matrix, widths=c(0.6,rep(2,6)), heights=c(rep(2,6),0.6))
par(mar=c(1,1,2,1))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(j, width=2, format="d", flag="0")
    current <- eval(parse(text=paste0("sim_g",I,"e",J,"_2")))
    plot(1:40, bt_array[i+1,j+1,,1], type="l", ylim=c(0,0.8), main=paste0("Gen ", i, "%, Epi ", j, "%"), xaxt="n", yaxt="n", xlab="", ylab="")
    lines(1:40, bt_array[i+1,j+1,,2], lty=2)
    lines(1:40, bt_array[i+1,j+1,,3], lty=2)
    lines(1:40, truth, col="red")
  }
}
par(mar=c(1,4,2,0))
for (i in 1:6) {
  plot(NULL, xlim=c(0,0), ylim=c(0,0.8), xaxt="n",frame=F,xlab="",ylab="Birth rate")
}
plot(NULL,xaxt="n",yaxt="n",xlab="",ylab="",main="",frame=F,xlim=c(0,0),ylim=c(0,0))
par(mar=c(4,1,0,1))
for (i in 1:6) {
  plot(NULL, xlim=c(1,40), ylim=c(0,0), yaxt="n",frame=F, xlab="Day",ylab="")
}
#dev.off()


#rmse
rmse <- function(sample, truth) {
  return(sqrt(mean((sample-truth)^2)))
}
rmse_array <- array(dim=c(6,6), dimnames=list(paste0("g",perc),paste0("e",perc)))
coverage_array <- array(dim=c(6,6), dimnames=list(paste0("g",perc),paste0("e",perc)))
ci_width_array <- array(dim=c(6,6), dimnames=list(paste0("g",perc),paste0("e",perc)))
for (i in 0:5) {
  I <- formatC(i, width=2, format="d", flag="0")
  for (j in 0:5) {
    J <- formatC(i, width=2, format="d", flag="0")
    rmse_array[i+1,j+1] <- rmse(sample=bt_array[i+1,j+1,,1], truth=truth)
    coverage_array[i+1,j+1] <- mean(bt_array[i+1,j+1,,2] <= truth & bt_array[i+1,j+1,,3] >= truth)
    ci_width_array[i+1,j+1] <- mean(bt_array[i+1,j+1,,3] - bt_array[i+1,j+1,,2])
  }
}
rmse_array
coverage_array
ci_width_array

apply(rmse_array[,-1],1,mean)
apply(rmse_array[,-1],2,mean)
apply(coverage_array[,-1],1,mean)
apply(coverage_array[,-1],2,mean)
apply(ci_width_array[,-1],1,mean)
apply(ci_width_array[,-1],2,mean)

plot(0:5, apply(rmse_array,1,mean), type="l", xlab="Percentage of data", ylab="RMSE", ylim=c(0,0.3))
lines(0:5, apply(rmse_array,2,mean), col="red")
legend("topright", legend=c("epidemic level", "genetic level"), col=c("black", "red"), lty=1)
