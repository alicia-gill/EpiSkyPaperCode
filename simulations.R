library(EpiSky)

perc <- c("00", "01", "02", "03", "04", "05")
perc_num <- as.numeric(perc)

stop_time <- 40
death_rate <- 0.1

peaked <- function(t) {
  stop_time <- 40
  peak_time <- 20
  min <- 0.1
  max <- 0.3
  
  if (t < peak_time) {
    slope <- (max-min)/(peak_time)
    return(min + slope*t)
  } else {
    slope <- (min-max)/(stop_time-peak_time)
    return(max + slope*(t-peak_time))
  }
}

set.seed(1)
epi <- epi_fn(birth_rate_fn = peaked, death_rate = death_rate, stop_time = stop_time, x0 = 1)
prev <- prevalence(epidemic = epi, stop_time = stop_time)

ptree <- phylo_tree(epi, stop_time)
set.seed(2)
sample05 <- sample_phylo(ptree = ptree, ptree_lag = 0, pi0 = 0, pi1 = 0.05)
set.seed(2)
sample04 <- sample_phylo(ptree = sample05$ptree, ptree_lag = sample05$ptree_lag, pi0 = 0, pi1 = 4/5)
set.seed(2)
sample03 <- sample_phylo(ptree = sample04$ptree, ptree_lag = sample04$ptree_lag, pi0 = 0, pi1 = 3/4)
set.seed(2)
sample02 <- sample_phylo(ptree = sample03$ptree, ptree_lag = sample03$ptree_lag, pi0 = 0, pi1 = 2/3)
set.seed(2)
sample01 <- sample_phylo(ptree = sample02$ptree, ptree_lag = sample02$ptree_lag, pi0 = 0, pi1 = 1/2)
sample00 <- sample_phylo(ptree = ptree, pi0 = 0, pi1 = 0)

set.seed(3)
noisy05 <- sample_prevalence(prev, 0.05)
set.seed(3)
noisy04 <- sample_prevalence(noisy05, 4/5)
set.seed(3)
noisy03 <- sample_prevalence(noisy04, 3/4)
set.seed(3)
noisy02 <- sample_prevalence(noisy03, 2/3)
set.seed(3)
noisy01 <- sample_prevalence(noisy02, 1/2)
noisy00 <- sample_prevalence(prev, 0)

library(doParallel)
registerDoParallel(cores=12)

iter <- 100000
max_time <- 8*60*60/3
sigma0 <- 0.05
reporting_prob0 <- 0.03
x0 <- 1
x0_prior <- "nbinom"
x0_mean <- 5
x0_var <- 50
n_particles <- NULL
min_n_particles <- 1
max_n_particles <- 1000
resampling_scheme <- "systematic"
backward_sim <- TRUE
print <- TRUE

chains <- foreach(task_id=1:36) %dopar% {
  i <- floor((task_id-1)/6)
  I <- perc[i+1]
  j <- ifelse(task_id %% 6 == 0, 6, task_id %% 6)
  J <- perc[j]
  
  set.seed(4)
  pmmh(iter = iter, max_time = max_time, sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0, death_rate = death_rate, ptree = eval(parse(text=paste0("sample",I,"$ptree"))), ptree_lag = eval(parse(text=paste0("sample",I,"$ptree_lag"))), sample_prevalence = eval(parse(text=paste0("noisy", J, ""))), x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var, n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles, resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
}

for (task_id in 1:36) {
  i <- floor((task_id-1)/6)
  I <- perc[i+1]
  j <- ifelse(task_id %% 6 == 0, 6, task_id %% 6)
  J <- perc[j]
  
  assign(paste0("sim_g", I, "e", J, "_2"), chains[[task_id]])
  save.image(file=file.path("~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/JRSSC Paper/RData", paste0("sim_g", I, "e", J, "_2.RData")))
}
