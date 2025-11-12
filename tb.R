task_id <- as.numeric(Sys.getenv("TASK_ID"))

library(EpiSky)
library(ape)

# ptree_raw <- read.tree("~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/Data/TB/eldholm2016.nwk")
# prev_raw <- read.csv("~/Documents/Papers/My papers/Bayesian Inference of Reproduction Number from Epidemiological and Genetic Data Using Particle MCMC/Data/TB/eldholm2016.csv", header=T)
# ptree_raw <- read.tree("C:/Users/alici/OneDrive/Documents/_PhD/Year 4/Data/TB/eldholm2016.nwk")
# prev_raw <- read.csv("C:/Users/alici/OneDrive/Documents/_PhD/Year 4/Data/TB/eldholm2016.csv", header=T)
ptree_raw <- read.tree("/storage/stats/maundl/Year4/6_tb/eldholm2016.nwk")
prev_raw <- read.csv("/storage/stats/maundl/Year4/6_tb/eldholm2016.csv", header=T)

#plot(ptree_raw, show.tip.label=F)
#axisPhylo()
#plot(prev_raw, type="o")

ptree <- sample_phylo(ptree_raw, pi0=1, pi1=1)$ptree
ptree$root_time <- 2009 - max(distToRoot(ptree))
#plot(ptree, show.tip.label=F)
#ape::axisPhylo()

prev_short <- rbind(data.frame("Year"=1995, "Incidence"=0), prev_raw)
prev_raw_long <- rbind(data.frame("Year"=1969:1995, "Incidence"=rep(0,27)), prev_raw)
#prev_raw_long

prev_null <- prev_raw_long
prev_null[,2] <- 0

iter <- 100000
max_time <- 47.5*60*60
sigma0 <- 0.1
reporting_prob0 <- 0.05
death_rate <- 1/3
x0_prior <- "nbinom"
n_particles <- NULL
min_n_particles <- 10000
max_n_particles <- 50000
resampling_scheme <- "systematic"
backward_sim <- TRUE
print <- TRUE

## Epi only
if (task_id %in% 1:5) {
  i <- ifelse(task_id %% 5 == 0, 5, task_id %% 5)
  x0 <- 200
  x0_mean <- 10000
  x0_var <- 380000000

  set.seed(i)
  chain <- pmmh(iter = iter, max_time = max_time,
                  sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                  death_rate = death_rate, ptree = NULL, sample_prevalence = prev_short,
                  x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                  n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                  resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)

  assign(paste0("tb_g0e1_",i), chain)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", paste0("tb_g0e1_",i,".RData")))
}

## Gen only
if (task_id %in% 6:10) {
  i <- ifelse(task_id %% 5 == 0, 5, task_id %% 5)
  x0 <- 5
  x0_mean <- 5
  x0_var <- 20

  set.seed(i)
  chain <- pmmh(iter = iter, max_time = max_time,
                  sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                  death_rate = death_rate, ptree = ptree, ptree_lag = 0, sample_prevalence = prev_null,
                  x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                  n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                  resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)

  assign(paste0("tb_g1e0_",i), chain)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", paste0("tb_g1e0_",i,".RData")))
}

## Epi+Gen long
if (task_id %in% 11:15) {
  i <- ifelse(task_id %% 5 == 0, 5, task_id %% 5)
  x0 <- 5
  x0_mean <- 7
  x0_var <- 56

  set.seed(i)
  chain <- pmmh(iter = iter, max_time = max_time,
                  sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                  death_rate = death_rate, ptree = ptree, ptree_lag = 0, sample_prevalence = prev_raw_long,
                  x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                  n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                  resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)

  assign(paste0("tb_g1e1_long_",i), chain)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", paste0("tb_g1e1_long_",i,".RData")))
}

## Epi+Gen short
if (task_id %in% 16:20) {
  i <- ifelse(task_id %% 5 == 0, 5, task_id %% 5)
  x0 <- 200
  x0_mean <- 210
  x0_var <- 211

  set.seed(i)
  chain <- pmmh(iter = iter, max_time = max_time,
                  sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                  death_rate = death_rate, ptree = ptree, ptree_lag = 0, sample_prevalence = prev_short,
                  x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                  n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                  resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)

  assign(paste0("tb_g1e1_short_",i), chain)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", paste0("tb_g1e1_short_",i,".RData")))
}
