task_id <- as.numeric(Sys.getenv("TASK_ID"))

library(EpiSky)
library(ape)

# ptree_raw <- read.tree( 'C:/Users/alici/OneDrive/Documents/_PhD/Year 4/Data/HIV/nc_refs_b.dater' )
# dates <- read.csv( 'C:/Users/alici/OneDrive/Documents/_PhD/Year 4/Data/HIV/nc_refs_b.dates.csv', stringsAs=FALSE )
# prev_paper <- read.csv( 'C:/Users/alici/OneDrive/Documents/_PhD/Year 4/Data/HIV/prevalence.csv', stringsAs=FALSE , skip = 4)
# prev_nida <- read.csv("~/_PhD/Year 4/Data/HIV/hiv_plot_data_rounded.csv", header=T)
ptree_raw <- read.tree( "/storage/stats/maundl/Year4/2_hiv/nc_refs_b.dater" )
dates <- read.csv( "/storage/stats/maundl/Year4/2_hiv/nc_refs_b.dates.csv", stringsAs=FALSE )
prev_paper <- read.csv( "/storage/stats/maundl/Year4/2_hiv/prevalence.csv", stringsAs=FALSE , skip = 4)

prev_paper <- prev_paper[,3:4]
prev_paper <- prev_paper[order(prev_paper$Year, decreasing = F),]
prev_paper$Cases <- as.numeric(gsub(",","",prev_paper$Cases))
prev_paper <- rbind(c(2009, 0), prev_paper)

sample100 <- sample_phylo(ptree_raw, ptree_lag=0, pi0=1, pi1=1)
ptree100 <- sample100$ptree
ptree_lag100 <- sample100$ptree_lag
# time <- floor(2019 - max(distToRoot(ptree100)))
time <- 1980

prev_paper_long <- rbind(data.frame("Year"=time:2008, "Cases"=rep(0,length(time:2008))), prev_paper)

prev_null <- prev_paper_long
prev_null[,2] <- 0


### ANALYSE DATA ###

iter <- 100000
max_time <- 47.5*60*60
sigma0 <- 0.05
reporting_prob0 <- 0.99
#x0 <- 3750
death_rate <- 1/10
x0_prior <- "nbinom"
n_particles <- NULL
min_n_particles <- 10000
max_n_particles <- 50000
resampling_scheme <- "systematic"
backward_sim <- TRUE
print <- TRUE

#epi only
if (task_id == 1) {
  x0 <- prev_paper[2,2]
  x0_mean <- 70000
  x0_var <- 800000000

  set.seed(1)
  hiv_g0e1_1 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = NULL, sample_prevalence = prev_paper,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g0e1_1.RData"))
}
if (task_id == 2) {
  x0 <- prev_paper[2,2]
  x0_mean <- 70000
  x0_var <- 800000000

  set.seed(2)
  hiv_g0e1_2 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = NULL, sample_prevalence = prev_paper,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g0e1_2.RData"))
}
if (task_id == 3) {
  x0 <- prev_paper[2,2]
  x0_mean <- 70000
  x0_var <- 800000000

  set.seed(3)
  hiv_g0e1_3 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = NULL, sample_prevalence = prev_paper,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g0e1_3.RData"))
}
if (task_id == 4) {
  x0 <- prev_paper[2,2]
  x0_mean <- 70000
  x0_var <- 800000000

  set.seed(4)
  hiv_g0e1_4 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = NULL, sample_prevalence = prev_paper,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g0e1_4.RData"))
}
if (task_id == 5) {
  x0 <- prev_paper[2,2]
  x0_mean <- 70000
  x0_var <- 800000000

  set.seed(5)
  hiv_g0e1_5 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = NULL, sample_prevalence = prev_paper,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g0e1_5.RData"))
}

#gen only
if (task_id == 6) {
  x0 <- 3750
  x0_mean <- 3230
  x0_var <- 28000

  set.seed(1)
  hiv_g1e0_1 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_null,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e0_1.RData"))
}
if (task_id == 7) {
  x0 <- 3750
  x0_mean <- 3230
  x0_var <- 28000

  set.seed(2)
  hiv_g1e0_2 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_null,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e0_2.RData"))
}
if (task_id == 8) {
  x0 <- 3750
  x0_mean <- 3230
  x0_var <- 28000

  set.seed(3)
  hiv_g1e0_3 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_null,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e0_3.RData"))
}
if (task_id == 9) {
  x0 <- 3750
  x0_mean <- 3230
  x0_var <- 28000

  set.seed(4)
  hiv_g1e0_4 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_null,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e0_4.RData"))
}
if (task_id == 10) {
  x0 <- 3750
  x0_mean <- 3230
  x0_var <- 28000

  set.seed(5)
  hiv_g1e0_5 <- pmmh(iter = iter, max_time = max_time,
                       sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                       death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_null,
                       x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                       n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                       resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e0_5.RData"))
}

# both long
if (task_id == 11) {
  x0 <- 3750
  x0_mean <- 3870
  x0_var <- 33000

  set.seed(1)
  hiv_g1e1_long_1 <- pmmh(iter = iter, max_time = max_time,
                            sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                            death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper_long,
                            x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                            n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                            resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_long_1.RData"))
}
if (task_id == 12) {
  x0 <- 3750
  x0_mean <- 3870
  x0_var <- 33000

  set.seed(2)
  hiv_g1e1_long_2 <- pmmh(iter = iter, max_time = max_time,
                            sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                            death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper_long,
                            x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                            n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                            resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_long_2.RData"))
}
if (task_id == 13) {
  x0 <- 3750
  x0_mean <- 3870
  x0_var <- 33000

  set.seed(3)
  hiv_g1e1_long_3 <- pmmh(iter = iter, max_time = max_time,
                            sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                            death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper_long,
                            x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                            n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                            resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_long_3.RData"))
}
if (task_id == 14) {
  x0 <- 3750
  x0_mean <- 3870
  x0_var <- 33000

  set.seed(4)
  hiv_g1e1_long_4 <- pmmh(iter = iter, max_time = max_time,
                            sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                            death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper_long,
                            x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                            n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                            resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_long_4.RData"))
}
if (task_id == 15) {
  x0 <- 3750
  x0_mean <- 3870
  x0_var <- 33000

  set.seed(5)
  hiv_g1e1_long_5 <- pmmh(iter = iter, max_time = max_time,
                            sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                            death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper_long,
                            x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                            n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                            resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_long_5.RData"))
}

# both short
if (task_id == 16) {
  x0 <- prev_paper[2,2]
  x0_mean <- 28000
  x0_var <- 200000

  set.seed(1)
  hiv_g1e1_short_1 <- pmmh(iter = iter, max_time = max_time,
                             sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                             death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper,
                             x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                             n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                             resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_short_1.RData"))
}
if (task_id == 17) {
  x0 <- prev_paper[2,2]
  x0_mean <- 28000
  x0_var <- 200000

  set.seed(2)
  hiv_g1e1_short_2 <- pmmh(iter = iter, max_time = max_time,
                             sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                             death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper,
                             x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                             n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                             resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_short_2.RData"))
}
if (task_id == 18) {
  x0 <- prev_paper[2,2]
  x0_mean <- 28000
  x0_var <- 200000

  set.seed(3)
  hiv_g1e1_short_3 <- pmmh(iter = iter, max_time = max_time,
                             sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                             death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper,
                             x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                             n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                             resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_short_3.RData"))
}
if (task_id == 19) {
  x0 <- prev_paper[2,2]
  x0_mean <- 28000
  x0_var <- 200000

  set.seed(4)
  hiv_g1e1_short_4 <- pmmh(iter = iter, max_time = max_time,
                             sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                             death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper,
                             x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                             n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                             resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_short_4.RData"))
}
if (task_id == 20) {
  x0 <- prev_paper[2,2]
  x0_mean <- 28000
  x0_var <- 200000

  set.seed(5)
  hiv_g1e1_short_5 <- pmmh(iter = iter, max_time = max_time,
                             sigma0 = sigma0, reporting_prob0 = reporting_prob0, x0 = x0,
                             death_rate = death_rate, ptree = ptree100, ptree_lag = ptree_lag100, sample_prevalence = prev_paper,
                             x0_prior = x0_prior, x0_mean = x0_mean, x0_var = x0_var,
                             n_particles = n_particles, min_n_particles = min_n_particles, max_n_particles = max_n_particles,
                             resampling_scheme = resampling_scheme, backward_sim = backward_sim, print = print)
  save.image(file=file.path("/storage/stats/maundl/Paper/RData", "hiv_g1e1_short_5.RData"))
}
