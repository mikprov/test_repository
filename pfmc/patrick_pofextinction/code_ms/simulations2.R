# simulations.R

# Overview
## Fig. 5 & 6

## basic session info & logging

source("./code_ms/functions.R")

# libraries
library(ggplot2)
library(reshape2)
# library(plyr)
library(data.table)
library(TSA)
library(biwavelet)
# library(grid)
library(compiler)
# load foreach package for parallel processing
library(foreach)
# set up backend to do parallel processing
library(doParallel)
# detectCores()
registerDoParallel() # defaults to half of the cores 

popSimPSvaryCmp <- cmpfun(popSimPSvary)

parSimCmp <- function(dt, simLen = 1024, survMean) {
  for (i in 1:length(freqCont)) {
    for (k in sigPSmult) {
      surv <- matrix(NA, nrow = nrow(noiseList[[i]]), ncol = ncol(noiseList[[i]]))
      surv[1:(nrow(surv)-simLen), ] <- noiseList[[i]][1:(nrow(surv) - simLen), ] * 0.01 + survMean
      surv[(nrow(surv)-(simLen - 1)):nrow(surv), ] <- noiseList[[i]][(nrow(surv)-(simLen - 1)):nrow(surv), ] * k + survMean
      surv[surv > 1] <- 1
      surv[surv < 0] <- 0
      for (l in alphaMult) {
        for (m in EQsp) {
          dt[i = (sigPSmult_c == k & alphaMult_c == l & EQsp_c == m), 
             j = freqCont[i] := list(melt(popSimPSvaryCmp(rand_surv = surv, 
                                                          surv1 = surv1, 
                                                          surv2 = surv2, 
                                                          surv3 = surv3, 
                                                          EQsp = m, 
                                                          wanted_frac = wanted_frac,
                                                          alpha_scale = l)[[2]])[,3])]
        }
      }
    }
  }
  return(dt)
} # end of parSimCmp()



system("hostname")  # record name of the machine
date() # record the date
sessionInfo() # documents the version of R and any included packages, 
# to reproduce the environment that yielded the results

# set random seed
r_seed <- 64
set.seed(r_seed)

# prelims

# todays_date <-  Sys.Date()
diss_dir_name <- file.path(".", "output_ms", "sim_results", "dissertation_files")
if(!exists(diss_dir_name)) dir.create(diss_dir_name, recursive = TRUE)

support_dir_name <- file.path(".", "output_ms", "sim_results", "support_files")
if (!exists(support_dir_name)) dir.create(support_dir_name, recursive = TRUE)

##  params

# "fixed" parameters for subsequent noise generation and population modeling
reps <- 1000
N <- 1424
freqCont <- c("white", "p34", "pgt10", "p34gt10", "one_over_f")
surv1 <- 0.02 # first year ocean survival
surv2 <- surv3 <- 0.8 # ocean survival in later years
wanted_frac <- 0.5 # parameter to determine the fraction spawning at age-3 and 4 (equilibrium).

# Before generating the random variables of noise to force survival
# determine the mean survival level based on the (low) survival rate
# the produces stock collapse [1/SPR == alpha] and then set survival
# rates at incrementally higher levels.

# Run simulations that approach the persistence threshold for each alpha level
# by dropping meanPS from some arbitrary distance above 1/SPR collapse 

# calculate the fraction spawning early
delta_e <- calc_de(wanted_frac, surv3)
# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS = 1)  
alphaMult <- 4 # multiplier of replacement line slope to set alpha [from 2-10]
alpha <- alphaMult * 1/spr
sps_crash <- 1/{alpha*surv1*surv2*{delta_e + surv3*{1 - delta_e}}}

#sps_adds <- seq(0.05, 0.3, by = 0.05) # 0.05 increments in survival
sps_mults <- seq(1.1, 1.5, by = 0.1) #seq(1.1, 1.5, by = 0.2)


# create white noise

# create "reps" number of white noise time series of length N using random sine wave approach
# with selected frequency contents
white_n <- customFR2ts(N = N, # number of time steps
                       reps = reps,
                       r_seed = r_seed,
                       amp = mk_white(N)) # mean = 0, variance/sd = 1

# Filtered/selected bandwidths

# Bandpass period 3-4 (frequencies 1/4 to 1/3) in white noise signals.

rsin_34_n <- customFR2ts(N = N, # number of time steps
                         reps = reps,
                         r_seed = r_seed,
                         amp = mk_rsin(N, highF=1/3, lowF=1/4)) 

rsin_gt10_n <- customFR2ts(N = N,
                           reps = reps,
                           r_seed = r_seed,
                           amp = mk_rsin(N, lowF=0, highF=1/10) ) 

red_beta_1 <- customFR2ts(N = N, # number of time steps
                          reps = reps,
                          r_seed = r_seed,
                          amp = mk_1_over_f_beta(N, beta = .5)) 

# Band-pass greater than period 10 and period 3-4 (frequencies lower than 0.1 plus period 3-4)
# By adding the period 10 and greater noise to the period 3-4 only noise (divide by 2) 
rsin_34_gt10_n <- matrix(NA, nrow = N, ncol = reps)
for (i in 1:reps) {
  rsin_34_gt10_n[,i] <- (rsin_gt10_n[,i] + rsin_34_n[,i])/2
}
rsin_34_gt10_n <- apply(rsin_34_gt10_n, 2, scale)

# Store simulated noise sets in a list for use in simualtions (easier indexing)
noiseList <- list(noise_white = white_n,
                  noise_34 = rsin_34_n,
                  noise_gt10 = rsin_gt10_n,
                  noise_34gt10 = rsin_34_gt10_n,
                  noise_red_beta_1 = red_beta_1)

# rm individual sets of noise
# rm(white_n, rsin_34_n, rsin_gt10_n, rsin_34_gt10_n, red_beta_1, ar_noise_list)

## Parameterize for extinction cals

sigPSmult <-  seq(.1, 0.5, by = 0.1) # seq(.1, 0.5, by = 0.2)  seq(.1, 0.6, by = 0.025)
alphaMult <- 4
meanPS <- c(0.3, 0.325, 0.35, 0.40, 0.5)
EQsp <- 7500

meanPS_r <- rep(meanPS, each = length(sigPSmult)*length(alphaMult)*length(EQsp)*reps*N)
sigPSmult_r <- rep(sigPSmult, each = length(alphaMult)*length(EQsp)*reps*N, times = length(meanPS))
alphaMult_r <- rep(alphaMult, each = length(EQsp)*reps*N, times = length(meanPS)*length(sigPSmult))
EQsp_r <- rep(EQsp, each = reps*N, times = length(meanPS)*length(sigPSmult)*length(alphaMult))
reps_r <- rep(1:reps, each = N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp))
n_r <- rep(1:N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp)*reps)


storageP <- data.table::data.table(meanPS_c = meanPS_r,
                                   sigPSmult_c = sigPSmult_r,
                                   alphaMult_c = alphaMult_r,
                                   EQsp_c = EQsp_r,
                                   reps_c = reps_r,
                                   N = n_r,
                                   white = 0) 

data.table::setkey(storageP, meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)

# split storageP by meanPS so each component of the
# then I'll put all the lists back together at the end of the simulation
french <- vector("list", length = length(meanPS)) # random name, another version I used italian named in honor of LWB

for (i in 1:length(meanPS)) {
  french[[i]] <- storageP[i = meanPS_c == meanPS[i]]
}

# tmp <- PopProjPSvar(survPS = runif(100, 0.6, 0.8), surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale = 4)

system.time(
  storage <- foreach(h = 1:length(meanPS),.packages="reshape2") %dopar% {
    parSimCmp(french[[h]], simLen = 1024, survMean = meanPS[h])
  }
)

storage <- rbindlist(storage)
rm(storageP, french)

## Fig_5 (actually 6). Tabulate 
# A function wrapper format for data.tables

# qeLev <- 100
# run_length <- 4
# Call the function
storage_sub <- copy(storage[ i = N > 400])

# CALCULATE QE YEAR ####

# A function wrapper format for data.tables
calcQEyr <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)]
}

qe_lev <- 20
sb_qeyr <- calcQEyr(dt = storage_sub, expr = list( as.integer(JA_consec(white, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(p34, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(pgt10, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(p34gt10, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(one_over_f, run_length = 4, qeLev = qe_lev)))) 

# need to name the new columns
setnames(sb_qeyr, c("V1", "V2", "V3", "V4", "V5"),
         c("white", "p34", "pgt10", "p34gt10", "one_over_f"))

spectra_names <- c(
  `white` = "White",
  `p34` = "Cohort Frequencies",
  `pgt10` = "Low Frequencies",
  `p34gt10` = "Both",
  `one_over_f` = "1/f"
)

## Make figures
pdf(file.path(".", "output_ms", "Fig_6_Surv_Freq_QE_time_Dist_lowSurv.pdf"), width = 5, height = 8)
qet_tmp_sb_m <- as.data.table(melt(copy(sb_qeyr[ i = sigPSmult_c == 0.4 & N > 400 & meanPS_c == 0.30]), id = c(1:5)))
# qet_tmp_sb_m$variable <- as.character(qet_tmp_sb_m$variable)

x <- ggplot(qet_tmp_sb_m, aes(x = value)) + 
  geom_histogram() + 
  facet_grid(variable ~ . , labeller = as_labeller(spectra_names)) +
  xlim(c(0, 100)) +
  xlab("Time (Year)") +
  ylab("Count") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) + 
  theme(strip.background = element_rect(fill="white"))
print(x)

dev.off()


# Fig 6
# CALCULATE PQE ####

calc_pQE <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c)]
}

sb_pQE <- calc_pQE(sb_qeyr, list(white_qe = length(which(!is.na(white)))/1000,
                                 p34_qe = length(which(!is.na(p34)))/1000,
                                 pgt10_qe = length(which(!is.na(pgt10)))/1000,
                                 p34gt10_qe = length(which(!is.na(p34gt10)))/1000,
                                 one_over_f_qe = length(which(!is.na(one_over_f)))/1000))

spectra_names_qe <- c(
  `white_qe` = "White",
  `p34_qe` = "Cohort Frequencies",
  `pgt10_qe` = "Low Frequencies",
  `p34gt10_qe` = "Both",
  `one_over_f_qe` = "1/f"
)


pdf(file.path(".", "output_ms", "Fig_5_sigma_vs_pQE2row.pdf"), width = 12, height = 9)
pqe_tmp_m <- melt(copy(sb_pQE[ i = alphaMult_c == alphaMult & meanPS_c %in% c(meanPS[1], meanPS[5])]), id = c(1:4))
x <- ggplot(pqe_tmp_m, aes(x = sigPSmult_c, y = value)) + 
  geom_line(colour = "gray30", size = 1.2) + 
  facet_grid(meanPS_c ~ variable, labeller = labeller(variable = spectra_names_qe)) +
  theme(axis.text = element_text(size = 8, colour = "black")) +
  ylab("Probability of Quasi-Extinction") +
  xlab(expression(paste(sigma))) +
  theme_bw()  +
  theme(strip.text.y = element_text(angle = 0)) + 
  theme(strip.background = element_rect(fill="white"))
print(x)

dev.off()


