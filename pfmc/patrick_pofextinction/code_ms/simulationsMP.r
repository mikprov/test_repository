# simulations.R

# Overview
## Fig. 2: Frequency response of Spring Run Chinook salmon to white 
## noise at different levels of mean prespawning survival 
## (0.275 - just above collapse, 0.5, and 0.8).

## basic session info & logging

source("C:/Users/provo/Documents/GitHub/pfmc/patrick_pofextinction/code_ms/functions.R")
#install.packages(c("data.table","ggplot2","reshape2","fields","TSA","RSEIS","biwavelet","foreach","doParallel"))

# libraries
library(ggplot2)
library(fields)
library(TSA)
library(RSEIS)
library(biwavelet)
library(compiler)

# trying to fix problems with data.table and reshape2 conflicts
#unloadNamespace('reshape2')
unloadNamespace('data.table')
unloadNamespace('reshape2')
install.packages("data.table",type="source",dependencies = TRUE)

library(data.table)
library(reshape2)

# load foreach package for parallel processing
library(foreach)
# set up backend to do parallel processing
library(doParallel)
registerDoParallel() # defaults to half of the cores 

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

# parameters for subsequent noise generation and population modeling
# number of simulations for each parm combination
reps <- 1e3
# length of each simulation
len_sim <- 2^10
burn_in <- 1e2
n_freqs <- 2^9
phasein_len <- 0 # 2^4
N <- burn_in + len_sim + phasein_len

# Parm combination
# Frequency response of "population" at three levels of survival for white noise
alphaMult <- 4
EQsp <- 7500 #equilibrium level (beta in the BH that I've been using)
# Set mean survival and bounds/range 
meanPS <- as.character(c(0.275, 0.5, 0.8))
sigPSmult <- as.character(seq(0.1, 0.5, by = 0.1))

freqCont <- c("white", "p34", "pgt10", "p34gt10", "one_over_f")
surv1 <- 0.02 # first year ocean survival
surv2 <- surv3 <- 0.8 # ocean survival in later years
wanted_frac <- 0.5 # parameter to determine the fraction spawning at age-3 and 4 (equilibrium).

# determine the mean survival level based on the (low) survival rate
# the produces stock collapse [1/SPR == alpha] and then set survival
# rates at incrementally higher levels.

# Run simulations that approach the persistence threshold for each alpha level
# by dropping meanPS from some arbitrary distance above 1/SPR collapse 

# calculate the fraction spawning early
delta_e <- calc_de(wanted_frac, surv3)
# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS = 1)  
# multiplier of replacement line slope to set alpha 
alphaMult <- 4 
alpha <- alphaMult * 1/spr

sps_crash <- 1/{alpha*surv1*surv2*{delta_e + surv3*{1 - delta_e}}}

sps_mults <- seq(1.1, 1.5, by = 0.1) #seq(1.1, 1.5, by = 0.2)

# create white noise

# create "reps" number of white noise time series of length N using random sine wave approach
# with selected frequency contents
white_n <- customFR2ts(N = N, # number of time steps
                       reps = reps,
                       r_seed = r_seed,
                       amp = mk_white(freq = n_freqs),
                       freq = n_freqs) # mean = 0, variance/sd = 1

# Filtered/selected bandwidths

# Bandpass period 3-4 (frequencies 1/4 to 1/3) in white noise signals.

# rsin_34_n <- customFR2ts(N = N, # number of time steps
#                          reps = reps,
#                          r_seed = r_seed,
#                          amp = mk_rsin(N, highF=1/3, lowF=1/4)) 

rsin_34_n <- customFR2ts(N = N, # number of time steps
                         reps = reps,
                         r_seed = r_seed,
                         amp = mk_rsin2(N, highF=1/3, lowF=1/4, freq = n_freqs),
                         freq = n_freqs) 

rsin_gt10_n <- customFR2ts(N = N,
                           reps = reps,
                           r_seed = r_seed,
                           amp = mk_rsin2(N, lowF=0, highF=1/10, freq = n_freqs),
                           freq = n_freqs) 

red_beta_1 <- customFR2ts(N = N, # number of time steps
                          reps = reps,
                          r_seed = r_seed,
                          amp = mk_1_over_f_beta(N, beta = .5, freq = n_freqs),
                          freq = n_freqs) 

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

## Create data.table to store simulation results

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

# pre-compile functions to improve run time.
popSimPSvaryCmp <- cmpfun(popSimPSvary) 
#not exactly sure what is going on in cmpfun()
#this site might help if I wanted to read more: 
#https://www.r-statistics.com/2012/04/speed-up-your-r-code-using-a-just-in-time-jit-compiler/

# Ask Patrick about parSimCmp
parSimCmp <- function(dt, #a subsection of the storage df, split by mean pre-spawn surv
                      noise_list, #simulated noise sets (each a 1124x1000 matrix)
                      freq_cont = c("white", "p34", "pgt10", "p34gt10", "one_over_f"),#noise sets names in order as list
                      sim_len = len_sim, #yrs in one simulation w/o burn in time
                      burn_in = burn_in, #burn in time, not included in simulation length
                      phasein_len = phasein_len, #0, ignore for now
                      surv_mean, #one of the three levels
                      sd_surv, #could be one of 5 sdstevs
                      alpha_mult = 4,
                      EQ_sp = EQsp, #beta
                      frac_wanted = 0.5) { #function starts below...
  for (i in 1:length(freq_cont)) { #for each noise set (white, period 3-4, etc)
      for (k in sd_surv) { #step through each sdstev level...
        #create new pre-spawn survival matrix that has mean of 0.275, 0.5, or 0.8
        #for each level of sdstdev
        surv <- make_surv_mat(noise_dat = noise_list[[i]], 
                              mean_surv = as.numeric(surv_mean), 
                              sd_surv = as.numeric(k),
                              sim_len = sim_len,
                              phasein_len = phasein_len,
                              burn_in = burn_in) 
        for (l in alpha_mult) { #alpha_mult is always 4 in Patrick's analysis 
          for (m in EQ_sp) { #equilibrium is always 7500 in Patrick's analysis
            #subset parts of 
            dt[i = (sigPSmult_c == k & alphaMult_c == l & EQsp_c == m), 
               j = freq_cont[i] := list(melt(popSimPSvaryCmp(rand_surv = surv, 
                                                             surv1 = surv1, 
                                                             surv2 = surv2, 
                                                             surv3 = surv3, 
                                                             EQsp = m, 
                                                             wanted_frac = frac_wanted,
                                                             alpha_scale = l)[[2]])[,3])]
          }
        }
      } 
  }
  return(dt)
} # end of parSimCmp()

# run simulations using foreach framework to send jobs to multiple cores

system.time(
  storage <- foreach(h = 1:length(meanPS),.packages="reshape2") %dopar% {
    parSimCmp(dt = french[[h]], 
              noise_list = noiseList,
              freq_cont = c("white", "p34", "pgt10", "p34gt10", "one_over_f"),  
              sim_len = len_sim,
              burn_in = burn_in,
              phasein_len = phasein_len,
              surv_mean = as.numeric(meanPS[h]), 
              sd_surv = as.numeric(sigPSmult), 
              alpha_mult = 4,
              EQ_sp = EQsp,
              frac_wanted = 0.5)
  }
) 

storage <- rbindlist(storage)

rm(storageP, french, red_beta_1, rsin_34_gt10_n, rsin_34_n, rsin_gt10_n,
   white_n, alphaMult_r, EQsp_r, meanPS_r, n_r, reps_r, 
   sigPSmult_r)

# Fig 2: Plot frequency response at 3 survival levels 
setEPS()
# postscript(file.path(".", "output_ms", "Fig_2_white_noise_popFreqResp_with_TimeSeries_CV.ps"), width = 6.85, height = 6.85)

. <- "C:/Users/provo/Documents/GitHub/pfmc/patrick_pofextinction/code_ms"
tiff(file.path(".", "output_ms", "Fig2revised_100.tiff"), units="in",width = 6.85, height = 6.85, res=300)
#postscript(file.path(".", "output_ms", "Fig2revised_100.ps"), width = 6.85, height = 6.85)
# pdf(file.path(".", "output_ms", "Fig_2_white_noise_popFreqResp_with_TimeSeries_CV.pdf"), width = 6.85, height = 6.85)
old <- par(mar = c(4,5,1,1))
plot_dat <- storage[ i = N > (burn_in + phasein_len) & sigPSmult_c == "0.1"]
plotMeanFR_DTmany(plot_dat, N = len_sim, surv = as.numeric(meanPS[1]), scale = "CV", yaxis_lim = c(0,4))
linesMeanFR_DTmany(plot_dat, N = len_sim, surv = as.numeric(meanPS[2]), line_color = "black", scale = "CV")
linesMeanFR_DTmany(plot_dat, N = len_sim, surv = as.numeric(meanPS[3]), line_color = "grey30", scale = "CV")

legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)
par(old) 
dev.off()


## Fig 3 Summary plots of noise signals
plot_idx <- 2

# setEPS()
# postscript(file.path(".", "output_ms", "Fig_3_summaryFreqContTS_Noise.ps"), width = 6.85, height = 6.85)
tiff(file.path(".", "output_ms", "Fig3revised_100.tiff"),units="in", res=300, width = 6.85, height = 6.85)
#postscript(file.path(".", "output_ms", "Fig3revised_100.ps"), width = 6.85, height = 6.85)
# pdf(file.path(".", "output_ms", "Fig_3_summaryFreqContTS_Noise.pdf"), width = 6.85, height = 6.85)
plot_gen_freq_wvlt(noise = noiseList,
                   burn_in_pd = (burn_in + phasein_len),
                   num_rows2plt = 200,
                   n = plot_idx, 
                   J1 = trunc((log(32/(2 * 1))/log(2))/0.01))
dev.off()

## Fig 4. Summary plot of spawning female abundance

setEPS()
# postscript(file.path(".", "output_ms", "/Fig_4_summaryFreqContTS_SpawningFemales.ps"), width = 6.85, height = 6.85)
#postscript(file.path(".", "output_ms", "/Fig4revised_100.ps"), width = 8, height = 8)
tiff(file.path(".", "output_ms", "/Fig4revised_100.tiff"), units="in",res=300,width = 8, height = 8)
# pdf(file.path(".", "output_ms", "/Fig_4_summaryFreqContTS_SpawningFemales.pdf"), width = 8, height = 6)
plot_surv_spawn_ts(spawners = storage,
                   noise = noiseList,
                   burn_in_pd = burn_in,
                   phasein_len = phasein_len,
                   num_rows2plt = 200,
                   sim_len = len_sim,
                   meanSurv = "0.5",
                   sigPSmult = "0.2", 
                   n = plot_idx,
                   J1 = trunc((log(32/(2 * 1))/log(2))/0.01))
dev.off()

# Quasi-extinction metrics
QE_sim_len <- 100
upper_lim <- burn_in + phasein_len + QE_sim_len

storage_sub <- copy(storage[ i = N > (burn_in + phasein_len) & N <= upper_lim  ])

# Fig. 6 CALCULATE QE YEAR -- histogram of QE year ####

# A function wrapper format for data.tables
calcQEyr <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)]
}

qe_lev <- 100
sb_qeyr <- calcQEyr(dt = storage_sub, expr = list( as.integer(JA_consec(white, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(p34, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(pgt10, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(p34gt10, run_length = 4, qeLev = qe_lev)),
                                                   as.integer(JA_consec(one_over_f, run_length = 4, qeLev = qe_lev)))) 

# set names of the new columns
setnames(sb_qeyr, c("V1", "V2", "V3", "V4", "V5"),
         c("white", "p34", "pgt10", "p34gt10", "one_over_f"))

spectra_names <- c(
  `white` = "White",
  `p34` = "Cohort Frequencies",
  `pgt10` = "Low Frequencies",
  `p34gt10` = "Both",
  `one_over_f` = "1/f^0.5"
)

## Make figures
setEPS()
# postscript(file.path(".", "output_ms", "Fig_6_Surv_Freq_QE_time_Dist_lowSurv.ps"), width = 6.85, height = 8)
#postscript(file.path(".", "output_ms", "Fig6revised_100.ps"), width = 6.85, height = 8)
tiff(file.path(".", "output_ms", "Fig6revised_100.tiff"), units="in",res=300, width = 6.85, height = 8)
# pdf(file.path(".", "output_ms", "Fig_6_Surv_Freq_QE_time_Dist_lowSurv.pdf"), width = 6.85, height = 8)
# qet_tmp_sb_m <- as.data.table(melt(copy(sb_qeyr[ i = sigPSmult_c == "0.4" & N > 400 & meanPS_c == "0.3"]), id = c(1:5)))
qet_tmp_sb_m <- as.data.table(melt(copy(sb_qeyr[ i = sigPSmult_c == "0.4" & meanPS_c == "0.275"]), id = c(1:5)))
x <- ggplot(qet_tmp_sb_m, aes(x = value)) + 
  geom_histogram(bins = 25) + 
  facet_grid(variable ~ . , labeller = as_labeller(spectra_names)) +
  xlim(c(0, 100)) +
  xlab("Time (Year)") +
  ylab("Count") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) + 
  theme(strip.background = element_rect(fill="white")) 
print(x)

dev.off()


# Fig 5
# CALCULATE PQE ####

calc_pQE <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c)]
}

sb_pQE <- calc_pQE(sb_qeyr, list(white_qe = length(which(!is.na(white)))/reps,
                                 p34_qe = length(which(!is.na(p34)))/reps,
                                 pgt10_qe = length(which(!is.na(pgt10)))/reps,
                                 p34gt10_qe = length(which(!is.na(p34gt10)))/reps,
                                 one_over_f_qe = length(which(!is.na(one_over_f)))/reps))

spectra_names_qe <- c(
  `white_qe` = "White",
  `p34_qe` = "Cohort\nFrequencies",
  `pgt10_qe` = "Low Frequencies",
  `p34gt10_qe` = "Both",
  `one_over_f_qe` = "1/f^0.5"# "1/f"
)

#setEPS()
# postscript(file.path(".", "output_ms", "Fig_5_sigma_vs_pQE2row.ps"),
#            width = 6.85, height = 5.08)
#postscript(file.path(".", "output_ms", "Fig5revised_100.ps"), width = 6.85, height = 5.08)
tiff(file.path(".", "output_ms", "Fig5revised_100.tiff"), units="in",res=300,width = 6.85, height = 5.08)
# pdf(file.path(".", "output_ms", "Fig_5_sigma_vs_pQE2row.pdf"), width = 6.85, height = 5.08)
pqe_tmp_m <- melt(copy(sb_pQE[ i = alphaMult_c == alphaMult & meanPS_c %in% c(meanPS[1], meanPS[2])]), id = c(1:4))
x <- ggplot(pqe_tmp_m, aes(x = as.numeric(sigPSmult_c), y = value)) + 
  geom_line(colour = "gray30", size = 1.2) + 
  facet_grid(meanPS_c ~ variable, labeller = labeller(variable = spectra_names_qe)) +
  theme(axis.text = element_text(size = 8, colour = "black")) +
  ylab("Probability of Quasi-Extinction") +
  xlab(expression(paste("Standard Deviation of Environment (",  sigma, ", ", y^{-1}, ")" ))) +
  theme_bw()  +
  theme(strip.text.y = element_text(angle = 0)) + 
  theme(strip.background = element_rect(fill="white"))
print(x)

dev.off()


