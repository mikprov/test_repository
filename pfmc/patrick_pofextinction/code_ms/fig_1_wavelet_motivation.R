# Figure 1. Motivation: Summer survival rates of Butte Creek spring-run Chinook 
# salmon are projected to decline as summer stream temperatures increase and flows 
# decrease from climate change. (a) Example time series of summer adult survival 
# rates from Thompson et al. (2012). (b) Periods of environmental variability equal 
# to the mean generation time (approx. 3-4 years) of these salmon are evident in 
# the wavelet power spectrum from the middle 2060s to the late 2080s. The white 
# dashed line denotes the cone of influence that indicates the region where which 
# edge effects distort the spectrum. Thick black contour lines denote regions 
# where variance is significantly greater (α = 0.05) than a red-noise process with 
# the same lag-1 autocorrelation as the original time series (Torrence and Compo 1998).

## basic session info & logging

system("hostname")  # record name of the machine
date() # record the date
sessionInfo() # documents the version of R and any included packages, 
# to reproduce the environment that yielded the results

# set random seed
r_seed <- 64
set.seed(r_seed)

## libraries
library(biwavelet)
library(RSEIS)

## data

load_data_SALMOD_series <- function(management = c("BAU", "NoDiversion", "ColdWater", "ForeCast", "ForecastError", 
                                                   "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade"),
                                    emissions = c("A2","B1"),
                                    global_mods = c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")) {
  
  return()
}

## parameters



## wrangle

get_simulation_data_wavelet_plots <- function(management, # c("BAU", "NoDiversion")
                                              emissions, # c("A2","B1")
                                              global_mods) # c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1") 
{
  # get Thompson et al. prespawn survival SALMOD simulation data to generate
  # the motivating wavelet plot
  # params to generate this plot
  # management = "NoDiversion"
  # emissions = "B1"
  # global_mods = "mpiecham5"
  sim_path <- file.path(".", "data_ms", "simulation_results")
  
  mgmtScenarios <- management 
  climateScenarios <- emissions
  
  salmod_out <- read.delim(file.path(sim_path, mgmtScenarios, 
                                     paste0(climateScenarios, "_", global_mods), 
                                     "SALMODsumOutMerge.txt"), sep = ",")
  
  mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
  mortFW$global_mods <- global_mods
  years <- 2010:2099
  prespSurv <-  with(mortFW, data.frame(year = years, gcms = global_mods, PreSpawn = allMortSF/AFem))
  eggSurv <- with(mortFW, data.frame(year = years, gcms = global_mods, Egg = FryGrad/Eggs))
  frySurv <- with(mortFW, data.frame(year = years, gcms = global_mods, Fry = FryExit/FryGrad))
  eggSurv$Egg[which(is.na(eggSurv$Egg))] <- 0
  frySurv$Fry[which(is.na(frySurv$Fry))] <- 0
  
  return(list(prespSurv, eggSurv, frySurv))
}

wavelet_plot <- function(df, 
                         management,
                         emissions, 
                         global_mods) {
  mgmtScenarios <- management 
  climateScenarios <- emissions
  gcms <- global_mods 
  gcmScenarios <- numeric(length(climateScenarios)*length(global_mods))
  
  # number of scales minus - 1
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) 
  
  # Calculate and visualize wavelet spectra
  
  out_wavelet <- file.path(".", "output_ms")
  if (!dir.exists(out_wavelet)) { dir.create(out_wavelet, recursive = TRUE) }
  
  # path4file <- file.path(out_wavelet, paste0("Fig_1_wavelet_", mgmtScenarios, "_", climateScenarios, "_", global_mods, ".pdf"))
  # pdf(file = path4file, width = 5.07, height = 5.07)

  # path4file <- file.path(out_wavelet, paste0("Fig_1_wavelet_", mgmtScenarios, "_", climateScenarios, "_", global_mods, ".ps"))
  path4file <- file.path(out_wavelet, paste0("Fig1.ps"))
  setEPS()
  postscript(file = path4file, width = 5.07, height = 5.07)

  # Time series
  par(fig = c(0, 1, 0.5, 1))
  old <- par(mar = c(2, 4, 1, 1) )
  plot(df$year, df$PreSpawn, type = "l", lwd = 2, col = "darkgrey",
       axes = FALSE, ylab = "Survival Rate", 
       xlim = c(2010, 2100), xaxs = "i")
  axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
       lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
       las = 1, cex.axis = 1, tck = 0.02)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       las = 2, cex.axis = 1, tck = 0.02)
  box(lwd = 2)
  # title(paste0("GCM: ", gcms, " Prespawn survival emissions scenario: ", 
  #              climateScenarios, " \n Water managament alternative: ", 
  #              mgmtScenarios))
  par(old)
  
  # Wavelet Power Spectrum
  par(fig= c(0, 1, 0, 0.5), new = TRUE)
  old <- par(mar = c(4, 4, 1, 1) )
  test_cwt <- wt(cbind(df$year, detrend(as.numeric(scale(df$PreSpawn)))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100), 
       xlab = "Year", ylab = "Period", lty.coi = 1, lwd.coi = 2)
  axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
       lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
       las = 1, cex.axis = 1, tck = 0.02)
  par(old)
  
  dev.off()
}

dat <- get_simulation_data_wavelet_plots(management = "NoDiversion",
                                         emissions = "B1", 
                                         global_mods = "mpiecham5")
ps_dat <- dat[[1]]

wavelet_plot(df = ps_dat,
             management = "NoDiversion",
             emissions = "B1",
             global_mods = "mpiecham5")



# wavelet_SALMOD_series <- function(management = "NoDiversion",
#                                     emissions = c("A2","B1"),
#                                     global_mods = c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1") 
#   ) {
#     
#     # Function that extracts a specific survival time series (prespawn, egg and fry) 
#     # from SALMOD output files for each climate and selected water management scenario.
#     # Then it prints a pdf with plots of the ts, WPS, global WPS and scale-averaged WPS
#     
#     sim_path <- file.path(".", "data", "simulation_results")
#     
#     # Suite of mgmtScenarios: 
#     #"ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
#     #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
#     mgmtScenarios <- management 
#     climateScenarios <- emissions
#     gcms <- global_mods 
#     gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
#     
#     J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1; I translated from your matlab code
#     
#     for (m in 1:length(mgmtScenarios))  {
#       for (c in 1:length(climateScenarios)) {
#         
#         for (k in 1:length(gcms)) {
#           salmod_out <- read.delim(file.path(sim_path, mgmtScenarios[m], 
#                                              paste0(climateScenarios[c],"_", gcms[k]), "SALMODsumOutMerge.txt"), 
#                                    sep = ",")
#           mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
#           mortFW$gcms <- gcms[k]
#           
#           if(k == 1) mortFWgcms <- mortFW else mortFWgcms <- rbind(mortFWgcms, mortFW)
#         } # end of k loop
#         
#         prespSurv <-  with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, PreSpawn = allMortSF/AFem))
#         eggSurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Egg = FryGrad/Eggs))
#         frySurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Fry = FryExit/FryGrad))
#         
#         prespSurvC <- dcast(prespSurv, year ~ gcms, value.var = "PreSpawn")
#         eggSurv$Egg[which(is.na(eggSurv$Egg))] <- 0
#         eggSurvC <- dcast(eggSurv, year ~ gcms, value.var = "Egg")
#         frySurv$Fry[which(is.na(frySurv$Fry))] <- 0
#         frySurvC <- dcast(frySurv, year ~ gcms, value.var = "Fry")
#         
#         # For each cliamte scenario and water management alternative, plot the
#         # the time series and wavelet spectra for each gcm for
#         # prespawn survival, egg survival and fry survival
#         # Calculate wavelet spectra
#         # Visualize spectra
#         dims <- dim(prespSurvC)
#         out_wavelet <- file.path(".", "output", "wavelet_plots")
#         if (!dir.exists(out_wavelet)) { dir.create(out_wavelet, recursive = TRUE) }
#         
#         path4file <- file.path(out_wavelet, paste0(mgmtScenarios[m], "_", climateScenarios[c], ".pdf"))
#         pdf(file = path4file)
#         
#         for (i in 2:dims[2]) {
#           #i <- 2
#           # Time series
#           par(fig= c(0,0.65,0.6,1))
#           old <- par(mar = c(3, 3, 3, 1) )
#           plot(2010:2099, prespSurvC[,i], type = "l", lwd = 2, col = "darkgrey",
#                axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
#           axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                las = 1, cex.axis = 1, tck = 0.02)
#           axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
#                lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
#                las = 2, cex.axis = 1, tck = 0.02)
#           box(lwd = 2)
#           title(paste("GCM: ", gcms[i-1], " Prespawn survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))
#           par(old)
#           
#           # Wavelet Power Spectrum
#           par(fig= c(0,0.65,0.3,0.6), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(prespSurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#           plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
#           axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                las = 1, cex.axis = 1, tck = 0.02)
#           par(old)
#           # Scale-averaged Power Spectrum - periods of greatest variability
#           par(fig= c(0,0.65,0,0.3), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
#           par(old)
#           # Global Wavelet Power Spectrum
#           par(fig= c(.65,1,0.3,.6), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           yrange <- NULL 
#           y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
#           FR <- log(rowSums(test_cwt$power.corr))
#           plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
#           axis(1)
#           axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
#           box()
#           par(old)
#         }
#         for (i in 2:dims[2]) {
#           #i <- 2
#           par(fig= c(0,0.65,0.6,1))
#           old <- par(mar = c(3, 3, 3, 1) )
#           # plot egg surv ts
#           plot(2010:2099, eggSurvC[,i], type = "l", lwd = 2, col = "darkgrey",
#                axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
#           axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                las = 1, cex.axis = 1, tck = 0.02)
#           axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
#                lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
#                las = 2, cex.axis = 1, tck = 0.02)
#           box(lwd = 2)
#           title(paste("GCM: ", gcms[i-1], " Eggs survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))
#           par(old)
#           # plot wavelet power spectrum of egg survival
#           par(fig= c(0,0.65,0.3,0.6), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(eggSurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#           plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
#           axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                las = 1, cex.axis = 1, tck = 0.02)
#           par(old)
#           # Scale-averaged Power Spectrum - periods of greatest variability
#           par(fig= c(0,0.65,0,0.3), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
#           par(old)
#           # Global Wavelet Power Spectrum
#           par(fig= c(.65,1,0.3,0.6), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           yrange <- NULL 
#           y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
#           FR <- log(rowSums(test_cwt$power.corr))
#           plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
#           axis(1)
#           axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
#           box()
#           par(old)
#         }
#         
#         for (i in 2:dims[2]) {
#           #i <- 2
#           # plot fry survival time series 
#           par(fig= c(0,0.65,0.6,1))
#           old <- par(mar = c(3, 3, 3, 1) )
#           plot(2010:2099, frySurvC[,i], type = "l", lwd = 2, col = "darkgrey",
#                axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
#           axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                las = 1, cex.axis = 1, tck = 0.02)
#           axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
#                lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
#                las = 2, cex.axis = 1, tck = 0.02)
#           box(lwd = 2)
#           title(paste("GCM: ", gcms[i-1], " Fry survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))        
#           par(old)
#           # plot Wavelet power spectrum
#           par(fig= c(0,0.65,0.3,0.6), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(frySurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#           plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
#           axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
#                las = 1, cex.axis = 1, tck = 0.02)
#           par(old)
#           # Scale-averaged Power Spectrum - periods of greatest variability
#           par(fig= c(0,0.65,0,0.3), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
#           par(old)
#           # Global Wavelet Power Spectrum
#           par(fig= c(.65,1,0.3,.6), new = TRUE)
#           old <- par(mar = c(3, 3, 1, 1) )
#           yrange <- NULL 
#           y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
#           FR <- log(rowSums(test_cwt$power.corr))
#           plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
#           axis(1)
#           axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
#           box()
#           par(old)
#         }
#         
#         
#         dev.off()
#         
#       } # end of c climate emissions loop
#     } # end of m water management loop
#     
#   } # end of wavelet_SALMOD_series()
# 
# 
# ## plot 
# 
# wavelet_SALMOD_series() # Plots ts, wavelet power spectrum, plus the frequency-averaged Global WPS and
# # the time-averaged Scaled WPS - need to check on how to treat/scale magnitudes 
