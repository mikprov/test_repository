# RUN
# by: mikaela Mikaelast
# last edited: Aug 1, 2019
# ===================================================================


library(tidyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(directlabels)
library(grid)
library(lattice)
library(tidyverse)
library(tidycensus)
library(ggrepel)
library(biwavelet)
library(data.table)
#install.packages("ggpubr")
library(ggpubr)
# *************************************** #
# Plan:
# (1) run script to calculate life table information and load parms df
# (2) assemble Leslie matrices for different F values
# (3) simulate time series 
# (4) calculate spectra



# *************************************** #
# (1) run script to get life table information, parms df, and functions
source("C:/Users/Mikaela/Documents/GitHub/pfmc/src/life_table_calculations.r")


# *************************************** #
# (2) assemble Leslie matrices for different F values

# make all LEPs equal
conLEP = 5 #1.1
# adjust fecundities by this much
adjFec = round(1/conLEP,digits=1)

# Find F values that correspond to Exploitation Index values (BEFORE LEP=1.1)
# EI_F_before1.1_L <- as.list(rep(NA,length=length(parms$spp)))
# names(EI_F_before1.1_L) <- parms$spp
# system.time(for (i in 1:length(EI_F_before1.1_L)){ #for each spp
#   if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes") {#if a problem, then
#     EI_F_before1.1_L[[i]] <- "problem"} #temporarily save in problemhold
#   else {
#     targetF <- calculate_F_from_EI(EI=c(0,0.4,0.6,0.8,0.95),
#                                maxage=as.numeric(parms[parms$spp == parms$spp[i],]$maxage),
#                                L_a=La_list[[i]],
#                                B_a=Ba_list[[i]],
#                                propmat_a=prop_list[[i]],
#                                eggs_a=eggs_list[[i]],
#                                vul_a=vul_list[[i]],
#                                M=as.numeric(parms[parms$spp == parms$spp[i],]$M),
#                                Fvals_many=seq(from=0,to=100,by=0.1))
#     targetF$spp <- rep(parms$spp[i],length=length(targetF[,1]))
#     EI_F_before1.1_L[[i]] <- targetF
#   }
# })
# EI_F_before1.1 <- bind_rows(EI_F_before1.1_L)

# ***
#parms <- parms[parms$problem_spp=="no",]
Fvals=seq(from=0,to=2,by=0.001)
A3d_list = as.list(rep(NA,length(parms$spp))) #Leslie arrays in list
LTABLE_list <- as.list(rep(NA,length(parms$spp))) #LTABLE arrays in list
names(A3d_list) <- parms$spp
names(LTABLE_list) <- parms$spp
eigan1 <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
eigan12 <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
LEPs <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
constLEP <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
EI_F_LEP1.1_L <- as.list(rep(NA,length=length(parms$spp)))
names(EI_F_LEP1.1_L) <- parms$spp


system.time(for (i in 1:length(parms$spp)){ #for each spp
  
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes") {#if a problem, then
    problemshold <- parms$spp[i]} #temporarily save in problemhold
  else {
  # define parms for spp i
  M = as.numeric(parms[parms$spp == parms$spp[i],]$M)
  maxage = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  
  # create empty arrays for Leslie and LTABLE, vectors for eigens 
  A3d <- array(NA,c(maxage,maxage,length(Fvals)))
  LTABLE3d <- array(NA,c(maxage,15,length(Fvals))) #12 cols
  e1 <- rep(NA,length=length(Fvals))
  e2 <- rep(NA,length=length(Fvals))
  e12 <- rep(NA,length=length(Fvals))
  lep <- rep(NA,length=length(Fvals))
  
  # assemble Leslie for each value of Fvals
  for(f in 1:length(Fvals)){
    out <- assemble_Leslie(maxage=maxage,
                           L_a=La_list[[i]],
                           B_a=Ba_list[[i]],
                           propmat_a=prop_list[[i]],
                           eggs_a=eggs_list[[i]],
                           vul_a=vul_list[[i]],
                           M=M,
                           Fval=Fvals[f]) #F=0, change this if want fishing
    LTABLE3d[,,f] <- data.matrix(out$LTABLE) #store LTABLE in 3d array
    
    # calculate LEP at F=0 (use this value of LEP to standardize unfished fecundity)
    lep_unfished = sum(out$LTABLE$LEP_a)
    
    # Asim = adjust fecundities so LEP equal for all species - simulation matrix
    Asim <- out$A
    Asim[1,] <- (Asim[1,]/(lep_unfished*adjFec))
    A3d[,,f] <- Asim
    
    # Aeig = adjust fecundities and multiply by k - eigen analysis matrix
    Aeig <- out$A
    Aeig[1,] <- (Aeig[1,]/(lep_unfished*adjFec)) # *Fvals[f]
    
    # re-calculate new LEP w/o F (should equal conLEP)
    term <- rep(NA,length=maxage)
    for(a in 1:maxage){term[a] <- Asim[1,a]*out$LTABLE$Survship[a]} #fecundity*survivorship
    constLEP[f,i] <- sum(term) #plotting 'term' is the spawning distribution over age
    rm(term,a)
    
    # calculate LEP w/F
    term <- rep(NA,length=maxage)
    for(a in 1:maxage){term[a] <- Asim[1,a]*out$LTABLE$Survship_Fing[a]} #fec*survshp w/F
    LEPs[f,i] <- sum(term)
    rm(term,a)
    
    # extract lambda1, lambda2/lambda1
    e1[f] = extract_first_eigen_value(Aeig)
    e2[f] = extract_second_eigen_value(Aeig)
    e12[f] = e2[f] / e1[f]
    
    
  } #closes F loop
  EI_F_LEP1.1_L[[i]] <- calculate_F_from_EI_Leslie(A=A3d[,,1],
                             vul_a=vul_list[[i]],
                             Fval=Fvals,
                             M=as.numeric(parms$M[i]),
                             EI=c(0,0.4,0.6,0.8))
  #dimnames(A3d)[[3]] <-as.character(Fvals) #label the 3rd dim by Fvals
  saveFs <- which(Fvals %in% EI_F_LEP1.1_L[[i]]$matching_F) #where are Fs of interest?
  A3dsave <- A3d[,,saveFs] #save only Leslie matrices that match Fs
  
  A3d_list[[i]] <- A3dsave 
  LTABLE_list[[i]] <- LTABLE3d
  eigan1[,i] <- e1
  eigan12[,i] <- e12
  LEPs[,i] <- lep
  #clean up before next i in loop
  #rm(A3d,LTABLE3d,out,M,maxage) 
  } #closes 'else' bracket for non-problem species.
  print(i)
}) #closes spp loop
rm(i,f,M,maxage,problemshold,lep,A3d,A3dsave,saveFs) #clean up


# *************************************** #
# () Calculate metrics from spawning distribution
# *************************************** #
# script below generates df 'spawndistmetrics' which has metric info
# keep only non-problem spp, convert numeric variables to numbers
# add metric info to ts.data
source("C:/Users/Mikaela/Documents/GitHub/pfmc/src/spawning_distributions.r")
spawndistmetrics <- spawndistmetrics[spawndistmetrics$problem_spp == "no",]
spawndistmetrics$mode_age <- as.numeric(as.character(spawndistmetrics$mode_age))
spawndistmetrics$sd_mode <- as.numeric(as.character(spawndistmetrics$sd_mode))
spawndistmetrics$cvs_mode <- as.numeric(as.character(spawndistmetrics$cvs_mode))
spawndistmetrics$maxage <- as.numeric(as.character(spawndistmetrics$maxage))
spawndistmetrics$spp <- as.character(spawndistmetrics$spp)

#tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/spawning_distribution_v2.tiff', units="in", width=12, height=9, res=300)
#do.call(grid.arrange,c(psub,ncol=2,left="Pr(spawning)"))
#dev.off()


peak_max <- ggplot(data=spawndistmetrics,aes(x=maxage,y=mode_age)) +
  geom_point() +
  geom_text_repel(data=spawndistmetrics,
                  aes(label = spp),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE)+
  theme_classic() + ylab("Peak spawning age") + xlab("Maximum age") 
  
sd_max <- ggplot(data=spawndistmetrics,aes(x=maxage,y=sd_mode)) +
  geom_point() +
  geom_text_repel(data=spawndistmetrics,
                  aes(label = spp),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE)+
  theme_classic() + ylab("Spawning distribution Stdev") + xlab("Maximum age") 
rm(peak_max,sd_max)

std_peak <- ggplot(data=spawndistmetrics,aes(x=mode_age,y=sd_mode)) +
  geom_point() +
  geom_text_repel(data=spawndistmetrics,
                  aes(label = spp),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE)+
  theme_classic() + ylab("Spawning distribution Stdev") + xlab("Peak spawning age")

cv_peak <- ggplot(data=spawndistmetrics,aes(x=mode_age,y=cvs_mode)) +
  geom_point() +
  geom_text_repel(data=spawndistmetrics,
                  aes(label = spp),
                  segment.color = "grey",
                  size = 2,
                  na.rm = TRUE)+
  theme_classic() + ylab("Spawning distribution CV") + xlab("Peak spawning age")

rm(std_peak,cv_peak)

# *************************************** #
# (3) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# Generate 3d Leslie arrays for EIlevels
FLEPlevels = c(1,0.8,0.6,0.35,0.3)
EIlevels=round((1-FLEPlevels)/0.9,digits=2)
#EIlevels=c(0,0.3,0.5,0.76,0.77,0.78) #last number corresponds to FLEP=0.3
#FLEPlevels=1-(0.9*EIlevels)
#FLEPlevels=round(1-(EIlevels*0.9),digits=2)

A3d_list <- generate3dLeslie(spp=parms[parms$problem_spp=="no",]$spp,
                             Fvals=Fvals,
                             EIlevels=EIlevels)



# Generate noise time series 
set.seed(62) #Set the seed so that every white ts uses same random sequence
timesteps=5000
rm_first_timesteps = 2000
# Noise signals
white <- rnorm((timesteps-2),mean=0,sd=1) # white noise
enso_ts_reps <- read.csv(file="C:/Users/Mikaela/Documents/GitHub/pfmc/parms_data/simulated_MEI.csv",header=TRUE) 
ensoS_ts_reps <- read.csv(file="C:/Users/Mikaela/Documents/GitHub/pfmc/parms_data/simulated_MEI_lowfreq.csv",header=TRUE) 
ensoF_ts_reps <- read.csv(file="C:/Users/Mikaela/Documents/GitHub/pfmc/parms_data/simulated_MEI_highfreq.csv",header=TRUE) 

# make sure same length as white
enso1 <- enso_ts_reps[,10][1:length(white)]
ensoF <- ensoF_ts_reps[,10][1:length(white)]
ensoS <- ensoS_ts_reps[,10][1:length(white)]
noisenames <- c("white","ENSO","ENSOfast","ENSOslow")
noisetypes <- list(white,enso1,ensoF,ensoS)
names(noisetypes) <- noisenames
# check total variance in signals
enso1 <- (enso1-mean(enso1))/sd(enso1)
ensoF <- (ensoF-mean(ensoF))/sd(ensoF)
ensoS <- (ensoS-mean(ensoS))/sd(ensoS)
# mean(enso1)/sd(enso1)
# mean(ensoF)/sd(ensoF)
# mean(ensoS)/sd(ensoS)
# sd(enso1)
# sd(ensoF)
# sd(ensoS)
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/noise_white.tiff', units="in", width=5, height=4, res=300) 
# plot(white[1:100],type="l",ylim=c(-3,3),main="white",xlab="Year",ylab="")
# dev.off()
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/noise_enso.tiff', units="in", width=5, height=4, res=300) 
# plot(enso1[1:100],type="l",ylim=c(-3,3),main="ENSO (historical)",xlab="Year",ylab="")
# dev.off()
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/noise_ensoF.tiff', units="in", width=5, height=4, res=300) 
# plot(ensoF[1:100],type="l",ylim=c(-3,3),main="ENSO 2x",xlab="Year",ylab="")
# dev.off()
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/noise_ensoS.tiff', units="in", width=5, height=4, res=300) 
# plot(ensoS[1:100],type="l",ylim=c(-3,3),main="ENSO 0.5x",xlab="Year",ylab="")
# dev.off()

# load simulation model:
source("C:/Users/Mikaela/Documents/GitHub/pfmc/src/simulation_model.r")

# Generate table of ts for plotting 
# note: alpha > 1/LEP. Since I set LEP=5 for all spp, alpha must be > 1/5
ts.L <- as.list(rep(NA,length=length(noisenames)))

for(n in 1:length(noisenames)){
  ts <- simulate_spp(noise = noisetypes[[n]],
                        timesteps = timesteps,
                        rm_first_timesteps = rm_first_timesteps, #burn in time
                        alphas = 1, #for BH
                        betas = 5000, #for BH
                        span.multiplier = 1, #not sure I need this yet
                        EIlevels = EIlevels, #vector of EI levels, 1st is OFL--match
                        OFL = max(EIlevels), #FLEP associated with OFL
                        List_of_Leslie3d = A3d_list, #should have F vals associated with FLEPlevels
                        spp = parms[parms$problem_spp == "no",]$spp,
                        sig_r=0.5)
  ts$noise <- rep(noisenames[n],length=length(ts[,1]))
  ts.L[[n]] <- ts   }
# convert list into one data frame
ts.data <- bind_rows(ts.L,id=NULL)
# add column with corresponding FLEP
ts.data$FLEP <- round(1-(ts.data$EI*0.9),digits=2)
rm(ts.L,n,ts)
spp = parms[parms$problem_spp == "no",]$spp
# -------------------------------------------------------
# Plot ts w/white and ENSO noise and threshold
# ts is with harvest rate such that FLEP=0.4 or 0.45 
# threshold = mean level of eggs/recruits/N when FLEP=0.3
# #ts.data$EI <- as.character(ts.data$EI)
# 
# multiplots <- as.list(rep(NA,length(spp)))
# names(multiplots) <- spp
#spp = parms[parms$spp == "Darkblotched",]$spp
#minSSB = 15000
#maxSSB = 50000
# for(n in 1:length(noisenames)){
#   multiplots <- as.list(rep(NA,length(spp)))
#   names(multiplots) <- spp
#   dfts <- ts.data[ts.data$FLEP==0.35 &
#                 ts.data$variable=="Nsize" &
#                 ts.data$noise==noisenames[n] &
#                 ts.data$year %in% rm_first_timesteps:(timesteps-2600),]
#   dfth <- ts.data[ts.data$FLEP==0.3 &
#                 ts.data$variable=="Nsize" &
#                 ts.data$noise==noisenames[n] &
#                 ts.data$year %in% rm_first_timesteps:(timesteps-2600),]
#   minSSB <- if(dfth$means[1] < min(dfts$value)){dfth$means[1]}else{min(dfts$value)}
#   maxSSB <- max(dfts$value)


  
# -------------------------------------------------------
# Calculate probability of exceeding overfishing limit
# for each variable type in both types of noise
probs <- array(NA,c(length(unique(ts.data$noise)),
                    length(unique(ts.data$variable)),
                    length(unique(ts.data$spp)))) 

for(s in 1:length(unique(ts.data$spp))) {
  for(n in 1:length(unique(ts.data$noise))) {
    for(v in 1:length(unique(ts.data$variable))) {
      #subset the right time series
      vals <- ts.data[ts.data$spp==unique(ts.data$spp)[s] &
        ts.data$variable==unique(ts.data$variable)[v] &
        ts.data$noise==unique(ts.data$noise)[n] &
        ts.data$FLEP==0.35 &
        ts.data$year %in% rm_first_timesteps:(timesteps-2600),]$value
      #pull out the right threshold value
      threshold <- ts.data[ts.data$spp==unique(ts.data$spp)[s] &
        ts.data$variable==unique(ts.data$variable)[v] &
        ts.data$noise==unique(ts.data$noise)[n] &
        ts.data$FLEP==0.3,]$means[1] #only need 1 value

      # calc probability: # of timesteps below threshold/total timesteps
      probs[n,v,s] <- sum(vals < threshold) / length(vals)
    }
  }
}
rm(s,n,v,threshold,vals)
variable_type <- unique(ts.data$variable)
ndf <- as.list(rep(NA,length(unique(ts.data$noise))))
for(n in 1:length(unique(ts.data$noise))){
  df <- as.data.frame(probs[n,,]) #select rec,egg,N ts for all spp in n noise
  df <- bind_cols(as.data.frame(variable_type),
                  as.data.frame(rep(unique(ts.data$noise)[n],length(df[,1]))),
                  df)
  colnames(df) <-  c("variable","noise",spp)
  df <- df %>% gather(key="spp",value="probs",3:ncol(df))
  ndf[[n]] <- df}
probdf <- bind_rows(ndf,id=NULL)
rm(n,df,ndf,probs)
# fill in probdf with age structure characteristics
probdf$peak <- spawndistmetrics[match(probdf$spp,spawndistmetrics$spp),"mode_age"]
probdf$cvs <- spawndistmetrics[match(probdf$spp,spawndistmetrics$spp),"cvs_mode"]
probdf$maxage <- spawndistmetrics[match(probdf$spp,spawndistmetrics$spp),"maxage"]
probdf$stdev <- spawndistmetrics[match(probdf$spp,spawndistmetrics$spp),"sd_mode"]



# # Plot p(OFL) vs max age, peak, CV, and stdev of spawning biomass distribution
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/agestructure_vs_probs_Nsize.tiff', units="in", width=8, height=11, res=300) 
# p1 <- ggplot(data=probdf[probdf$variable=="Nsize" ,],aes(x=maxage,y=probs,color=spp)) +
#   geom_point() + facet_grid(vars(noise)) + ggtitle("(1a) Fraction of SSB time series\n that is below mean recruitment\n when FLEP=0.3 (the assumed OFL) vs max age") +
#   theme(plot.title = element_text(size=8))
# p2 <- ggplot(data=probdf[probdf$variable=="Nsize" ,],aes(x=maxage,y=probs,color=noise)) +
#   geom_point()  + ggtitle("(1b) Similar plot as 1a, but color-coded by noise")+
#   theme(plot.title = element_text(size=8))
# 
# p3 <- ggplot(data=probdf[probdf$variable=="Nsize" ,],aes(x=cvs,y=probs,color=spp)) +
#   geom_point() + facet_grid(vars(noise)) + ggtitle("(2a) Fraction of SSB time series\n that is below mean recruitment\n when FLEP=0.3 (the assumed OFL) vs CV")+
#   theme(plot.title = element_text(size=8))
# p4 <- ggplot(data=probdf[probdf$variable=="Nsize" ,],aes(x=cvs,y=probs,color=noise)) +
#   geom_point()  + ggtitle("(2b) Similar plot as 2a, but color-coded by noise")+
#   theme(plot.title = element_text(size=8))
# 
# p5 <- ggplot(data=probdf[probdf$variable=="Nsize" ,],aes(x=peak,y=probs,color=spp)) +
#   geom_point() + facet_grid(vars(noise)) + ggtitle("(3a) Fraction of SSB time series\n that is below mean recruitment\n when FLEP=0.3 (the assumed OFL) vs peak age")+
#   theme(plot.title = element_text(size=8))
# p6 <- ggplot(data=probdf[probdf$variable=="Nsize" ,],aes(x=peak,y=probs,color=noise)) +
#   geom_point()  + ggtitle("(3b) Similar plot as 3a, but color-coded by noise")+
#   theme(plot.title = element_text(size=8))
# 
# pplots <- list(p1,p2,p3,p4,p5,p6)
# do.call(grid.arrange,c(pplots,ncol=2))
# dev.off()

# # -------------------------------------------------------
# # Test: do the raw timeseries of biomass and timeseries
# # when CV has been standardized to 1 look the same?
# # ANSWER: Yes. Note=can't plot log() scale
# tempts <- ts.data[ts.data$spp==spp[1] &
#                     ts.data$FLEP==1 &
#                     ts.data$variable=="recruits" &
#                     ts.data$noise %in% c("white") &
#                     ts.data$year %in% rm_first_timesteps:(timesteps-2),]
# values.mean0 <- (tempts$value-mean(tempts$value))/sd(tempts$value)
# tempts.mean0 <- as.data.frame(cbind(tempts$year,values.mean0))
# colnames(tempts.mean0) <- c("year","value")
# head(tempts.mean0)
# temp1 <- tempts.mean0
# temp1 <- tempts.mean0
# # ggplot(data=tempts.mean0,aes(x=year,y=value)) +
# #       geom_line() +
# #       
# #       ylab("Spawning stock biomass") +
# #       theme_bw() 
# sptest <- spec.pgram(temp1$value,spans=c(m,m),
#                         taper=0.1,plot=FALSE,demean=TRUE)
# d1 <- as.data.frame(cbind(sptest$freq,sptest$spec))
# colnames(d1) <- c("freq","sp")
# sptest <- spec.pgram(temp12$value,spans=c(m,m),
#                         taper=0.1,plot=FALSE,demean=TRUE)
# d2 <- as.data.frame(cbind(sptest$freq,sptest$spec))
# colnames(d2) <- c("freq","sp")
# 
# ggplot(data=d1,aes(x=freq, y=sp)) +
#   geom_line() + scale_y_log10()
#   geom_line(data=d1,aes(x=freq,y=sp))+
#   scale_y_log10()
# -------------------------------------------------------
# Calculate spectra in time series
# 1) create df w/spec for different noises
ts.for.spec <- ts.data
ts.for.spec$spec <- rep(NA,length(ts.for.spec[,1]))
ts.for.spec <- ts.for.spec[ts.for.spec$year %in% rm_first_timesteps:(timesteps-2),]
span.multiplier= 1
# setting 'span' - a vector of odd integers to specify the smoothers
# sq root of timeseries lgth, rounded
# make it odd, if the square root is even
tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-2)))) 
if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} 
m = m * span.multiplier

df.sppL <- as.list(rep(NA,length(unique(ts.for.spec$spp))))
names(df.sppL) <- unique(ts.for.spec$spp)
# for each species...
for (s in 1:length(unique(ts.for.spec$spp))){
  df.EIL <- as.list(rep(NA,length(unique(ts.for.spec$EI))))
  names(df.EIL) <- unique(ts.for.spec$EI)
  # and for each exploitation level...
  for(e in 1:length(unique(ts.for.spec$EI))){
    df.varL <- as.list(rep(NA,length(unique(ts.for.spec$variable))))
    names(df.varL) <- unique(ts.for.spec$variable)
    # and for each type of output data (eggs, recruits, N)...
    for(v in 1:length(unique(ts.for.spec$variable))){
      df.noiseL <- as.list(rep(NA,length(unique(ts.for.spec$noise))))
      names(df.noiseL) <- unique(ts.for.spec$noise)
      # and for each type of noise...
      # calculate the spectrum in the time series
      for(n in 1:length(unique(ts.for.spec$noise))){
        ts <- ts.for.spec[ts.for.spec$spp==unique(ts.for.spec$spp)[s] &
                      ts.for.spec$EI==unique(ts.for.spec$EI)[e] &
                      ts.for.spec$variable==unique(ts.for.spec$variable)[v] &
                      ts.for.spec$noise==unique(ts.for.spec$noise)[n],]$value
        
        # ---
        # Chunk 1: Run if I want to standardize CV in ts
        ts.mod <- (ts-mean(ts))/sd(ts) 
        # ---
        # Chunk 2: Run if I DON'T want to standardize CV
        #ts.mod <- ts
        # ---
        sp = spec.pgram(ts.mod,spans=c(m,m),
                        taper=0.1,plot=FALSE,demean=TRUE)
        #divide each value of spec by the area under the curve
        specN <- as.numeric(sp$spec/sum(sp$spec)) # sum(specN)=1
        spec <- sp$spec # sum(spec) not equal 1
        freq <- sp$freq
        noise <- rep(unique(ts.for.spec$noise)[n],length(freq))
        df=data.frame(cbind(spec,specN,freq,noise))
        df.noiseL[[n]]=df
        print(unique(ts.for.spec$noise)[n])
      }#close noise loop
      dfv <- bind_rows(df.noiseL,.id=NULL)
      dfv$variable <- rep(unique(ts.for.spec$variable)[v],length(dfv[,1]))
      df.varL[[v]]<-dfv
      print(unique(ts.for.spec$variable)[v])
    }#close variable loop
    dfe <- bind_rows(df.varL,.id=NULL)
    dfe$EI <- rep(unique(ts.for.spec$EI)[e],length(dfe[,1]))
    df.EIL[[e]]<-dfe
    print(unique(ts.for.spec$EI)[e])
  }#close EI loop
  dfs <- bind_rows(df.EIL,.id=NULL)
  dfs$spp <- rep(unique(ts.for.spec$spp)[s],length(dfs[,1]))
  df.sppL[[s]]<-dfs
  print(unique(ts.for.spec$spp)[s])
}
ts.spec <- bind_rows(df.sppL,.id=NULL) 
ts.spec$spec <- as.numeric(ts.spec$spec)
ts.spec$specN <- as.numeric(ts.spec$specN)
ts.spec$freq <- as.numeric(levels(ts.spec$freq))[ts.spec$freq]
# add spp age structure characteristics to table
ts.spec$peak <- spawndistmetrics[match(ts.spec$spp,spawndistmetrics$spp),"mode_age"]
ts.spec$maxage <- spawndistmetrics[match(ts.spec$spp,spawndistmetrics$spp),"maxage"]
ts.spec$cvs <- spawndistmetrics[match(ts.spec$spp,spawndistmetrics$spp),"cvs_mode"]
ts.spec$stdev <- spawndistmetrics[match(ts.spec$spp,spawndistmetrics$spp),"sd_mode"]
rm(dfe,dfs,df.sppL,df.EIL,dfv,df.varL)

# # 2) calculate spec of ENSO
# ensosp <- spec.pgram(x=enso1[rm_first_timesteps:(timesteps-2)],
#                        spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
# # check length of spec for enso, should be same length as spec for spp
# if(length(ensosp$spec)==length(ts.spec[ts.spec$spp=="Sardine" &
#                                    ts.spec$EI==unique(ts.spec$EI)[1] &
#                                    ts.spec$noise==unique(ts.spec$noise)[1] &
#                                    ts.spec$variable==unique(ts.spec$variable)[1],1]) ){print("yes")}else{print("no")}
# # plot ENSO frequency content
# plot(x=ensosp$freq,y=ensosp$spec,xaxs = "i", yaxs = "i",
#      type = "l",log="y",yaxt="n",ylab="")
# 
# # 3) ENSO spectrum x pop spectrum in white noise = pop spectrum in ENSO
# multiplied.spec <- ts.spec[!ts.spec$noise=="ENSO",]
# multiplied.spec$noise <- rep("Spec w/white x ENSO",length(multiplied.spec[,1]))
# multiplied.spec$spec <- rep(NA,length(multiplied.spec[,1]))
# head(multiplied.spec)
# 
# for(s in 1:length(unique(ts.spec$spp))){
#   for(e in 1:length(unique(ts.spec$EI))){
#     for(v in 1:length(unique(ts.spec$variable))){
#       wh.spec <- ts.spec[ts.spec$spp==unique(ts.spec$spp)[s] &
#                          ts.spec$EI==unique(ts.spec$EI)[e] &
#                          ts.spec$variable==unique(ts.spec$variable)[v] &
#                          ts.spec$noise=="white",]$spec
#       multiplied.spec[multiplied.spec$spp==unique(multiplied.spec$spp)[s] &
#                       multiplied.spec$EI==unique(multiplied.spec$EI)[e] &
#                       multiplied.spec$variable==unique(multiplied.spec$variable)[v],]$spec <- wh.spec * ensosp$spec
#     }  } } 
# rm(s,e,v)
# multiplied.spec <- as.data.frame(multiplied.spec)
# ts.spec.combo <- rbind(ts.spec,multiplied.spec)
# 
# # 4) Plot overlapping spectrum
# overlap.specL <- as.list(rep(NA,length(unique(multiplied.spec$spp))))
# for(s in 1:length(unique(ts.spec.combo$spp))){
#   overlap.specL[[s]] <- ggplot(data=ts.spec.combo[ts.spec.combo$spp==unique(ts.spec.combo$spp)[s] &
#     ts.spec.combo$EI==0 &
#     ts.spec.combo$variable=="recruits" &
#     ts.spec.combo$noise %in% c("ENSO","Spec w/white x ENSO"),],
#     aes(x=freq,y=spec,color=noise)) +
#     geom_line() + scale_color_grey() +  
#     scale_y_log10(limits=c(min(ts.spec.combo$spec),max(ts.spec.combo$spec)),
#                   breaks=c(1e1,1e7)) +
#     ggtitle(unique(ts.spec.combo$spp)[s]) + 
#     theme_classic() + ylab("") + xlab("") +
#     theme(legend.title = element_blank(),legend.position=c(0.7,0.9))   }#end ggplot loop
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/Overlapping_spectra.tiff', units="in", width=7, height=13, res=300) 
# do.call(grid.arrange,c(overlap.specL,ncol=2,left="Sensitivity", bottom="Frequency"))
# dev.off()
# 
# # 5) Plot overlapping SSB spectra in white & ENSO 
# overlap.specL <- as.list(rep(NA,length(unique(ts.spec$spp))))
# names(overlap.specL) <- unique(ts.spec$spp)
# maxspec <- max(ts.spec[ts.spec$EI==0 &
#     ts.spec$variable=="Nsize" &
#     ts.spec$noise %in% c("white","ENSO"),]$specN)
# minspec <- min(ts.spec[ts.spec$EI==0 &
#     ts.spec$variable=="Nsize" &
#     ts.spec$noise %in% c("white","ENSO"),]$specN)
# for(s in 1:length(unique(ts.spec$spp))){
#   ts <- ts.spec[ts.spec$spp==unique(ts.spec$spp)[s] &
#     ts.spec$EI==0 &
#     ts.spec$variable=="Nsize" &
#     ts.spec$noise %in% c("white","ENSO"),]
#   overlap.specL[[s]] <- ggplot(data=ts,aes(x=freq,y=specN,linetype=noise)) +
#     geom_line() + scale_color_grey() +
#     scale_y_continuous(limits=c(minspec,maxspec),trans="log10") +
#     ggtitle(unique(ts.spec$spp)[s]) +
#     theme_classic() + ylab("") + xlab("") +
#     theme(legend.title = element_blank(),legend.position=c(0.7,0.9))
# }#end ggplot loop
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/Overlapping_spectra_whiteENSO_all10.tiff', units="in", width=6, height=12, res=300) 
# do.call(grid.arrange,c(overlap.specL,ncol=2,left="Sensitivity", bottom="Frequency"))
# dev.off()
# rm(overlap.specL,s,ts,maxspec,minspec)

# # 5) Plot overlapping SSB spectra in ENSO, ENSOfast, ENSOslow 
# overlap.specL <- as.list(rep(NA,length(unique(ts.spec$spp))))
# names(overlap.specL) <- unique(ts.spec$spp)
# maxspec <- max(ts.spec[ts.spec$EI==0 &
#     ts.spec$variable=="Nsize" &
#     ts.spec$noise %in% c("ENSO","ENSOfast","ENSOslow"),]$specN)
# minspec <- min(ts.spec[ts.spec$EI==0 &
#     ts.spec$variable=="Nsize" &
#     ts.spec$noise %in% c("ENSO","ENSOfast","ENSOslow"),]$specN)
# for(s in 1:length(unique(ts.spec$spp))){
#   ts <- ts.spec[ts.spec$spp==unique(ts.spec$spp)[s] &
#     ts.spec$EI==0 &
#     ts.spec$variable=="Nsize" &
#     ts.spec$noise %in% c("ENSO","ENSOfast","ENSOslow"),]
#   overlap.specL[[s]] <- ggplot(data=ts,aes(x=freq,y=specN,linetype=noise)) +
#     geom_line() + scale_color_grey() +
#     scale_y_continuous(limits=c(minspec,0.021),trans="log10") +
#     ggtitle(unique(ts.spec$spp)[s]) +
#     theme_classic() + ylab("") + xlab("") +
#     theme(legend.title = element_blank(),legend.position=c(0.7,1))
# }#end ggplot loop
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/Overlapping_spectra_ENSOfs_all10.tiff', units="in", width=8, height=10, res=300) 
# do.call(grid.arrange,c(overlap.specL,ncol=3,left="log10", bottom="Frequency",top="Overlapping SSB spectrum with ENSO, ENSO-fast, and ENSO-slow \n(SSB CV is constant across spp)"))
# dev.off()
# rm(overlap.specL,s,ts,maxspec,minspec)

# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/Overlapping_spectra_whiteENSO_all10_baseR.tiff', units="in", width=6, height=12, res=300) 
# par(mfrow=c(6,2),mai=c(0.2,0.3,0.2,0.2))
# for(i in 1:length(spp)){
#   plot(x=ts.spec[ts.spec$spp==spp[i] &
#                      ts.spec$noise=="ENSO" &
#                      ts.spec$variable=="Nsize" &
#                      ts.spec$EI==0,]$freq,
#          y=log10(ts.spec[ts.spec$spp==spp[i] &
#                      ts.spec$noise=="ENSO" &
#                      ts.spec$variable=="Nsize" &
#                      ts.spec$EI==0,]$specN),
#          type = "l",lty="solid",ylab="",xlab="frequency")
#   lines(x=ts.spec[ts.spec$spp==spp[i] &
#                      ts.spec$noise=="white" &
#                      ts.spec$variable=="Nsize" &
#                      ts.spec$EI==0,]$freq,
#          y=log10(ts.spec[ts.spec$spp==spp[i] &
#                      ts.spec$noise=="white" &
#                      ts.spec$variable=="Nsize" &
#                      ts.spec$EI==0,]$specN),
#          type = "l",lty="dashed",ylab="",xlab="frequency")
# 
# }
# dev.off()

# # -------------------------------------------------------
# # Plot wavelet multipanel plots - exports 2 images for each noise
# for(n in 1:length(noisenames)){
#   
#   plot_wvlt_CVnotstd_AUC1(burn_in_pd = rm_first_timesteps, #remove from beginning of ts
#                    num_rows2plt = 1000, #number of years to plot
#                    timesteps = timesteps,
#                    span.multiplier = 1,
#                    noise = noisetypes[[n]],
#                    noisename = noisenames[n],
#                    ts.data = ts.data, #timeseries dataset
#                    ts.spec = ts.spec, #spectrum dataset
#                    spawning_dist_data=spawning_dist_data)
#   print(noisenames[n])
#   
# }

# -------------------------------------------------------
# Calculate total variance in time series (stdev/mean)
# Use this to calculate absolute variance at high and low frequencies
fillList <- as.list(rep(NA,length(unique(ts.data$spp)))) 
names(fillList) <- unique(ts.data$spp)
system.time(for(s in 1:length(unique(ts.data$spp))){
  fillsdm <- array(NA,c(length(unique(ts.data$variable)),
                    length(unique(ts.data$noise)),
                    length(unique(ts.data$EI))))
  for(e in 1:length(unique(ts.data$EI))){
    for(n in 1:length(unique(ts.data$noise))){
      for(v in 1:length(unique(ts.data$variable))){
        stdev <- sd(ts.data[ts.data$spp==unique(ts.data$spp)[s] &
                            ts.data$EI==unique(ts.data$EI)[e] &
                            ts.data$noise==unique(ts.data$noise)[n] &
                            ts.data$variable==unique(ts.data$variable)[v],]$value[rm_first_timesteps:(timesteps-2)])
        m <- ts.data[ts.data$spp==unique(ts.data$spp)[s] &
                            ts.data$EI==unique(ts.data$EI)[e] &
                            ts.data$noise==unique(ts.data$noise)[n] &
                            ts.data$variable==unique(ts.data$variable)[v],]$mean[1]
        fillsdm[v,n,e] <- stdev/m
      }
    }
  }
  fillList[[s]] <- fillsdm
  print(unique(ts.data$spp)[s])
})

list.of.each.spp <- as.list(rep(NA,length(unique(ts.data$spp))))
names(list.of.each.spp) <- unique(ts.data$spp)
for(s in 1:length(fillList)){
  df <- fillList[[s]]
  list.of.each.EI <- as.list(rep(NA,length(unique(ts.data$EI))))
  names(list.of.each.EI) <- unique(ts.data$EI)
  
  for(e in 1:length(unique(ts.data$EI))){
    dfe <- as.data.frame(df[,,e])
    list.of.each.noise <- as.list(rep(NA,length(unique(ts.data$noise))))
    names(list.of.each.noise) <- unique(ts.data$noise)
    
    for(n in 1:length(unique(ts.data$noise))){
      dfn <- dfe[,n]
      d <- as.data.frame(cbind(unique(ts.data$variable),dfn,
                               rep(unique(ts.data$noise)[n])))
      colnames(d) <- c("variable","variance","noise")
      list.of.each.noise[[n]] <- d
      print(unique(ts.data$noise)[n])
    }
    dfnoise <- bind_rows(list.of.each.noise,.id=NULL)
    EI <- rep(unique(ts.data$EI)[e],length(dfnoise[,1]))
    dfnoise <- cbind(EI,dfnoise)
    list.of.each.EI[[e]] <- dfnoise
    print(unique(ts.data$EI)[e])
  }
  dfei <- bind_rows(list.of.each.EI,.id=NULL)
  spp <- rep(names(fillList)[s],length(dfei[,1]))
  dfei <- cbind(spp,dfei)
  list.of.each.spp[[s]] <- dfei
  print(unique(ts.data$spp)[s])
}
totalvar <- bind_rows(list.of.each.spp,.id=NULL)
# fill in totalvar with age structure characteristics
totalvar$peak <- spawndistmetrics[match(totalvar$spp,spawndistmetrics$spp),"mode_age"]
totalvar$cvs <- spawndistmetrics[match(totalvar$spp,spawndistmetrics$spp),"cvs_mode"]
totalvar$maxage <- spawndistmetrics[match(totalvar$spp,spawndistmetrics$spp),"maxage"]
totalvar$stdev <- spawndistmetrics[match(totalvar$spp,spawndistmetrics$spp),"sd_mode"]
totalvar$variance <- as.numeric(totalvar$variance)
rm(list.of.each.spp,dfei,spp,dfnoise,list.of.each.EI,list.of.each.noise,dfn)

# -------------------------------------------------------
# Calculate fraction & amount of variance at high & low freq

# 2 definitions for high and low frequencies: choose one!
# 1) threshold for high/low frequencies is 1/(2T) 
spawndistmetrics$AUCthresholdT <- 1/(spawndistmetrics$mode_age*2)
AUCthreshold_ordered_by_peak <- spawndistmetrics %>% arrange(mode_age) %>% pull(AUCthresholdT)
spp_ordered_by_peak <- spawndistmetrics %>% arrange(mode_age) %>% pull(spp)

# 2) threshold could be the same for all pops:
AUCthreshold10 <- rep(0.1,length=length(spawndistmetrics[,1]))
spawndistmetrics <- cbind(spawndistmetrics,AUCthreshold10)

# Fraction of H/L - AUC equals 1
# Amount of H/L - AUC equals amount, sum 'spec' for freq 
# Amount of H/L - multiply totalvar by fraction of AUC
AUC_amount_highL <- as.list(rep(NA,length=length(noisenames)))
AUC_amount_lowL <- as.list(rep(NA,length=length(noisenames)))
AUC_fraction_highL <- as.list(rep(NA,length=length(noisenames)))
AUC_fraction_lowL <- as.list(rep(NA,length=length(noisenames)))
TOT_highL <- as.list(rep(NA,length=length(noisenames)))
TOT_lowL <- as.list(rep(NA,length=length(noisenames)))
AUC_totalamt_L <- as.list(rep(NA,length=length(noisenames)))

names(AUC_amount_highL) <- noisenames
names(AUC_amount_lowL) <- noisenames
names(AUC_fraction_highL) <- noisenames
names(AUC_fraction_lowL) <- noisenames
names(AUC_totalamt_L) <- noisenames
names(TOT_highL) <- noisenames
names(TOT_lowL) <- noisenames

# frequency increament
freq <- min(ts.spec$freq)

for (n in 1:length(noisenames)){ #for each noise type...
  
  #empty dfs for high var at high freq
  AUC_amount_high <- rep(NA,length=length(spp_ordered_by_peak))
  AUC_fraction_high <- rep(NA,length=length(spp_ordered_by_peak))
  TOT_high <- rep(NA,length=length(spp_ordered_by_peak))
  #empty dfs for low var at high freq
  AUC_amount_low <- rep(NA,length=length(spp_ordered_by_peak))
  AUC_fraction_low <- rep(NA,length=length(spp_ordered_by_peak))
  TOT_low <- rep(NA,length=length(spp_ordered_by_peak))
  
  AUC_totalamt <- rep(NA,length=length(spp_ordered_by_peak))
  
  for (s in 1:length(spp_ordered_by_peak)){ #for each spp...
    
    # ** Check inequality sign for percent high or low variability **
    # multiply ht of each bar under curve by the width (one unit of x, freq)
    AUC_amount_high[s] <- sum(ts.spec[ts.spec$variable == "Nsize" 
                                  & ts.spec$EI==0.72
                                  & ts.spec$spp==spp_ordered_by_peak[s]
                                  & ts.spec$noise==noisenames[n]
                                  & ts.spec$freq > AUCthreshold_ordered_by_peak[s],]$spec)
                               # > is for amount of high, =< is for amount of low
    AUC_amount_low[s] <- sum(ts.spec[ts.spec$variable == "Nsize" 
                                  & ts.spec$EI==0.72
                                  & ts.spec$spp==spp_ordered_by_peak[s]
                                  & ts.spec$noise==noisenames[n]
                                  & ts.spec$freq <= AUCthreshold_ordered_by_peak[s],]$spec)
    AUC_fraction_high[s] <- sum(ts.spec[ts.spec$variable == "Nsize" 
                                  & ts.spec$EI==0.72
                                  & ts.spec$spp==spp_ordered_by_peak[s]
                                  & ts.spec$noise==noisenames[n]
                                  & ts.spec$freq > AUCthreshold_ordered_by_peak[s],]$specN)
    
    AUC_fraction_low[s] <- sum(ts.spec[ts.spec$variable == "Nsize" 
                                  & ts.spec$EI==0.72
                                  & ts.spec$spp==spp_ordered_by_peak[s]
                                  & ts.spec$noise==noisenames[n]
                                  & ts.spec$freq <= AUCthreshold_ordered_by_peak[s],]$specN)
    
    AUC_totalamt[s] <- sum(ts.spec[ts.spec$variable == "Nsize" 
                                  & ts.spec$EI==0.72
                                  & ts.spec$spp==spp_ordered_by_peak[s]
                                  & ts.spec$noise==noisenames[n],]$spec)
    
    TOT_high[s] <- (totalvar[totalvar$spp==spp_ordered_by_peak[s] 
                             & totalvar$EI==0.72 
                             & totalvar$variable=="Nsize" 
                             & totalvar$noise==noisenames[n],]$variance)*AUC_fraction_high[s]
    
    TOT_low[s] <- (totalvar[totalvar$spp==spp_ordered_by_peak[s]
                            & totalvar$EI==0.72 
                            & totalvar$variable=="Nsize"
                            & totalvar$noise==noisenames[n],]$variance)*AUC_fraction_low[s]
  }
  #store percents for each type of noise
  AUC_amount_highL[[n]] <- AUC_amount_high 
  AUC_amount_lowL[[n]] <-  AUC_amount_low 
  AUC_fraction_highL[[n]] <- AUC_fraction_high 
  AUC_fraction_lowL[[n]] <- AUC_fraction_low
  TOT_highL[[n]] <- TOT_high
  TOT_lowL[[n]] <- TOT_low
  AUC_totalamt_L[[n]] <- AUC_totalamt
  
  print(noisenames[n])
}
rm(n,s)
AUC_amount_highdf <- data.frame(do.call(cbind,AUC_amount_highL))
AUC_amount_lowdf <- data.frame(do.call(cbind,AUC_amount_lowL))
AUC_fraction_highdf <- data.frame(do.call(cbind,AUC_fraction_highL))
AUC_fraction_lowdf <- data.frame(do.call(cbind,AUC_fraction_lowL))
TOT_highdf <- data.frame(do.call(cbind,TOT_highL))
TOT_lowdf <- data.frame(do.call(cbind,TOT_lowL))
AUC_totalamt_df <- data.frame(do.call(cbind,AUC_totalamt_L))

AUC_amount_highdf$spp <- spp_ordered_by_peak
AUC_amount_lowdf$spp <- spp_ordered_by_peak
AUC_fraction_highdf$spp <- spp_ordered_by_peak
AUC_fraction_lowdf$spp <- spp_ordered_by_peak
AUC_totalamt_df$spp <- spp_ordered_by_peak
TOT_highdf$spp <-spp_ordered_by_peak
TOT_lowdf$spp <-spp_ordered_by_peak

# convert dfs to long format
AUC_amount_highdflong <- AUC_amount_highdf %>% 
  gather(noise,value,1:length(noisenames)) %>% 
  mutate(AUCdes=rep("amt_high_spec"))

AUC_amount_lowdflong <- AUC_amount_lowdf %>% 
  gather(noise,value,1:length(noisenames)) %>% 
  mutate(AUCdes=rep("amt_low_spec") )

AUC_fraction_highdflong <- AUC_fraction_highdf %>% 
 gather(noise,value,1:length(noisenames)) %>% 
  mutate(AUCdes=rep("fra_high"))

AUC_fraction_lowdflong <- AUC_fraction_lowdf %>% 
 gather(noise,value,1:length(noisenames)) %>% 
  mutate(AUCdes=rep("fra_low"))

AUC_total_dflong <- AUC_totalamt_df %>% 
  gather(noise,value,1:length(noisenames)) %>% 
  mutate(AUCdes=rep("total_spec"))

TOT_high_dflong <- TOT_highdf %>%
  gather(noise,value,1:length(noisenames)) %>%
  mutate(AUCdes=rep("tot_high_sgm"))

TOT_low_dflong <- TOT_lowdf %>%
  gather(noise,value,1:length(noisenames)) %>%
  mutate(AUCdes=rep("tot_low_sgm"))

AUCdat <- rbind(AUC_amount_highdflong,
                AUC_amount_lowdflong,
                AUC_fraction_highdflong,
                AUC_fraction_lowdflong,
                AUC_total_dflong,
                TOT_high_dflong,
                TOT_low_dflong)
rm(AUC_amount_high,AUC_amount_highdf,AUC_amount_highdflong,AUC_amount_highL,
   AUC_amount_low,AUC_amount_lowdf,AUC_amount_lowdflong,AUC_amount_lowL,
   AUC_fraction_high,AUC_fraction_highdf,AUC_fraction_highdflong,AUC_fraction_highL,
   AUC_fraction_low,AUC_fraction_lowdf,AUC_fraction_lowdflong,AUC_fraction_lowL,
   AUC_totalamt_df,AUC_total_dflong,
   TOT_high,TOT_highdf,TOT_high_dflong,TOT_highL,
   TOT_low,TOT_lowdf,TOT_low_dflong,TOT_lowL)
# add age structure columns to AUCdat
AUCdat$peak <- spawndistmetrics[match(AUCdat$spp,spawndistmetrics$spp),"mode_age"]
AUCdat$maxage <- spawndistmetrics[match(AUCdat$spp,spawndistmetrics$spp),"maxage"]
AUCdat$cvs <- spawndistmetrics[match(AUCdat$spp,spawndistmetrics$spp),"cvs_mode"]

# # # -------------------------------------------------------
# # # Format data to make plots
# # 
# # # a) Is total CV different with ENSO? Higher or lower?
# # #    Plot: amt and fraction low fraction as a function of totalvar
# # df.tot <- totalvar %>% filter(noise %in% c("white","ENSO","ENSOfast","ENSOslow")) %>%
# #   filter(EI==0) %>% filter(variable=="Nsize")
# # df.amt <- AUCdat %>% filter(noise %in% c("white","ENSO","ENSOfast","ENSOslow")) %>%
# #   filter(AUCdes=="tot_low_sgm") %>% select(spp,noise,value)
# # names(df.amt)[names(df.amt)=="value"]<-"tot_low_sgm"
# # df.fra <- AUCdat %>% filter(noise %in% c("white","ENSO","ENSOfast","ENSOslow")) %>%
# #   filter(AUCdes=="fra_low") %>% select(spp,noise,value)
# # names(df.fra)[names(df.fra)=="value"]<-"fra_low"
# # df.splo <- AUCdat %>% filter(noise %in% c("white","ENSO","ENSOfast","ENSOslow")) %>%
# #   filter(AUCdes=="amt_low_spec") %>% select(spp,noise,value)
# # names(df.splo)[names(df.splo)=="value"]<-"amt_low_spec"
# # df.totAUC <- AUCdat %>% filter(noise %in% c("white","ENSO","ENSOfast","ENSOslow")) %>%
# #   filter(AUCdes=="total_spec") %>% select(spp,noise,value)
# # names(df.totAUC)[names(df.totAUC)=="value"]<-"total_spec"
# # df.combo <- full_join(df.tot,df.amt)
# # df.combo <- full_join(df.combo,df.fra)
# # df.combo <- full_join(df.combo,df.splo)
# # df.combo <- full_join(df.combo,df.totAUC)
# # 
# # 
# # p1 <- ggplot(data=df.combo,aes(x=noise,y=total_spec,shape=spp,color=maxage)) +
# #   geom_point() + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(data=df.combo,aes(x=noise,y=total_spec,group=spp,color=maxage)) +
# #   xlab("Noise") +
# #   ylab("Total AUC\n(in SSB spectrum)") + theme(legend.position="none")
# # 
# # p2 <- ggplot(data=df.combo) +
# #   geom_point(aes(x=noise,y=variance,shape=spp,color=maxage)) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(data=df.combo,aes(x=noise,y=variance,group=spp,color=maxage)) +
# #   ylab("Total variance\n(sigma/CV in SSB time series)") +
# #   xlab("Noise") + theme(legend.position="none")
# # 
# # p3 <- ggplot(data=df.combo) +
# #   geom_point(aes(x=noise,y=fra_low,shape=spp,color=maxage)) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(data=df.combo,aes(x=noise,y=fra_low,group=spp,color=maxage)) +
# #   ylab("Fraction of AUC at low freq\n(in SSB spectrum)") +
# #   xlab("Noise") + theme(legend.position="none")
# # 
# # p4 <- ggplot(data=df.combo) +
# #   geom_point(aes(x=noise,y=tot_low_sgm,shape=spp,color=maxage)) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(aes(x=noise,y=tot_low_sgm,group=spp,color=maxage)) +
# #   ylab("Amount of variance at low freq \n(sigma/mean)*(fraction of AUC at low freq)") +
# #   xlab("Noise") + theme(legend.position="none")
# # 
# # p5 <- ggplot(data=df.combo) +
# #   geom_point(aes(x=noise,y=amt_low_spec,shape=spp,color=maxage)) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(aes(x=noise,y=amt_low_spec,group=spp,color=maxage)) +
# #   ylab("Amount of variance at low freq \n(AUC at low freq in SSB spectrum)") +
# #   xlab("Noise")
# # 
# # 
# # 
# # tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/variance_change_with_noise4.tiff', units="in", width=10, height=6, res=300) 
# # p <- list(p1,p2,p3,p4)
# # do.call(grid.arrange,c(p,ncol=2))
# # dev.off()
# # tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/variance_change_with_noise1.tiff', units="in", width=8, height=5, res=300) 
# #  ggplot(data=df.combo) +
# #   geom_point(aes(x=noise,y=amt_low_spec,shape=spp,color=maxage)) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(aes(x=noise,y=amt_low_spec,group=spp,color=maxage)) +
# #   ylab("Amount of variance at low freq \n(AUC at low freq in SSB spectrum)") +
# #   xlab("Noise")+ theme(legend.direction = "vertical", legend.box = "horizontal")
# # dev.off()
# # 
# # 
# # 
# # # How does the AUC at low freq & fraction
# # # of AUC at low frequencies change when we go
# # # from white to ENSO?
# # p6 <- ggplot(data=df.combo) +
# #   geom_point(aes(x=total_spec,y=fra_low,shape=spp,color=noise),size=2) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(data=df.combo,aes(x=total_spec,y=fra_low,group=spp)) +
# #   ylab("Fraction of AUC at low freq") +
# #   xlab("Total AUC") + theme(legend.position="none")
# # 
# # p7 <- ggplot(data=df.combo) +
# #   geom_point(aes(x=total_spec,y=amt_low_spec,shape=spp,color=noise),size=2) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(data=df.combo,aes(x=total_spec,y=amt_low_spec,group=spp)) +
# #   ylab("Amount of AUC at low freq") +
# #   xlab("Total AUC") + theme(legend.position="none")
# # 
# # p8 <- ggplot(data=df.combo) +
# #   geom_point(aes(x=variance,y=tot_low_sgm,shape=spp,color=noise),size=2) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(data=df.combo,aes(x=variance,y=tot_low_sgm,group=spp)) +
# #   ylab("Amount of variance at low freq \n(sigma/mean)*(fraction of AUC at low freq)") +
# #   xlab("Total variance (sigma/mean)")+ theme(legend.direction = "vertical", legend.box = "horizontal")
# # 
# # tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/lowfreqvar_as_a_function_of_totalvar2.tiff', units="in", width=7, height=4, res=300) 
# # p <- list(p6,p7)
# # do.call(grid.arrange,c(p,ncol=2))
# # dev.off()
# # 
# # tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/lowfreqvar_as_a_function_of_totalvar1.tiff', units="in", width=7, height=5, res=300) 
# # ggplot(data=df.combo) +
# #   geom_point(aes(x=variance,y=tot_low_sgm,shape=spp,color=noise),size=2) + 
# #   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
# #   geom_line(data=df.combo,aes(x=variance,y=tot_low_sgm,group=spp)) +
# #   ylab("Amount of variance at low freq \n(sigma/mean)*(fraction of AUC at low freq)") +
# #   xlab("Total variance (sigma/mean)") + theme(legend.direction = "vertical", legend.box = "horizontal")
# # dev.off()
# # # How does pOFL depend on variability?
# # # Merge freq-dependent var & pOFL data frames (for plots)
# # # match on spp, noise (note: EI=0 for all right now)
# probdf.N <- probdf %>% filter(variable=="Nsize") %>% filter(noise %in% c("white","ENSO","ENSOfast","ENSOslow")) %>% select(spp,noise,probs)
# df.combo.p <- full_join(df.combo,probdf.N)
# head(df.combo.p)
# 
# p9 <- ggplot(data=df.combo.p) +
#   geom_point(aes(x=total_spec,y=probs,shape=spp,color=noise),size=2) + 
#   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
#   geom_line(data=df.combo.p,aes(x=total_spec,y=probs,group=spp)) +
#   ylab("Fraction of time below OFL") +
#   xlab("Total AUC\n(SSB spectrum)") + theme(legend.position="none")
# 
# p10 <- ggplot(data=df.combo.p) +
#   geom_point(aes(x=variance,y=probs,shape=spp,color=noise),size=2) + 
#   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
#   geom_line(data=df.combo.p,aes(x=variance,y=probs,group=spp)) +
#   ylab("Fraction of time below OFL") +
#   xlab("Total variance\n(sigma/mean in SSB)") + theme(legend.position="none")
# 
# p11 <- ggplot(data=df.combo.p) +
#   geom_point(aes(x=amt_low_spec,y=probs,shape=spp,color=noise),size=2) + 
#   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
#   geom_line(data=df.combo.p,aes(x=amt_low_spec,y=probs,group=spp)) +
#   ylab("Fraction of time below OFL") +
#   xlab("Amount of AUC at low freq\n(in SSB spectrum)") + theme(legend.position="none")
# 
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/pOFL_dependence_on_lowfreqvariance.tiff', units="in", width=8, height=6, res=300) 
# ggplot(data=df.combo.p) +
#   geom_point(aes(x=tot_low_sgm,y=probs,shape=spp,color=noise),size=2) + 
#   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
#   geom_line(data=df.combo.p,aes(x=tot_low_sgm,y=probs,group=spp)) +
#   ylab("Fraction of time below OFL") + 
#   xlab("Amount of variance at low freq \n(sigma/mean)*(fraction of AUC at low freq)") 
#  
# dev.off()
# 
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/pOFL_dependence_on_fractionlowfreq.tiff', units="in", width=8, height=5, res=300) 
# ggplot(data=df.combo.p) +
#   geom_point(aes(x=fra_low,y=probs,shape=spp,color=noise),size=3) + 
#   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
#   geom_line(data=df.combo.p,aes(x=fra_low,y=probs,group=spp,color=""),color="gray90") +
#   ylab("Fraction of time below OFL") +
#   xlab("Fraction of AUC at low freq\n(in SSB spectrum)") + theme(legend.direction = "vertical", legend.box = "horizontal") + theme_classic() +
#   geom_smooth(aes(x=fra_low,y=probs),color="gray50",method="lm",se=FALSE)
# dev.off()
# 
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/pOFL_dependence_on_variance6.tiff', units="in", width=7, height=8, res=300) 
# p <- list(p9,p10,p11,p12)
# do.call(grid.arrange,c(p,ncol=2))
# dev.off()
# 
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/pOFL_dependence_on_variance1.tiff', units="in", width=6, height=4, res=300) 
# ggplot(data=df.combo.p) +
#   geom_point(aes(x=fra_low,y=probs,shape=spp,color=noise),size=2) + 
#   scale_shape_manual(values=seq(from=0,to=length(unique(totalvar$spp)))) +
#   geom_line(data=df.combo.p,aes(x=fra_low,y=probs,group=spp)) +
#   ylab("Fraction of time below OFL") +
#   xlab("Fraction of AUC at low freq\n(in SSB spectrum)") + theme(legend.direction = "vertical", legend.box = "horizontal")
# dev.off()
# 
# # Plot: difference in pOFL between white and ENSO vs diff in fraction of low frequency variability
# diff.df <- data.frame(cbind(as.character(unique(probdf.N$spp)),
#                  rep(NA,length(unique(probdf.N$spp))),
#                  rep(NA,length(unique(probdf.N$spp))),
#                  rep(NA,length(unique(probdf.N$spp)))),stringsAsFactors = FALSE)
# colnames(diff.df) <- c("spp","noise_pair","lowfreqvar_diff","pofl_diff")
# for(s in 1:length(unique(df.combo.p$spp))){
#   t <- df.combo.p[df.combo.p$spp==unique(df.combo.p$spp)[s],]
#   # subtract low var w/white noise from low var w/ENSO
#   diff.df[diff.df$spp==unique(probdf.N$spp)[s],]$lowfreqvar_diff <- t[t$noise=="white",]$amt_low_spec - t[t$noise=="ENSO",]$amt_low_spec
#   # subtract pofl w/white noise from pofl w/ENSO
#   diff.df[diff.df$spp==unique(probdf.N$spp)[s],]$pofl_diff <- t[t$noise=="white",]$probs - t[t$noise=="ENSO",]$probs
#   # fill in noise-pair
#   diff.df[diff.df$spp==unique(probdf.N$spp)[s],]$noise_pair <- "white-enso" 
# }
# str(diff.df)
# diff.df$lowfreqvar_diff <- as.numeric(diff.df$lowfreqvar_diff)
# diff.df$pofl_diff <- as.numeric(diff.df$pofl_diff)
# diff.df$maxage <- rep(NA,length(diff.df[,1]))
# diff.df$maxage <- spawndistmetrics[match(diff.df$spp,spawndistmetrics$spp),"maxage"]
# 
# p14 <- ggplot(data=diff.df) + 
#   geom_point(aes(x=lowfreqvar_diff,y=pofl_diff,shape=spp,color=maxage),size=2) +
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   xlab("ENSO-white (amount AUC at low freq)") + ylab("ENSO-white (pOFL)") +
#   #geom_smooth(aes(x=lowfreqvar_diff,y=pofl_diff),method='lm',formula=y~x) +
#   ggtitle("Difference in pOFL between\nwhite and ENSO as a function\nof the difference in amount\nof AUC at low freq between\nwhite and ENSO") +
#   scale_shape_manual(name="Species",values=seq(from=0,to=length(unique(df$spp)))) 
# 
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/pOFL_depends_on_noise.tiff', units="in", width=6, height=4, res=300) 
# ggplot(data=df.combo.p) + 
#   geom_point(aes(x=noise,y=probs,shape=spp,color=maxage),size=2) +
#   xlab("ENSO-white (amount AUC at low freq)") + ylab("ENSO-white (pOFL)") +
#   geom_line(aes(x=noise,y=probs,group=spp,color=maxage)) +
#   scale_shape_manual(name="Species",values=seq(from=0,to=length(unique(df$spp)))) 
# dev.off()
# 
# head(ts.spec)
# head(ts.data)
# ggplot(data=ts.spec[ts.spec$spp=="Sardine" &
#                     ts.spec$noise %in% c("white","ENSO") &
#                     ts.spec$variable=="Nsize" &
#                     ts.spec$EI==0.72,]) + 
#   geom_line(aes(x=freq,y=specN,color=noise),size=1) +
#   xlab("Frequency") + ylab("Variance") +
#   theme_bw() + scale_y_log10()
# 
# ggplot(data=ts.data[ts.data$spp=="Sardine" &
#                     ts.data$noise %in% c("white","ENSO") &
#                     ts.data$variable=="Nsize" &
#                     ts.data$EI==0.72,]) + 
#   geom_line(aes(x=year,y=value,color=noise),size=1) +
#   xlab("Frequency") + ylab("Variance") +
#   theme_bw() + scale_y_log10()
# 
# ggplot(data=ts.data[ts.data$spp=="Sardine" &
#         ts.data$FLEP==0.35 &
#         ts.data$variable=="Nsize" &
#         ts.data$noise %in% c("white","ENSO") &
#         ts.data$year %in% rm_first_timesteps:(timesteps-2900),]) +
#   geom_line(aes(x=year,y=value,linetype=noise)) +
#   scale_y_log10() + 
#   geom_hline(data=ts.data[ts.data$spp=="Sardine" &
#         ts.data$FLEP==0.30 &
#         ts.data$variable=="Nsize" &
#         ts.data$noise %in% c("white","ENSO") &
#         ts.data$year %in% rm_first_timesteps:(timesteps-2900),],
#              aes(yintercept = means,linetype=noise)) +
#   scale_color_manual(values=c("black","grey")) +
#   theme_classic() +
#   ggtitle("Sardine") +
#   theme(plot.title = element_text(margin = margin(t = 10, b = -20))) 
# 
#  
#     tempthres <-ts.data[ts.data$spp==spp[s] &
#                           ts.data$FLEP==0.3 &
#                           ts.data$variable=="Nsize" &
#                           ts.data$noise %in% c("white") &
#                           #ts.data$noise==noisenames[n] &
#                           ts.data$year %in% rm_first_timesteps:(timesteps-2900),] 
#     
#     multiplots[[s]] <- ggplot(data=tempts,aes(x=year,y=value,color=noise)) +
#       geom_line() +
#       scale_y_log10() +
#       #scale_y_log10(limits=c(minSSB,maxSSB)) +
#      
#       ggtitle(paste(noisenames[n],"\n",spp[s],"(",spawndistmetrics[spawndistmetrics$spp==spp[s],]$mode_age,")","\nOFL=mean SSB when FLEP=0.3",", FLEP=0.35",sep="")) +
#       ylab("Spawning stock biomass") +
#      #facet_grid(noise ~ variable) +
#      #theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
#      #theme(axis.title.y = element_text(angle = 0)) +
#      theme(legend.position="none",
#            plot.title = element_text(margin = margin(t = 10, b = -20))) +
#       geom_hline(data=tempthres,
#              aes(yintercept = means)) 
#     
# # Plot: show fraction of low freq variance vs max age and peak age
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/fig1g_fralow_vs_maxage.tiff', units="in", width=6, height=4, res=300) 
# ggplot(data=df.combo[df.combo$EI==0 & 
#                      df.combo$variable=="Nsize" &
#                      df.combo$noise=="white",],aes(x=maxage,y=fra_low)) +
#   geom_point() + xlab("Maximum Age (yrs)") + ylab("Fraction of total variance\nin low frequencies") +
#   geom_label_repel(aes(label=spp),label.size=NA,) +
#   theme_classic()
# dev.off()
#   
# # Plot barplot: for each species (maybe only 5?) bar for pOFL in each noise
# tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/fig5_barplot.tiff', units="in", width=8, height=4, res=300) 
# noise_levels <- c("white","ENSO","ENSOfast","ENSOslow")
# df.combo.p$noise <- factor(df.combo.p$noise,levels=noise_levels)
# spawndistmetrics$spp_max <- paste(spawndistmetrics$spp," (",spawndistmetrics$maxage,")",sep="")
# spp_levels_maxage <- spawndistmetrics[order(spawndistmetrics$maxage),]$spp_max 
# df.combo.p$spp_max <- paste(df.combo.p$spp," (",df.combo.p$maxage,")",sep="")
# df.combo.p$spp_max <- factor(df.combo.p$spp_max, levels=spp_levels_maxage)
# ggplot(data=df.combo.p[df.combo.p$variable=="Nsize",],
#        aes(x=spp_max,y=probs, fill=noise)) +
#   geom_bar(stat = "identity",position = position_dodge()) +
#   theme_bw() + theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
#   scale_fill_grey() + ylab("Fraction of time below OFL") + xlab("")
# dev.off()