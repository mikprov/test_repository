# RUN
# by: mikaela provost
# last edited: Aug 1, 2019
# ===================================================================

library(tidyr)
library(dplyr)
library(gridExtra)
library(ggplot2)

# *************************************** #
# Plan:
# (1) run script to calculate life table information and load parms df
# (2) assemble Leslie matrices for different F values
# (3) simulate time series 
# (4) calculate spectra



# *************************************** #
# (1) run script to get life table information, parms df, and functions
source("C:/Users/provo/Documents/GitHub/pfmc/src/life_table_calculations.r")

# *************************************** #
# (2) assemble Leslie matrices for different F values

# make all LEPs equal
conLEP = 1.1
# adjust fecundities by this much
adjFec = round(1/conLEP,digits=1)


Fvals = seq(from=0,to=3,by=0.1) #instantaneous F rates
A3d_list = as.list(rep(NA,length(parms$spp))) #Leslie arrays in list
LTABLE_list <- as.list(rep(NA,length(parms$spp))) #LTABLE arrays in list
names(A3d_list) <- parms$spp
names(LTABLE_list) <- parms$spp
eigan1 <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
eigan12 <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
LEPs <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))
newLEP <- matrix(NA,nrow=length(Fvals),ncol=length(parms$spp))

for (i in 1:length(parms$spp)){ #for each spp
  
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes") {#if a problem, then
    problemshold <- parms$spp[i]} 
  else {
  # define parms for spp i
  M = as.numeric(parms[parms$spp == parms$spp[i],]$M)
  maxage = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  
  # empty arrays for Leslie and LTABLE, vectors for eigens
  A3d <- array(NA,c(maxage,maxage,length(Fvals)))
  LTABLE3d <- array(NA,c(maxage,12,length(Fvals))) #12 cols
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
                           Fval=Fvals[f])
    LTABLE3d[,,f] <- data.matrix(out$LTABLE) #store LTABLE in 3d array
    # calculate LEP at F
    lep[f] = sum(out$LTABLE$LEP_a)
    
    # adjust fecundities so LEP equal for all species
    A <- out$A
    A[1,] <- (A[1,]/(lep[f]*adjFec))
    A3d[,,f] <- A
    
    # re-calculate new LEP (should equal conLEP)
    term <- rep(NA,length=maxage)
    for(a in 1:maxage){term[a] <- A[1,a]*out$LTABLE$Survship[a]} #fecundity*survivorship
    newLEP[f,i] <- sum(term) #if you plot 'term', it's the spawning distribution over age
    
    # extract lambda1, lambda2/lambda1
    e1[f] = extract_first_eigen_value(out$A)
    e2[f] = extract_second_eigen_value(out$A)
    e12[f] = e2[f] / e1[f]
    
    
  } #closes F loop
  A3d_list[[i]] <- A3d 
  LTABLE_list[[i]] <- LTABLE3d
  eigan1[,i] <- e1
  eigan12[,i] <- e12
  LEPs[,i] <- lep
  #clean up before next i in loop
  #rm(A3d,LTABLE3d,out,M,maxage) 
  } #closes 'else' bracket for non-problem species.
  print(i)
} #closes spp loop
rm(i,f,a,M,maxage,problemshold,lep,term) #clean up


# *************************************** #
# (3) Simulate pops. Loop over Aarray list to simulate using different Leslie matrices 
# *************************************** #
# set params for simulation:
# note: alpha > 1/LEP. Since I set LEP=1.1 for all spp, alpha must be > 1/1.1
timesteps = 5000 #need this now to create
rm_first_timesteps = 2000
alpha = 1.2
betas = 1000
#sig_r = 0.3
span.multiplier = 1 # adjusting the span in spec.prgm()
alphas <- rep(alpha, length=length(parms$spp)) #alpha could be diff for pops

output.3d.list <- as.list(rep(NA,length=length(parms$spp))) #store timeseries here
names(output.3d.list) <- parms$spp

# load simulation model:
source("C:/Users/provo/Documents/GitHub/pfmc/src/simulation_model.r")

for (i in 1:length(A3d_list)) { #step through each pop
  
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes") {#if a problem, then
    problemshold <- parms$spp[i]} 
  else {
  
    Leslie3d = A3d_list[[i]] #select the 3d array of Leslie matricies
    # array dims: row=ts length, col=4 is number of ts (eggs,recruits,Nt,Nsize), depth=F vals
    output.matrix <- array(NA,c(timesteps-1,3,length(Fvals))) 
  
    for (f in 1:length(Fvals)) { #step through each Leslie matrix (for each F value)
    
      output = sim_model(A=Leslie3d[,,f], timesteps=timesteps, 
                       alpha=alphas[i], beta=betas, 
                       sig_r=as.numeric(parms$sigmaR[i]), initial_eggs=betas)
    
      length(output$Nsize) <- length(output$eggs) #trim Nsize ts vector, -1 elements
      output.matrix[,,f] <- do.call(cbind,output) #fill in array for pop i
      #colnames(output.matrix) <- names(output)
  }
  
  output.3d.list[[i]] <- output.matrix
  }#closes else bracket
  print(i)
}
rm(i,f,output.matrix,output,problemshold) #clean up

# At this point I have one important object:
# 1. [output.3d.list] a list of 3d arrays. Each array is timeseries output
#    from simulations at different F levels. 


# *************************************** #
# (5) Format output ts for plotting simulations using output.3d.list
# *************************************** #
variable_type <- c("eggs","recruits","Nsize")
# --- reorganize recruit timeseries data --- #
var.number <- which(variable_type == "recruits") # recruits
df.list <- as.list(rep(NA,length=length(parms$spp)))
names(df.list) <- parms$spp
for (i in 1:length(output.3d.list)) {
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes") {#if a problem, then
    problemshold <- parms$spp[i]} 
  else { #if not problem spp, reformat data...
  # ...to work with ggplot
  aa <- as.data.frame(output.3d.list[[i]][,var.number,])
  aa$year <- seq(from=1, to=length(aa[,1]),by=1)
  colnames(aa) <- c(Fvals,"year")
  aa1 <- aa %>% gather(Fvalue,value,1:length(Fvals))
  aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
  aa1$spp <- rep(parms$spp[i],length=length(aa[,1]))
  aa1$Fvalue <- as.numeric(as.character(aa1$Fvalue))
  df.list[[i]] <- aa1
  }
  print(i)
}
# keep only df in list for non-problem species
df.list <- df.list[parms[parms$problem_spp == "no",]$spp]
parms <- parms[parms$problem_spp == "no",]
# do some formatting on new parms
parms$maxage <- as.numeric(parms$maxage)

recruits.ts <- bind_rows(df.list,id=NULL)
str(recruits.ts)
rm(df.list,i,aa1,aa,var.number) #clean up
ts.data <- recruits.ts
ts.data$maxage <- parms[match(ts.data$spp,parms$spp),"maxage"]

# *************************************** #
# (6) Now that timeseries data is formated, let's plot! 
# *************************************** #
#plot recruitment - one plot per spp
str(ts.data)
#vector of spp names in order of max age
spp.ordered.maxage <- parms %>% arrange(maxage) %>% pull(spp)
p.forappendix <- as.list(rep(NA,length=length(spp.ordered.maxage)))
names(p.forappendix) <- spp.ordered.maxage

for (i in 1:length(spp.ordered.maxage)){
  dd <- ts.data[ts.data$variable == "recruits" & 
                  ts.data$spp == spp.ordered.maxage[i] &
                  ts.data$year %in% seq(from=rm_first_timesteps,to=(timesteps-100),by=1) &
                  ts.data$Fvalue %in% c(1),]
  dd$Fvalue <- as.factor(dd$Fvalue)
  #subtract mean from ts, and adjust so all values are positive (makes it easier to interpet)
  #dd$valueminusmean <- (dd$value-mean(dd$value))+500 
  
  p.forappendix[[i]] <- ggplot(dd,aes(x=year,y=value,color=Fvalue)) +
    xlab("") + ylab("") +
    geom_line() + theme_classic() +
    #scale_color_brewer(palette = "Reds") +
    ylim(c(150,1000)) +
    ggtitle(paste(spp.ordered.maxage[i]," maxage=",dd$maxage[1],sep="")) 
    print(i)
}
tiff(file='C:/Users/provo/Documents/GitHub/pfmc/results/timeseries_v1_F1.tiff', units="in", width=7, height=13, res=300) 
do.call(grid.arrange,c(p.forappendix,ncol=2,left="Recruits (before noise)", bottom="Year",
                       top="F=1,"))
dev.off()




# *************************************** #
# (6) Calculate frequency content from timeseries
# *************************************** #
# Plan:
# 1. Walk through each cod pop, do spectral analysis at beta levels
# 2. Store spec values for eggs, recruits, and Nsize

# 1. Walk through each cod pop, do spectral analysis at alpha levels
#sp.eggsL <- as.list(rep(NA,length=length(codNames))) #object for spec analysis  
sp.recruitLsm <- as.list(rep(NA,length=length(spp.ordered.maxage))) 
span.multiplier= 1
Fvalue = 1

for (i in 1:length(spp.ordered.maxage)){
  
  ts <- ts.data[ts.data$spp == spp.ordered.maxage[i],] #subset data for pop i 
  
  # setting 'span' - a vector of odd integers to specify the smoothers
  tmp <- ceiling(sqrt(length(1:(timesteps-rm_first_timesteps-1)))) #sq root of timeseries lgth, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  
  # --- spectral analysis on RECRUIT --- #
  spsaveL <- as.list(rep(NA,length=length(Fvalue)))
  names(spsaveL) <- Fvalue
  
  for (f in 1:length(Fvalue)){
    yy = ts %>% filter(Fvalue == 0) %>% select(value) %>% slice(rm_first_timesteps:(timesteps-2)) 
    #meanyy = mean(rawyy$value)
    #yy = rawyy$value-meanyy
    
    #sp = spec.pgram(yy$value,spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
    sp = spec.pgram(yy$value,spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
    spsaveL[[f]] = 2*sp$spec # 
    print(f)
  }
  spsave <- as.data.frame(do.call(cbind,spsaveL))
  spsave$freq <- sp$freq
  spsavelong <- spsave %>% gather(Fvalue, value, 1:length(Fvalue))
  spsavelong$spp.ordered.maxage <- rep(spp.ordered.maxage[i],length=length(spsavelong[,1]))
  sp.recruitLsm[[i]] <- spsavelong
  rm(spsave,yy,sp)
  print(i)
}

recsm <- bind_rows(sp.recruitLsm,id=NULL)
recsm$variable.type <- rep("recruits",length=length(recsm$freq))
specdatalong <- recsm
head(specdatalong)
specdatalong$maxage <- parms[match(specdatalong$spp.ordered.maxage,parms$spp),"maxage"]
str(specdatalong)
rm(recsm)

# *************************************** #
# (6) Plot spectral analysis 
# *************************************** #
# --- recruits spectra --- #
prec_sp <- list()

for (i in 1:length(spp.ordered.maxage)){
  dataforplot <- specdatalong[specdatalong$variable.type == "recruits" 
                              & specdatalong$spp.ordered.maxage==spp.ordered.maxage[i],]
  
  # store plots in list
  prec_sp[[i]] <- ggplot(data=dataforplot, aes(x=freq,y=value)) + 
    geom_line() + 
    scale_y_log10() +
    #ylim(c(-3,7)) +
    ggtitle(paste(spp.ordered.maxage[i],
                  " maxage=",dataforplot$maxage[1],
                  sep="")) + 
    theme_classic() + ylab("") + xlab("") #+ scale_y_log10(limits=c(0.01,56000)) 
  
}
names(prec_sp) <- spp.ordered.maxage
rm(i,dataforplot)

tiff(file='C:/Users/provo/Documents/GitHub/pfmc/results/spectra_v1_F1_sm55.tiff', units="in", width=7, height=13, res=300) 
do.call(grid.arrange,c(prec_sp,ncol=2,left="Variance",bottom="Frequency",top="F=0, smoother=55yr"))
dev.off()

tiff(file='C:/Users/provo/Documents/GitHub/pfmc/results/spectra_v1_oneplot_sm55.tiff', units="in", width=5, height=7, res=300)
# Plot spectra on one plot:
ggplot(data=specdatalong, aes(x=freq,y=value,group=spp.ordered.maxage,color=maxage)) + 
  geom_line() + 
  scale_y_log10() +
  theme_classic() + ylab("") + xlab("") +
  ggtitle("smoother=55")
dev.off()
tiff(file='C:/Users/provo/Documents/GitHub/pfmc/results/spectra_v2_oneplot_sm55.tiff', units="in", width=6, height=7, res=300)
ggplot(data=specdatalong, aes(x=freq,y=value,color=spp.ordered.maxage)) + 
  geom_line() + 
  scale_y_log10() +
  theme_classic() + ylab("") + xlab("") +
  ggtitle("smoother=55")
dev.off()
# *************************************** #
# (10) Area under curve
# *************************************** #
head(specdatalong)
# Plan:
# 1. calculate AUC at high and low frequencies (threshold=peak spawning age, others)
# 2. plot %AUClow and %AUChigh vs peak spawning age 
# 3. plot %AUClow and %AUChigh vs max age
# 4. plot %AUClow and %AUChigh vs spawning distribution sd & CV
AUC_greater <- rep(NA,length=length(codNames_ordered_by_peak))
AUC_less <- rep(NA,length=length(codNames_ordered_by_peak))
AUC_total <- rep(NA,length=length(codNames_ordered_by_peak))
AUCperlow <- rep(NA,length=length(codNames_ordered_by_peak))
AUCperhigh <- rep(NA,length=length(codNames_ordered_by_peak))

AUCthreshold <- c(0.05,0.1,0.2,1)
plottitles <- c("threshold freq=0.5","threshold freq=0.1","threshold freq=0.2","threshold freq=peak spawning age")

par(mfrow=c(4,3))
for (j in 1:length(AUCthreshold)){ #step through different AUC thresholds
  
for (i in 1:length(codNames_ordered_by_peak)){
  if(AUCthreshold[j] == 1) 
  {AUCthreshold[j] <- (1/eigentable[eigentable$codNames == codNames_ordered_by_peak[i],]$mode_age)} else
  {AUCthreshold[j]}
  
  AUC_greater[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                     & specdatalong$codNames==codNames_ordered_by_peak[i]
                                     & specdatalong$freq > AUCthreshold[j],]$value)
  
  AUC_less[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                  & specdatalong$codNames==codNames_ordered_by_peak[i]
                                  & specdatalong$freq <= AUCthreshold[j],]$value)
  if(j==3){AUC_less[i] <- AUC_less[i]*AUCthreshold[j]} else {AUC_less[i]}
  if(j==3){AUC_greater[i] <- AUC_greater[i]*(0.5-AUCthreshold[j])} else {AUC_greater[i]}
  
  AUC_total[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                   & specdatalong$codNames==codNames_ordered_by_peak[i],]$value)
  AUCperlow[i] <- AUC_less[i]/AUC_total[i]
  AUCperhigh[i] <- AUC_greater[i]/AUC_total[i]
  print(AUCperlow[i]+AUCperhigh[i])
}
plot(x=eigentable[order(eigentable$mode_age),]$mode, 
               xlab="peak spawning age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
plot(x=eigentable[order(eigentable$mode_age),]$max_ages, 
               xlab="max age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
plot(x=eigentable[order(eigentable$mode_age),]$cvs_mode, 
              xlab="Spawning biomass distribution CV",y=AUCperlow,
             ylab="Percent AUC at low freq",main=plottitles[j])
#plist <- as.list(p1peak,p2max,p3CV)
#pList[[j]] <- plist

}
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(x=eigentable[order(eigentable$mode_age),]$mode, 
     xlab="peak spawning age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable[order(eigentable$mode_age),]$max_ages, 
     xlab="max age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable[order(eigentable$mode_age),]$cvs_mode, 
     xlab="Spawning biomass distribution CV",y=AUC_total,
     ylab="Total AUC")
par(mfrow=c(1,1))

# *************************************** #
# (10b) Area under curve - remove pops with truncated distributions
# *************************************** #
keep <- c("Northsea","W_Baltic","Faroe","Celtic","Iceland","GB","GM",
              "cod2J3KL","cod3Ps")
codNames_ordered_by_peak_keep <- codNames_ordered_by_peak[codNames_ordered_by_peak %in% keep]
eigentable.trunc <- eigentable[eigentable$codNames %in% codNames_ordered_by_peak_keep,]

AUC_greater <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUC_less <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUC_total <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUCperlow <- rep(NA,length=length(codNames_ordered_by_peak_keep))
AUCperhigh <- rep(NA,length=length(codNames_ordered_by_peak_keep))

AUCthreshold <- c(0.1,0.2,1)
plottitles <- c("threshold freq=0.1","threshold freq=0.2","threshold freq=peak spawning age")

par(mfrow=c(4,3))
for (j in 1:length(AUCthreshold)){ #step through different AUC thresholds
  
  for (i in 1:length(codNames_ordered_by_peak_keep)){
    if(AUCthreshold[j] == 1) 
    {AUCthreshold[j] <- (1/eigentable[eigentable$codNames == codNames_ordered_by_peak_keep[i],]$mode_age)} else
    {AUCthreshold[j]}
    
    AUC_greater[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                               & specdatalong$codNames==codNames_ordered_by_peak_keep[i]
                                               & specdatalong$freq > AUCthreshold[j],]$value)
    
    AUC_less[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                            & specdatalong$codNames==codNames_ordered_by_peak_keep[i]
                                            & specdatalong$freq <= AUCthreshold[j],]$value)
    if(j==3){AUC_less[i] <- AUC_less[i]*AUCthreshold[j]} else {AUC_less[i]}
    if(j==3){AUC_greater[i] <- AUC_greater[i]*(0.5-AUCthreshold[j])} else {AUC_greater[i]}
    AUC_total[i] <- sum(freq[1]*specdatalong[specdatalong$variable.type == "recruit" 
                                             & specdatalong$codNames==codNames_ordered_by_peak_keep[i],]$value)
    AUCperlow[i] <- AUC_less[i]/AUC_total[i]
    AUCperhigh[i] <- AUC_greater[i]/AUC_total[i]
    print(AUCperlow[i]+AUCperhigh[i])
  }
  plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$mode, 
       xlab="peak spawning age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
  plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$max_ages, 
       xlab="max age",y=AUCperlow,ylab="Percent AUC at low freq",main=plottitles[j])
  plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$cvs_mode, 
       xlab="Spawning biomass distribution CV",y=AUCperlow,
       ylab="Percent AUC at low freq",main=plottitles[j])
  #plist <- as.list(p1peak,p2max,p3CV)
  #pList[[j]] <- plist
  
}
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$mode, 
     xlab="peak spawning age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$max_ages, 
     xlab="max age",y=AUC_total,ylab="Total AUC")
plot(x=eigentable.trunc[order(eigentable.trunc$mode_age),]$cvs_mode, 
     xlab="Spawning biomass distribution CV",y=AUC_total,
     ylab="Total AUC")
par(mfrow=c(1,1))

