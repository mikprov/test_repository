# PFMC functions
# similar to cod functions

# --------------------
# (1) functions for convert age to length (vonB & Schnute) 
vonbertgrowth <- function(maxage,t0,Linf,Kvonb,ages){
  ages <- seq(from=1,to=maxage,by=1)
  L_a = Linf * (1-exp(-Kvonb*(ages-t0)))
  return(L_a)
  rm(ages)
}

schnutegrowth <- function(maxage,t1,t2,L1,L2,Ksch,ages){
  ages <- seq(from=1,to=maxage,by=1)
  L_a = L1 + (L2-L1) * (1-exp(-Ksch*(ages-t1)))/(1-exp(-Ksch*(t2-t1)))
  return(L_a)
  rm(ages)
}

# --------------------
# (2) function for calculate prop mature at age
calc_mat_at_age <- function(maxage,slope,inflection){
  ages <- seq(from=1,to=maxage,by=1)
  propmat_a <- 1/(1+exp(slope*(ages-inflection)))
  return(propmat_a)
  rm(ages)
}

calc_mat_at_length <- function(lengths_at_age,slope,inflection){
  propmat_l_at_a <- 1/(1+exp(slope*(lengths_at_age-inflection)))
  return(propmat_l_at_a)
}

# --------------------
# (3) function for selectivity at age 
#vul_a

# --------------------
# (4) function for weight at age
calc_wt_at_age <- function(L_a,cmkga,cmkgb){
  B_a <- cmkga*(L_a)^cmkgb
  return(B_a)
}

# --------------------
# (5) functions for fecundity at age (2 versions)
eggsFUN1 <- function(B_a,intercept,slope){
  eggs_a = B_a*(intercept + slope*B_a) 
  return(eggs_a)}
eggsFUN2 <- function(B_a,intercept,slope){
  #if(parms[parms$spp == parms$spp[i],]$eggs_units == "eggs_per_gm") {#if units in gm, then
  #  eggs_a = (intercept + (slope*B_a) )*1000 #convert gm to kg
  #  } else { eggs_a = intercept + slope*B_a }
  eggs_a = intercept + slope*B_a
  return(eggs_a)}

# --------------------
# (6) assemble Leslie
assemble_Leslie <- function(maxage,L_a,B_a,propmat_a,eggs_a,vul_a,M,Fval){
  
  # create vector of ages
  Age=1:maxage 
  
  # assemble LTABLE: 
  LTABLE = data.frame(cbind(Age,L_a,B_a,eggs_a,propmat_a,vul_a))  
  LTABLE$M = rep(M,length=length(LTABLE[,1])) #fill in natural mortality
  LTABLE$Fval = rep(Fval,length=length(LTABLE[,1]))#fill in F mortality
  
  # empty Leslie matrix 
  A = matrix(0,length(Age),length(Age)) 
  
  # for each age, get survival at age from F&M
  for(j in 1:length(Age)){ # step through ages
    LTABLE$FISH[j] = LTABLE$vul_a[j]*Fval	# FISH: F rate = vulnerability-at-age * F  
    LTABLE$SURV_Fing[j] = exp(-(LTABLE$FISH[j]+LTABLE$M[j])) #SURV_Fing:fraction surviving at each age (F+M)
    LTABLE$SURV[j] = exp(-LTABLE$M[j]) #SURV: fraction surviving at each age (M) -> F=0
    
    # set up column for survivorship (amount or fraction present at age, w/ and w/o F)
    LTABLE$Survship_Fing = 0 
    LTABLE$Survship_Fing[1] = 1
    LTABLE$Survship = 0
    LTABLE$Survship[1] = 1}
  
  # for each row in LTABLE (..each age)
  for(k in 1:(nrow(LTABLE)-1)){ # step through ages  
    LTABLE$Survship_Fing[k+1] = LTABLE$Survship_Fing[k]*LTABLE$SURV_Fing[k] #fraction surv from 1 to a w F
    LTABLE$Survship[k+1] = LTABLE$Survship[k]*LTABLE$SURV[k]} #and w/o F
  
  # LEP at age (w/ & w/o F)
  LTABLE$LEP_a <- LTABLE$eggs_a*LTABLE$propmat_a*LTABLE$Survship
  LTABLE$LEP_a_Fing <- LTABLE$eggs_a*LTABLE$propmat_a*LTABLE$Survship_Fing
  
  # insert fecundity (maturity * eggs) in top row
  A[1,] = LTABLE$propmat_a*LTABLE$eggs_a 
  
  # insert survival on subdiagonal (has both M&F)
  for(u in 2:length(Age)-1){ 
    A[u+1,u]=LTABLE$SURV_Fing[u] }
  
  # -- export Leslie matrix (A) & LTABLE df for one F value
  return(list(A=A,LTABLE=LTABLE))
  
}


# --------------------
# (7) extract_first_eigen_value()
extract_first_eigen_value <- function(Lesliematrix){
  # get leading eigenvalue - check to make sure it's positive
  ev = eigen(Lesliematrix)
  # squaring the imaginary part is not necessary b/c the first 
  # eigen value has no imaginary part. 
  a.sq = (Re(ev$values[1]))^2 # square the real part of 1st eigenvalue
  #b.sq = (Im(ev$values[1]))^2 # square the imaginary part of 1st eigenvalue
  firstval = sqrt(a.sq)# + b.sq) # magnitude of 1st eigenvalue
  return(firstval)
  rm(a.sq) # remove from workpace
  #rm(b.sq) # remove from workpace
}

# --------------------
# (8) extract_second_eigen_value()
extract_second_eigen_value <- function(Lesliematrix){
  # get magnitude of second eigenvalue
  # think of real value on x-axis, imaginary value on y-axis
  # magnitude is the vector between them
  ev = eigen(Lesliematrix)
  a.sq = (Re(ev$values[2]))^2 # square the real part of 2nd eigenvalue
  b.sq = (Im(ev$values[2]))^2 # square the imaginary part of 2nd eigenvalue
  secondval = sqrt(a.sq + b.sq) # magnitude of 2nd eigenvalue
  return(secondval)
  rm(a.sq) # remove from workpace
  rm(b.sq) # remove from workpace
}


# --------------------
# (9) calculate LSB at age
calculate_LSB_at_age_by_F <- function(maxage,L_a,B_a,propmat_a,eggs_a,vul_a,M,Fvals){
  
  # test
  #maxage=15
  #L_a= La_list[[i]]
  #B_a= Ba_list[[i]]
  #propmat_a= prop_list[[i]]
  #eggs_a= eggs_list[[i]]
  #vul_a= vul_list[[i]]
  #M=0.4
  #Fvals=c(0,0.1,0.2,0.3)
  
  # create vector of ages
  Age=1:maxage 
  
  # assemble LTABLE: 
  LTABLE = data.frame(cbind(Age,L_a,B_a,eggs_a,propmat_a,vul_a))  
  #LTABLE$M = rep(M,length=length(LTABLE[,1])) #fill in natural mortality
  #LTABLE$Fval = rep(Fval,length=length(LTABLE[,1]))#fill in F mortality
  
  # survival at age for different F values
  surv_at_age_by_F <- matrix(NA,nrow=length(Age),ncol=length(Fvals))
  survivorship_by_F <- matrix(0, nrow=length(Age),ncol=length(Fvals))
  LEP_at_age_by_F <- matrix(0,nrow=length(Age),ncol=length(Fvals))
  for(f in 1:length(Fvals)){ # step through F values
    
    for(a in 1:length(Age)) { #step through ages
    Ftemp = LTABLE$vul_a[a]*Fvals[f]	# FISH: F rate = vulnerability-at-age * F  
    surv_at_age_by_F[a,f] = exp(-(Ftemp+M)) #SURV:fraction surviving at each age
    } #closes loop calculating survival at age
    
    # set up matrix for survivorship (amount or fraction present at age)
    survivorship_by_F[1,] <- 1
   
    # for each age calc survival from age 1 to age a
    for(k in 1:(nrow(survivorship_by_F)-1)){ # step through ages  
      survivorship_by_F[k+1,f] = survivorship_by_F[k,f]*surv_at_age_by_F[k,f]}
    
    # set up matrix for LEP at age
    # LEP at age = prop mat at age * eggs at age * surv from age 0 to age a
    LEP_at_age_by_F[,f] <- LTABLE$eggs_a*LTABLE$propmat_a*survivorship_by_F[,f]
    
    }#closes f loop
  
  # Get FLEP from LEP
  LEP <- colSums(LEP_at_age_by_F)
  FLEP <- LEP/LEP[1]
  FLEP_by_F <- as.data.table(cbind(Fvals,FLEP))
  rm(LEP,FLEP,LTABLE,surv_at_age_by_F,survivorship_by_F,Ftemp,LEP_at_age_by_F,f,a,k)
  return(FLEP_by_F)
}

# --------------------
# (10) calculate F values for exploitation index - using life table vectors
calculate_F_from_EI <- function(EI,maxage,L_a,B_a,propmat_a,eggs_a,vul_a,M,Fvals_many){
  #notes:
  #Fval = should be a sequence of F values (not just one)
  #testing:
  # EI=c(0,0.4,0.6,0.8,0.95)
  # maxage=as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  # L_a=La_list[[i]]
  # B_a=Ba_list[[i]]
  # propmat_a=prop_list[[i]]
  # eggs_a=eggs_list[[i]]
  # vul_a=vul_list[[i]]
  # M=as.numeric(parms[parms$spp == parms$spp[i],]$M)
  # Fvals_many=seq(from=0,to=100,by=0.1)
  
  FLEP_by_F <- calculate_LSB_at_age_by_F(maxage,
                          L_a,
                          B_a,
                          propmat_a,
                          eggs_a,
                          vul_a,
                          M,
                          Fvals)
  #What are corresponding FLEPs for EI values?
  FLEPs_of_interest <- as.data.frame((1-(EI*0.9)))
  colnames(FLEPs_of_interest) <- "FLEP_EI"
  FLEPs_of_interest$matching_F <- rep(NA,length=length(FLEPs_of_interest))
  
  for(e in 1:length(FLEPs_of_interest[,1])){
    diff <- as.data.frame(cbind(FLEP_by_F, (FLEP_by_F$FLEP-FLEPs_of_interest$FLEP_EI[e])))
    FLEPs_of_interest$matching_F[e] <- diff[which.min(abs(diff$V2)),]$Fvals
    FLEPs_of_interest$diff[e] <- diff[which.min(abs(diff$V2)),]$V2
  }
  
  return(FLEPs_of_interest)
  
}

# --------------------
# (11) calculate F values for exploitation index - using Leslie matrix
calculate_F_from_EI_Leslie <- function(A,vul_a,Fvals,M,EI){
  
  #testing
  # out <- assemble_Leslie(maxage=as.numeric(parms$maxage[18]),
  #                        L_a=La_list[[18]],
  #                        B_a=Ba_list[[18]],
  #                        propmat_a=prop_list[[18]],
  #                        eggs_a=eggs_list[[18]],
  #                        vul_a=vul_list[[18]],
  #                        M=as.numeric(parms$M[18]),
  #                        Fval=0) #F=0, change this if want fishing
  # lep <- sum(out$LTABLE$LEP_a)
  # Asim <- out$A
  # Asim[1,] <- (Asim[1,]/(lep*adjFec))
  # A=Asim
  # vul_a=vul_list[[18]]
  # Fvals=c(0)
  # M=as.numeric(parms$M[18])
  # 
  # create vector of ages
  Age=1:length(A[1,]) 
  
  # assemble LTABLE: 
  L = data.frame(cbind(Age,vul_a))  
  
  # survival at age for different F values
  surv_at_age_by_F <- matrix(NA,nrow=length(Age),ncol=length(Fvals))
  survivorship_by_F <- matrix(0, nrow=length(Age),ncol=length(Fvals))
  LEP_at_age_by_F <- matrix(0,nrow=length(Age),ncol=length(Fvals))
  for(f in 1:length(Fvals)){ # step through F values
    
    for(a in 1:length(Age)) { #step through ages
      Ftemp = L$vul_a[a]*Fvals[f]	# FISH: F rate = vulnerability-at-age * F  
      surv_at_age_by_F[a,f] = exp(-(Ftemp+M)) #SURV:fraction surviving at each age
    } #closes loop calculating survival at age
    
    # set up matrix for survivorship (amount or fraction present at age)
    survivorship_by_F[1,] <- 1
    
    # for each age calc survival from age 1 to age a
    for(k in 1:(nrow(survivorship_by_F)-1)){ # step through ages  
      survivorship_by_F[k+1,f] = survivorship_by_F[k,f]*surv_at_age_by_F[k,f]}
    
    # set up matrix for LEP at age
    # LEP at age = fecundity-in-top-row * surv from age 0 to age a
    # LEP_at_age_by_F is spawning biomass distribution
    LEP_at_age_by_F[,f] <-A[1,]*survivorship_by_F[,f]
    #check new LEP (only when F=0 should new LEP=1.1):
    #sum(LEP_at_age_by_F[,f])
    
  }#closes f loop
  
  # Get FLEP from LEP
  LEP <- colSums(LEP_at_age_by_F)
  FLEP <- LEP/LEP[1]
  FLEP_by_F <- data.table(cbind(Fvals,FLEP))
  
  #What are corresponding FLEPs for EI values?
  FLEPs_of_interest <- data.frame((1-(EI*0.9)))
  colnames(FLEPs_of_interest) <- "FLEP_EI"
  FLEPs_of_interest$matching_F <- rep(NA,length=length(FLEPs_of_interest))

  for(e in 1:length(FLEPs_of_interest[,1])){
    diff <- data.frame(cbind(FLEP_by_F, (FLEP_by_F$FLEP-FLEPs_of_interest$FLEP_EI[e])))
    FLEPs_of_interest$matching_F[e] <- diff[which.min(abs(diff$V2)),]$Fvals
    FLEPs_of_interest$diff[e] <- diff[which.min(abs(diff$V2)),]$V2
  }
  rm(LEP,FLEP,L,surv_at_age_by_F,survivorship_by_F,Ftemp,LEP_at_age_by_F,f,a,k,e,diff)
  return(FLEPs_of_interest)
  #return(FLEP_by_F)
}


# --------------------
# (12) plot wavelet
# this is the actual sample code to export plot (in Patrick's code)
# tiff(file.path(".", "output_ms", "Fig3revised_100.tiff"),units="in", res=300, width = 6.85, height = 6.85)
# #postscript(file.path(".", "output_ms", "Fig3revised_100.ps"), width = 6.85, height = 6.85)
# # pdf(file.path(".", "output_ms", "Fig_3_summaryFreqContTS_Noise.pdf"), width = 6.85, height = 6.85)
# plot_gen_freq_wvlt(noise = noiseList,
#                    burn_in_pd = (burn_in + phasein_len),
#                    num_rows2plt = 200,
#                    n = plot_idx, 
#                    J1 = trunc((log(32/(2 * 1))/log(2))/0.01))
# dev.off()

plot_gen_freq_wvlt <- function(noise = df.list, 
                               burn_in_pd = rm_first_timesteps, #remove from beginning of ts
                               num_rows2plt = 1000, #number of years to plot
                               n = 1, #for indexing column n in noise df?
                               J1 = trunc((log(32/(2 * 1))/log(2))/0.01),
                               timesteps = timesteps,
                               span.multiplier = 1,
                               spp =  ) {
  # ---
  test = subset(parms[parms$problem_spp=="no",]$spp,select=c("spp","maxage"))
  test$maxage <- as.numeric(test$maxage)
  test <- test[order(test$maxage),]
  spp <- test$spp
  df.list <- df.list[spp] #set order
  line_d <- 2 #which margin to start counting at zero, 2=left
  ylim <- c(0,2)
  tsylim <- c(3000,4500)
  spwn_ylim <- c(0,0.35)
  spwn_xlim <- c(0,110)
  rows2plot <- (burn_in_pd+1):(burn_in_pd+num_rows2plt)
  n_rows <- seq(from=1,to=length(rows2plot))
  old <- par(mar = c(3,2,2,1), cex = .7)
  # set Daniell smoother for frequency response plot
  tmp <- ceiling(sqrt(length(1:(timesteps-burn_in_pd-1)))) #sq root of timeseries lgth, rounded
  if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  m = m * span.multiplier
  
  # mtext()
  #side 2 is left
  #las 1:horizontal 2:perpendicular to axis 3:vertical
  #at position along y axis
  #line on which MARgin line, starting at 0 counting outwards
  # ---
  
  tiff(file='C:/Users/provo/Documents/GitHub/pfmc/results/wvlet_panels_1-7_frqrsp_rotate_55_v2.tiff', units="in", width=8, height=11, res=300) 
  
  par(mfrow=c(7,4),mai=c(0.2,0.3,0.2,0.2))
  
  # --------------
  # 1) White noise
  # Label 
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "white noise", cex = 1.3)
  
  # Plot recruit time series
  plot(x = n_rows,
       y = whitenoise[rows2plot], cex=0.8, ylab="recruit",
       type = "l", xaxs = "i", yaxs = "i")
  #mtext("a", side = 2, las = 1, at = max(whitenoise), line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, whitenoise[rows2plot]), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("b", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  # Plot frequency response
  #m = 55
  p_spec <- spec.pgram(x=whitenoise,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="")
  #mtext("c", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  
  # --------------
  # 2) Sardine
  i=1
  
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
          x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
          y0=0,y1=0.33,lty=2,col="red",lwd=2)
  
  # Plot recruit time series
  #par(fig=c(0.1, 0.4, 0.80, 1), new = TRUE)
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("d", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  #par(fig=c(0.4, 0.7, 0.80, 1), new = TRUE)
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("e", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  # Plot frequency response
  #par(fig=c(0.7, 1, 0.8, 1), new = TRUE)
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("f", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  
  # --------------
  # 3) Bluefin
  i=2
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
       #mtext("g", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("h", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("i", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 4) Cabezon
  i=3
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("j", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("k", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("l", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 5) Canary
  i=4
  # Plot Label
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("m", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("n", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("o", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 6) Chilipepper
  i=5
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("p", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("q", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("r", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 7) Darkblotched
  i=6
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("s", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("t", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
 
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("u", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  dev.off()
  
  # --------
  tiff(file='C:/Users/provo/Documents/GitHub/pfmc/results/wvlet_panels_8-14_frqrsp_rotate_55_v2.tiff', units="in", width=8, height=11, res=300) 
  
  par(mfrow=c(7,4),mai=c(0.2,0.3,0.2,0.2))
  # 8) Dover sole
  i=7
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,
       label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                   "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("v", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==parms$spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("w", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("x", side = 2, las = 1, at = 0.5, line = line_d, cex = 1.2)
  
  # 9) English
  i=8
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,
       label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                   "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("y", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==parms$spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("z", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("aa", side = 2, las = 1, at = 0.5, line = line_d, cex = 1.2)
  
  # 10) Lingcod
  i=9
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,
       label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                   "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("bb", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==parms$spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("cc", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("dd", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 11) Hake
  i=10
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,
       label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                   "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("ee", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==parms$spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("ff", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("gg", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 12) Petrale
  i=11
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,
       label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                   "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("hh", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==parms$spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("ii", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("jj", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 13) Sablefish
  i=12
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,
       label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                   "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("kk", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==parms$spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("ll", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("mm", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  # 14) Widow
  i=13
  # Label 
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxs = "i", cex=0.8, 
       xlab="Age",ylab="p(spawning)",xlim=spwn_xlim)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,sep=""),
                               cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = df.list[[i]][df.list[[i]]$EI == 0 & df.list[[i]]$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxs = "i", cex=0.8, ylab="recruit")
  #mtext("nn", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  # plot wavelet power spectrum
  white.wt <- wt(cbind(1:num_rows2plt, 
                       ts.data[ts.data$spp==parms$spp[i] & 
                                 ts.data$EI==0 & 
                                 ts.data$year %in% rows2plot,]$value), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  #mtext("oo", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # Plot frequency response
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  p_spec <- spec.pgram(x=df.list[[i]][df.list[[i]]$EI == 0,]$value,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="",ylim=c(2,270000),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("pp", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  dev.off()
  
  
  # # plot wavelet power spectrum
  # par(fig=c(0.55, 0.9, 0, 0.2), new = TRUE)
  # p_one_over_f.wt <- wt(cbind(1:num_rows2plt, noiseList[[5]][rows2plot,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  # plot(p_one_over_f.wt, xlab = "Year")
  # mtext("j", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  # 
  # par(fig=c(0.4, 1, 0, 1), new = TRUE)
  # image.plot(legend.only = TRUE, zlim = c(0,64), add = FALSE) 
  # par(old)
  # 
} # end plot_gen_freq_wvlt
