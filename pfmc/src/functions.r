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
calculate_F_from_EI_Leslie <- function(A,vul_a,Fval,M,EI){
  
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
  surv_at_age_by_F <- matrix(NA,nrow=length(Age),ncol=length(Fval))
  survivorship_by_F <- matrix(0, nrow=length(Age),ncol=length(Fval))
  LEP_at_age_by_F <- matrix(0,nrow=length(Age),ncol=length(Fval))
  for(f in 1:length(Fval)){ # step through F values
    
    for(a in 1:length(Age)) { #step through ages
      Ftemp = L$vul_a[a]*Fval[f]	# FISH: F rate = vulnerability-at-age * F  
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
  FLEP_by_F <- data.table(cbind(Fval,FLEP))
  
  #What are corresponding FLEPs for EI values?
  FLEPs_of_interest <- data.frame((1-(EI*0.9)))
  colnames(FLEPs_of_interest) <- "FLEP_EI"
  FLEPs_of_interest$matching_F <- rep(NA,length=length(FLEPs_of_interest))

  for(e in 1:length(FLEPs_of_interest[,1])){
    diff <- data.frame(cbind(FLEP_by_F, (FLEP_by_F$FLEP-FLEPs_of_interest$FLEP_EI[e])))
    FLEPs_of_interest$matching_F[e] <- diff[which.min(abs(diff$V2)),]$Fval
    FLEPs_of_interest$diff[e] <- diff[which.min(abs(diff$V2)),]$V2
  }
  rm(LEP,FLEP,L,surv_at_age_by_F,survivorship_by_F,Ftemp,LEP_at_age_by_F,f,a,k,e,diff)
  return(FLEPs_of_interest)
  #return(FLEP_by_F)
}




# --------------------
# (13) Simulate population
# Simulate a population given the: noise, harvest level (which FLEP?), 
# threshold FLEP & calc mean for threshold level. Return: N, recruits,
# eggs time series in large data frame for plotting.
# Note: subset spp, List of Leslie3d to non-problem spp


simulate_spp <- function(noise,
                         timesteps,
                         rm_first_timesteps, #burn in time
                         alphas, #for BH
                         betas, #for BH
                         span.multiplier, #not sure I need this yet
                         EIlevels, #vector of EI levels, 1st is OFL--match
                         OFL, #FLEP associated with OFL
                         List_of_Leslie3d, #should have F vals associated with FLEPlevels
                         spp,
                         sig_r) {
  #Testing
  # noise = white
  # timesteps = 5000
  # rm_first_timesteps = 2000 #burn in time
  # alphas = 1.2 #for BH
  # betas = 5000 #for BH
  # span.multiplier = 1 #not sure I need this yet
  # EIlevels = EIlevels #vector of EI levels, 1st is OFL--match
  # OFL = EIlevels[3] #FLEP associated with OFL
  # List_of_Leslie3d = A3d_list #should have F vals associated with FLEPlevels
  # spp = parms[parms$problem_spp == "no",]$spp[1:3]
  # sig_r=0.3
  
  # Loop over Aarray list to simulate using different Leslie matrices
  # set params for simulation:
  # Note: A3d_list --> each Leslie matrix in array is associated with F
  # value that matches EI value. 
  output.3d.list <- as.list(rep(NA,length=length(List_of_Leslie3d))) #store timeseries here
  names(output.3d.list) <- spp
  for (i in 1:length(List_of_Leslie3d)) { #step through each pop
    Leslie3d = List_of_Leslie3d[[i]] #select the 3d array of Leslie matricies
    # array dims: row=ts length, col=3 is number of ts (eggs,recruits,Nsize), depth=F vals
    output.matrix <- array(NA,c(timesteps-2,3,length(Leslie3d[1,1,]))) 
     for (f in 1:length(Leslie3d[1,1,])) { #step through each Leslie matrix (for each F value)
      output = sim_model(A=Leslie3d[,,f], timesteps=timesteps, 
                           alpha=alphas, beta=betas, 
                           initial_eggs=betas,sig_r,
                           noise=noise) #specify noise
      output.matrix[,,f] <- do.call(cbind,output) #fill in array for pop i
      }#close f loop
   output.3d.list[[i]] <- output.matrix #store simulation output in list
   print(i)
  }#closes i loop for spp
  
  # *************************************** #
  # Format output ts for plotting simulations using output.3d.list
  # At this point I have one important object:
  # 1. [output.3d.list] a list of 3d arrays. Each array is timeseries output
  #    from simulations at different F levels. 
  variable_type <- c("eggs","recruits","Nsize")
  # --- reorganize recruit timeseries data --- #
  var.number <- which(variable_type == "recruits") # recruits
  df.list <- as.list(rep(NA,length=length(spp)))
  names(df.list) <- spp
  for (i in 1:length(output.3d.list)) {
    aa <- as.data.frame(output.3d.list[[i]][,var.number,])
    aa$year <- seq(from=1, to=length(aa[,1]),by=1)
    colnames(aa) <- c(EIlevels,"year")
    aa1 <- aa %>% gather(EI,value,1:length(EIlevels))
    aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
    aa1$spp <- rep(spp[i],length=length(aa[,1]))
    aa1$EI <- as.numeric(as.character(aa1$EI))
    aa1$means <- rep(0,length=length(aa1[,1]))
    means <- rep(0,length=length(unique(aa1$EI)))
    for(f in 1:length(unique(aa1$EI))){
      means[f] <- mean(aa1[aa1$EI==unique(aa1$EI)[f] &
                             aa1$year %in% rm_first_timesteps:(timesteps-2),]$value)
      aa1[aa1$EI==unique(aa1$EI)[f],]$means <- rep(means[f],
                                                   length=length(aa1[aa1$EI==unique(aa1$EI)[f],1]))}
    df.list[[i]] <- aa1
    }
  recruits.ts <- bind_rows(df.list,id=NULL)
  recruits.ts$peak <- spawndistmetrics[match(recruits.ts$spp,spawndistmetrics$spp),"mode_age"]
  recruits.ts$cvs <- spawndistmetrics[match(recruits.ts$spp,spawndistmetrics$spp),"cvs_mode"]
  recruits.ts$maxage <- spawndistmetrics[match(recruits.ts$spp,spawndistmetrics$spp),"maxage"]
  recruits.ts$stdev <- spawndistmetrics[match(recruits.ts$spp,spawndistmetrics$spp),"sd_mode"]
  rm(i,f,means,aa1,aa,var.number,df.list) #clean up
  
  # --- reorganize egg timeseries data --- #
  var.number <- which(variable_type == "eggs") # recruits
  df.list <- as.list(rep(NA,length=length(spp)))
  names(df.list) <- spp
  for (i in 1:length(output.3d.list)) {
    aa <- as.data.frame(output.3d.list[[i]][,var.number,])
    aa$year <- seq(from=1, to=length(aa[,1]),by=1)
    colnames(aa) <- c(EIlevels,"year")
    aa1 <- aa %>% gather(EI,value,1:length(EIlevels))
    aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
    aa1$spp <- rep(spp[i],length=length(aa[,1]))
    aa1$EI <- as.numeric(as.character(aa1$EI))
    aa1$means <- rep(0,length=length(aa1[,1]))
    means <- rep(0,length=length(unique(aa1$EI)))
    for(f in 1:length(unique(aa1$EI))){
      means[f] <- mean(aa1[aa1$EI==unique(aa1$EI)[f] &
                             aa1$year %in% rm_first_timesteps:(timesteps-2),]$value)
      aa1[aa1$EI==unique(aa1$EI)[f],]$means <- rep(means[f],
                                                   length=length(aa1[aa1$EI==unique(aa1$EI)[f],1]))}
    df.list[[i]] <- aa1
    }
  egg.ts <- bind_rows(df.list,id=NULL)
  egg.ts$peak <- spawndistmetrics[match(egg.ts$spp,spawndistmetrics$spp),"mode_age"]
  egg.ts$cvs <- spawndistmetrics[match(egg.ts$spp,spawndistmetrics$spp),"cvs_mode"]
  egg.ts$maxage <- spawndistmetrics[match(egg.ts$spp,spawndistmetrics$spp),"maxage"]
  egg.ts$stdev <- spawndistmetrics[match(egg.ts$spp,spawndistmetrics$spp),"sd_mode"]
  rm(i,f,means,aa1,aa,var.number,df.list) #clean up
  
  # --- reorganize N timeseries data --- #
  var.number <- which(variable_type == "Nsize") # recruits
  df.list <- as.list(rep(NA,length=length(spp)))
  names(df.list) <- spp
  for (i in 1:length(output.3d.list)) {
    aa <- as.data.frame(output.3d.list[[i]][,var.number,])
    aa$year <- seq(from=1, to=length(aa[,1]),by=1)
    colnames(aa) <- c(EIlevels,"year")
    aa1 <- aa %>% gather(EI,value,1:length(EIlevels))
    aa1$variable <- rep(variable_type[var.number],length=length(aa1[,1]))
    aa1$spp <- rep(spp[i],length=length(aa[,1]))
    aa1$EI <- as.numeric(as.character(aa1$EI))
    aa1$means <- rep(0,length=length(aa1[,1]))
    means <- rep(0,length=length(unique(aa1$EI)))
    for(f in 1:length(unique(aa1$EI))){
      means[f] <- mean(aa1[aa1$EI==unique(aa1$EI)[f] &
                             aa1$year %in% rm_first_timesteps:(timesteps-2),]$value)
      aa1[aa1$EI==unique(aa1$EI)[f],]$means <- rep(means[f],
                                                   length=length(aa1[aa1$EI==unique(aa1$EI)[f],1]))}
    df.list[[i]] <- aa1
  }
  N.ts <- bind_rows(df.list,id=NULL)
  N.ts$peak <- spawndistmetrics[match(N.ts$spp,spawndistmetrics$spp),"mode_age"]
  N.ts$cvs <- spawndistmetrics[match(N.ts$spp,spawndistmetrics$spp),"cvs_mode"]
  N.ts$maxage <- spawndistmetrics[match(N.ts$spp,spawndistmetrics$spp),"maxage"]
  N.ts$stdev <- spawndistmetrics[match(N.ts$spp,spawndistmetrics$spp),"sd_mode"]
  
  #combine all 3 ts dfs
  ts.data <- bind_rows(recruits.ts,egg.ts,N.ts)
  return(ts.data)
  rm(i,f,means,aa1,aa,var.number,df.list,recruits.ts,egg.ts,N.ts) #clean up
  
}


# (14) Generate list of 3d Leslie arrays for different f levels
# ***

generate3dLeslie <- function(spp,Fvals,EIlevels){
  # create empty matrices
  A3d_list = as.list(rep(NA,length(spp))) #Leslie arrays in list
  LTABLE_list <- as.list(rep(NA,length(spp))) #LTABLE arrays in list
  names(A3d_list) <- spp
  names(LTABLE_list) <- spp
  eigan1 <- matrix(NA,nrow=length(Fvals),ncol=length(spp))
  eigan12 <- matrix(NA,nrow=length(Fvals),ncol=length(spp))
  LEPs <- matrix(NA,nrow=length(Fvals),ncol=length(spp))
  constLEP <- matrix(NA,nrow=length(Fvals),ncol=length(spp))
  EI_F_LEP1.1_L <- as.list(rep(NA,length=length(spp)))
  names(EI_F_LEP1.1_L) <- spp
  
  system.time(for (i in 1:length(spp)){ #for each spp
    # define parms for spp i
    M = as.numeric(parms[parms$spp == spp[i],]$M)
    maxage = as.numeric(parms[parms$spp == spp[i],]$maxage)
    
    # create empty arrays for Leslie and LTABLE, vectors for eigens 
    A3d <- array(NA,c(maxage,maxage,length(Fvals)))
    LTABLE3d <- array(NA,c(maxage,15,length(Fvals))) 
    e1 <- rep(NA,length=length(Fvals))
    e2 <- rep(NA,length=length(Fvals))
    e12 <- rep(NA,length=length(Fvals))
    lep <- rep(NA,length=length(Fvals))
    
    # assemble Leslie for each value of Fvals
    for(f in 1:length(Fvals)){
      out <- assemble_Leslie(maxage=maxage,
                             L_a=La_list[spp][[i]],
                             B_a=Ba_list[spp][[i]],
                             propmat_a=prop_list[spp][[i]],
                             eggs_a=eggs_list[spp][[i]],
                             vul_a=vul_list[spp][[i]],
                             M=M,
                             Fval=Fvals[f]) #F=0, change this if want fishing
      LTABLE3d[,,f] <- data.matrix(out$LTABLE) #store LTABLE in 3d array
      # calculate LEP at F=0 (use this value of LEP to standardize unfished fecundity)
      lep_unfished = sum(out$LTABLE$LEP_a)
      # Asim: adjust fecundities so LEP equal for all species - simulation matrix
      Asim <- out$A
      Asim[1,] <- (Asim[1,]/(lep_unfished*adjFec))
      A3d[,,f] <- Asim
      # Aeig: adjust fecundities and multiply by k - eigen analysis matrix
      Aeig <- out$A
      Aeig[1,] <- (Aeig[1,]/(lep_unfished*adjFec)) # *Fvals[f]
      
      # re-calculate adjusted LEP w/o F (should equal conLEP)
      term <- rep(NA,length=maxage)
      for(a in 1:maxage){term[a] <- Asim[1,a]*out$LTABLE$Survship[a]} #fecundity*survivorship
      constLEP[f,i] <- sum(term) #plotting 'term' is the spawning distribution over age
      rm(term,a)
      
      # calculate adjusted LEP w/F
      term <- rep(NA,length=maxage)
      for(a in 1:maxage){term[a] <- Asim[1,a]*out$LTABLE$Survship_Fing[a]} #fec*survshp w/F
      LEPs[f,i] <- sum(term)
      rm(term,a)
      
      # extract lambda1, lambda2/lambda1
      e1[f] = extract_first_eigen_value(Aeig)
      e2[f] = extract_second_eigen_value(Aeig)
      e12[f] = e2[f] / e1[f]    } #closes F loop
   
      EI_F_LEP1.1_L[[i]] <- calculate_F_from_EI_Leslie(A=A3d[,,1],
                                                       vul_a=vul_list[spp][[i]],
                                                       Fval=Fvals,
                                                       M=M,
                                                       EI=EIlevels)
      saveFs <- which(Fvals %in% EI_F_LEP1.1_L[[i]]$matching_F) #where are Fs of interest?
      A3dsave <- A3d[,,saveFs] #save only Leslie matrices that match Fs
      
      A3d_list[[i]] <- A3dsave 
      LTABLE_list[[i]] <- LTABLE3d
      eigan1[,i] <- e1
      eigan12[,i] <- e12
      LEPs[,i] <- lep
      #clean up before next i in loop
      rm(A3d,A3dsave,LTABLE3d,e1,e12,lep,out,M,maxage,saveFs) 
      print(spp[i]) #which spp did I finish?
    
  }) #closes spp loop
  return(A3d_list)
}#closes function
