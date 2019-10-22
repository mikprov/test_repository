# Run the simulation model
# by: mikaela provost
# last edited: Feb 19, 2019
# ===================================================================

library(tidyr)
library(dplyr)
library(gridExtra)
# ---
# load functions
source("C:/Users/Mikaela/Documents/GitHub/pfmc/src/functions.r")
# load parms df
parms = read.csv("C:/Users/Mikaela/Documents/GitHub/pfmc/parms_data/pfmc_parms.csv",
                 header=TRUE,stringsAsFactors = FALSE)
parms = as.data.frame(parms)
# notes on parms: all cols may be 'chr' because of words used in columns for numbers. 

# *************************************** #
# Plan
# (1) convert age to length (vonB & Schnute) 
# (2) convert length(cm) to weight (kg) at age 
# (3) calculate eggs at age 
# (4) calculate proportion mature at age
# (5) calculate vulnerability-at-age to fishing 
# (6) plot eggs vs age (or length if necessary) and 
#     compare with stock assessment. Make sure matches
# *************************************** #
# This script calculates life table infomation needed
# to assemble Leslie matrices for different F values:
# L_a, B_a, propmat_a, eggs_a, vul_a
# maxage & M are stored in parms df
# assemble_Leslie(maxage,L_a,B_a,propmat_a,eggs_a,vul_a,M)


# *************************************** #
# (1) convert age to length (vonB & Schnute) 
La_list <- as.list(rep(NA,length=length(parms$spp)))
names(La_list) <- parms$spp #list of cm-at-age vectors

for(i in 1:length(parms$spp)){ #for each spp
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes"){ #if this is a problem species
    #fill in zeros into L_a vector
    La_list[[i]] <- 0
    
  } else { #if not a problem species, use growth equation to convert age to length
    
  if(parms[parms$spp == parms$spp[i],]$growthFUN == "vonb"){ #if vonb equation
    #define the vonb parameters maxage,t0,Linf,Kvonb
    maxage = as.numeric(parms[parms$spp==parms$spp[i],]$maxage)
    t0 = as.numeric(parms[parms$spp==parms$spp[i],]$tknot)
    Linf = as.numeric(parms[parms$spp==parms$spp[i],]$Linf)
    Kvonb = as.numeric(parms[parms$spp==parms$spp[i],]$Kvonb)
    #run the vonb model to get cm per age
    ages <- seq(from=1,to=maxage,by=1)
    L_a <- vonbertgrowth(maxage=maxage,t0=t0,Linf=Linf,Kvonb=Kvonb,ages=ages)
    rm(maxage,t0,Linf,Kvonb,ages)
    
  } else {
    #define schnute parameters maxage,t1,t2,L1,L2,Ksch
    maxage = as.numeric(parms[parms$spp==parms$spp[i],]$maxage)
    t1 = as.numeric(parms[parms$spp==parms$spp[i],]$t1)
    t2 = as.numeric(parms[parms$spp==parms$spp[i],]$t2)
    L1 = as.numeric(parms[parms$spp==parms$spp[i],]$L1)
    L2 = as.numeric(parms[parms$spp==parms$spp[i],]$L2)
    Ksch = as.numeric(parms[parms$spp==parms$spp[i],]$Ksch)
    #run the schnute model to get cm per age
    ages <- seq(from=1,to=maxage,by=1)
    L_a <- schnutegrowth(maxage=maxage,t1=t1,t2=t2,L1=L1,L2=L2,Ksch=Ksch,ages=ages)
    rm(t1,t2,L1,L2,Ksch,ages,maxage)}
    
    # store L_a vector
    La_list[[i]] <- L_a
  }
  print(i)
}
rm(i,L_a)


# *************************************** #
# (2) convert length(cm) to weight (kg) at age 
Ba_list <- as.list(rep(NA,length=length(parms$spp)))
names(Ba_list) <- parms$spp #list of kg-at-age vectors 

for(i in 1:length(parms$spp)){
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes"){ #if species is a problem spp, then
    Ba_list[[i]] <- 0 #store zero in list
    } else { #if not a problem species
  Ba_list[[i]] <- calc_wt_at_age(
    L_a=La_list[[i]],
    cmkga=as.numeric(parms[parms$spp == parms$spp[i],]$cmkga),
    cmkgb=as.numeric(parms[parms$spp == parms$spp[i],]$cmkgb)) }
  }
rm(i)


# *************************************** #
# (3) calculate eggs at age 
eggs_list <- as.list(rep(NA,length=length(parms$spp)))
names(eggs_list) <- parms$spp
for(i in 1:length(parms$spp)){
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes"){ #if species is a problem spp, then
    #store age vector of zeros in list
    eggs_list[[i]] <- 0
    
  } else { #if not a problem species, then use eggs functions to get eggs-at-age vector
    
  if(parms[parms$spp == parms$spp[i],]$eggsFUN == "fun1") {
    #define parms for fecFUN1 B_a,eggs_intercept,eggs_slope
    B_a = Ba_list[[i]]
    eggs_intercept = as.numeric(parms[parms$spp==parms$spp[i],]$eggsintercept)
    eggs_slope = as.numeric(parms[parms$spp==parms$spp[i],]$eggsslope)
    eggs <- eggsFUN1(B_a=B_a,intercept=eggs_intercept,slope=eggs_slope)
  } else {
    #define parms for fecFUN2
    B_a = Ba_list[[i]] 
    eggs_intercept = as.numeric(parms[parms$spp==parms$spp[i],]$eggsintercept)
    eggs_slope = as.numeric(parms[parms$spp==parms$spp[i],]$eggsslope)
    eggs <- eggsFUN2(B_a=B_a,intercept=eggs_intercept,slope=eggs_slope)
  }
  eggs_list[[i]] <- eggs #store vector of eggs in list
  rm(B_a,eggs_intercept,eggs_slope,eggs)
  }
}
rm(i)


# *************************************** #
# (4) calculate proportion mature at age
prop_list <- as.list(rep(NA,length=length(parms$spp)))
names(prop_list) <- parms$spp

for(i in 1:length(parms$spp)){
  
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes"){ #if species is a problem spp, then
    #store age vector of zeros in list
    prop_list[[i]] <- 0 
    
  } else { #if not a problem species, then prop mature function to get propmat-at-age vector
    
    if(parms[parms$spp == parms$spp[i],]$propmatFUN == "age"){ #if prop mat is in terms of age, then
      #use the calc_mat_at_age function as is, based on vector of ages
      maxage=as.numeric(parms[parms$spp==parms$spp[i],]$maxage)
      slope=as.numeric(parms[parms$spp==parms$spp[i],]$propmat_slope)
      inflection=as.numeric(parms[parms$spp==parms$spp[i],]$propmat_inflection)
      #calculate proportion mature at age, store in list
      prop_list[[i]] <- calc_mat_at_age(maxage=maxage,slope=slope,inflection=inflection)
    } else {
    
    if(parms[parms$spp == parms$spp[i],]$propmatFUN == "length"){ #if prop mat is in terms of length, then
     #used calc_mat_at_length, feed vector of lengths at age. Define inputs for function:
      lengths_at_age = La_list[[i]]
      slope=as.numeric(parms[parms$spp==parms$spp[i],]$propmat_slope)
      inflection=as.numeric(parms[parms$spp==parms$spp[i],]$propmat_inflection)
      prop_list[[i]] <- calc_mat_at_length(lengths_at_age=lengths_at_age,slope=slope,inflection=inflection)}
    } 
  
  
  }
}
rm(i,maxage,slope,inflection,lengths_at_age)
prop_list$Bluefin <- c(0,0,0.2,0.5,1,1,rep(1,length=20))


# *************************************** #
# (5) calculate vulnerability-at-age to fishing 
# for now, using proportion mature
vul_list <- prop_list


