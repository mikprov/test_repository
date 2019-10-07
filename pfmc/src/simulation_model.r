# Simulation model
# by: Mikaela Provost



# ===================================================================
# 1) define sim model
sim_model <- function(A,timesteps,alpha,beta,sig_r,initial_eggs,noise) {
  # testing
  # A = A3d_list[[2]][,,1]
  # timesteps = 5000
  # alpha = alphas 
  # beta = betas
  # sig_r = 0.3
  # initial_eggs = betas
  # noise = whitenoise

  age_at_mat <- 1 # set inequality to min threshold for maturity
  maxage = length(A[,1])
  ages = length(seq(1,maxage,by=1))
  N0 = c(initial_eggs, rep(0,ages-1)) #vector, length=num of ages, first age is initial eggs
  
  Nt = matrix(0,ages,timesteps) #empty matrix to store output: rows=ages, cols=time
  #Initialize vector of population sizes with extra 
  #rows for egg production (top value in age vector) & 
  #recruitment before variability (output from BH)
  
  Nt[,1] = N0 #Put in initial values
  Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Leslie to get 2nd age vector
  eggs = c(initial_eggs, rep(NA,timesteps-1)) #will save egg production here, this will be the input in BH
  
  recruits0 = eggs[1]/( (1/alpha) + (eggs[1]/beta) )
  recruits = c(recruits0, rep(0, timesteps-1)) #will save recruits here (output from BH)
  
   for(t in 1:(timesteps-2)){ #step through time
     #Save egg production, this is new number of age 1 individuals
     eggs[t+1] = Nt[1,t+1]
     #save recruits, treat egg production as spawners in BH
     recruits[t+1] = eggs[t+1]/( (1/alpha) + (eggs[t+1]/beta) ) #((alpha*eggs[t+1])/(1+beta*eggs[t+1]))
     #replace age 1 with recruits from BH, add noise
     Nt[1,t+1] = (recruits[t+1])*exp(sig_r*noise[t+1])
     #perform population projection for one time step
     Nt[,t+2] = A %*% Nt[,t+1]
  
   }
   #Nsize = colSums(Nt) #total population size, sums num at age for each timestep
  Nsize = colSums(Nt[age_at_mat:maxage,]) #sum rows that correspond to spawning adults
  #N_t = Nt[1,][2:(timesteps-1)] #Nt is number of spawning adults, aka eggs, after environmental noise
  eggs = eggs[2:(timesteps-1)] #egg production
  recruits = recruits[2:(timesteps-1)] #recruits produced via BH, before influence of environmental noise
  Nsize = Nsize[-c(4999,5000)] #
  return(list(eggs=eggs, recruits=recruits, Nsize=Nsize))

  
  # Nt[,1] = N0 #Put in initial values
  # Nt[,2] = A %*% Nt[,1]  # multiply initial age vector with Jacobian to get 2nd age vector
  # eggs = c(initial_eggs, rep(NA,timesteps-2)) #will save egg production here, this will be the input in BH
  # eggs[2] = Nt[1,2]
  # recruits1 = eggs[1]/( (1/alpha) + (eggs[1]/beta) )
  # recruits = c(0,rep(NA, timesteps-2)) #will save recruits here (output from BH)
  # recruits[2] = recruits1 #we are ignoring recorinding recruits at time 0
  # 
  # 
  # for(t in 1:(timesteps-2)){ #step through time
  #   #Save egg production, this is new number of age 1 individuals
  #   eggs[t+1] = Nt[1,t+1] 
  #   #save recruits, treat egg production as spawning biomass in BH
  #   recruits[t+1] = eggs[t]/( (1/alpha) + (eggs[t]/beta) ) #((alpha*eggs[t+1])/(1+beta*eggs[t+1])) 
  #   #replace age 1 in numbers-at-age matrix with recruits from BH, add noise    
  #   Nt[1,t+1] = (recruits[t+1])*exp(sig_r*noise[t]) 
  #   #perform population projection for one time step   
  #   Nt[,t+2] = A %*% Nt[,t+1] 
  # }
  # Nsize = colSums(Nt[age_at_mat:maxage,]) #sum rows that correspond to spawning adults
  # eggs = eggs[1:(timesteps-1)] #egg production
  # recruits = recruits[1:(timesteps-1)] #recruits produced via BH, before influence of environmental noise
  # return(list(eggs=eggs, recruits=recruits, Nsize=Nsize))
}

# ===================================================================


