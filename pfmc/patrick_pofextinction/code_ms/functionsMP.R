# Population vulnerability modelling ####

# Determine the delta_e (fraction of spawning females spawning at age-3)
# returns a plot of N3/(N3+N4) vs delta_e for a given age-3 survival
# Also returns the de for the the wanted fraction (N3/(N3+N4)) and age-3 survival
calc_de <- function(surv3, wanted_frac, plot_wanted = "no") {
  # For Butte Creek chinook salmon this function calculates the delta_e
  # value that returns the "wanted fraction", which is the fraction of 
  # spawners that return at age three assuming equilibrium conditions
  # for a given/assumed value of ocean survival from age-3 to age-4 (surv3)
  # The number of age-3 spawners at EQ is: surv1*surv2*delta_e
  # The number of age-3 spawners at EQ is: surv1*surv2*(1-delta_e)*surv3
  # 
  # Limited data from coded wire tag returns of spawning fish suggest that 
  # BC spring-run chinook spawn at ages 3 and 4, with sometimes as few as 20% 
  # spawning at age-3 in a given year to as much as 80% spawning at age-3.
  #
  #
  # It also returns a plot of the fraction of age-3 and age-4 to total spawnwers
  # [age-3 + age-4 spawners] vs. delta_e
  
  n <- 100
  r3_S <- r4_S <- numeric(n)
  de <- (1:n)/n
  
  for (i in 1:n) {
    r3_S[i] <- de[i]/(de[i] + surv3*(1-de[i]))
    r4_S[i] <- (surv3*(1-de[i]))/(de[i] + surv3*(1-de[i]))
  }
  
  if (plot_wanted == "yes") {
    plot(de, r3_S, type = "l", col = "slateblue", lwd = 3,
         xlab = "Fraction spawning at age-3 (delta_e)",
         ylab = "Number of Age-a Spawners/Total Number of Spawners")
    lines(de, r4_S, type = "l", col = "red", lwd = 3)
    abline(h=0.5, lty = 2)
    abline(v=0.5, lty = 2)
    title(paste("Age-3 ocean survival = ", surv3, sep = ""))
    legend("right", c("a=3","a=4"), # puts text in the legend 
           lty=c(1,1), # gives the legend appropriate symbols (lines)
           lwd=c(2.5,2.5),
           col=c("slateblue","red"),
           cex = .7,
           bty = "n") # gives the legend lines the correct color and width
  }
  
  desired_de <- wanted_frac*surv3/(wanted_frac*(surv3 - 1) + 1)
  #print("The delta_e that you want is:")
  #print(desired_de)
  return(desired_de)
} # end bracket of calc_de()

# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
SPR_srcs <- function(surv1, surv2, surv3, delta_e, survPS = 1) {
  SPR <- surv1*surv2*delta_e*survPS + surv1*surv2*surv3*(1-delta_e)*survPS
  return(SPR)
} # end bracket for SPR_scrs

# plot how SPR changes with early ocean survival and fraction of age-3 spawners
plotSPR_de_surv1 <- function (sprDF) {
  ggplot(sprDF, aes(x = delta_e, y = SPR, colour = factor(surv1))) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("black","grey20", "grey40", "grey70"), name = "Early Ocean \nSurvival") +
    ylab("Spawners Per Recruit") +
    xlab("Fraction of spawners returning at age-3") +
    xlim(c(0, 1)) +
    ylim(c(0, 0.08)) +
    theme_bw() +
    theme(panel.grid.major = theme_line(colour = NA), panel.grid.minor = theme_line(colour = NA)) 
}


# function to calculate beta of Beverton-Holt (where alpha/beta = max recruits)
# based on a given alpha, equilibrium spawning stock size (EQsp) and SPR level

calc_beta <-function(alpha, EQsp, spr) {
  # For a given SPR, alpha and EQsp, we can calculate beta, for Butte Creek 
  # spring run chinook salmon
  # To calculate beta for multiple alpha values, just pass a vecter of alphas
  # to the alpha argument
  #if (length(EQsp) == 1 & all(EQsp == "infinity"))  beta <- 0 else beta <- (alpha * spr - 1)/EQsp
  
  beta <- (alpha * spr - 1)/EQsp
  if (any(beta < 0)) warning("beta value(s) is (are) negative - must be positive")
  return(beta)
}

# function to plot SR relationships investigated in population model

plotStockRecruit <- function(alpha, beta, EQsp, delta_e, spr, make_plot) {
  
  spawners <- seq(0, max(EQsp)*2, len = 200)
  
  if (length(alpha) > 1 & length(EQsp) > 1) stop("Only alpha or beta (EQsp) can be varied in a function call")
  if (make_plot == TRUE) {
    if (length(alpha) > 5) warning("only 1 to 5 values of alpha may be plotted at a time")
    if (length(EQsp) > 5) warning("only 1 to 5 values of EQsp may be plotted at a time")
  }
  #if (is.empty(alpha) | is.empty(beta)) stop("need to enter alpha and beta")
  
  # Set plotting colors
  if (length(alpha)  == 2 | length(EQsp) == 2) cols <- c("blue", "red")
  if (length(alpha)  == 3 | length(EQsp) == 3) cols <- c("blue", "black", "red")
  if (length(alpha)  == 4 | length(EQsp) == 4) cols <- c("blue", "dodgerblue", "plum", "red")
  if (length(alpha)  == 5 | length(EQsp) == 5) cols <- c("blue", "dodgerblue", "black", "plum",  "red")
  if (length(alpha) == 1 & length(EQsp) == 1) cols <- "black" 
  #   if (length(alpha) > 5) cols <- rainbow(n)
  #   if (length(EQsp) > 5) cols <- rainbow(n)
  # makes a plot of underlying SR relationship for population projections
  # if only one alpha and one beta
  
  if (length(alpha) > 1 & length(EQsp) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(alpha))
    alpha_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(alpha)) {
        recruits[i,j] <- alpha[j]*spawners[i]/(1 + beta[j]*spawners[i])
        alpha_names[j] <- paste("alpha_",toString(signif(alpha[j], 2)), sep = "")
      }
    }
    
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(alpha)+1)] <- alpha_names
    srm <- melt(sr, id.vars = "spawners")
    names(srm)[2:3] <- c("Alpha", "recruits")
    spawn4spr <- (alpha*spr - 1)/beta
    
    p <- ggplot(data = srm, aes(spawners, recruits, colour = Alpha) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1/spr, colour = "black", linetype = 2) +
      geom_point(data = data.frame(x = spawn4spr, y = alpha*spawn4spr/(1+beta*spawn4spr)), aes(x = x, y = y), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(recruits)*1.1))
  }
  
  if (length(EQsp) > 1 & length(alpha) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(EQsp))
    beta_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(beta)) {
        recruits[i,j] <- alpha*spawners[i]/(1 + beta[j]*spawners[i])
        beta_names[j] <- paste("EQsp_",toString(EQsp[j]), sep = "")
      }
    }
    
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(beta)+1)] <- beta_names
    srm <- melt(sr, id.vars = "spawners")
    names(srm)[2:3] <- c("Beta", "recruits")
    spawn4spr <- (alpha*spr - 1)/beta
    
    p <- ggplot(data = srm, aes(spawners, recruits, colour = Beta) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1/spr, colour = "black", linetype = 2) +
      geom_point(data = data.frame(x = spawn4spr, y = alpha*spawn4spr/(1+beta*spawn4spr)), aes(x = x, y = y), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(recruits)*1.1))
    
  }
  
  if (length(beta) == 1 & length(alpha) == 1) {
    recruits <- alpha*spawners/(1 + beta*spawners)
    # make df for ggplot call
    sr <- data.frame(spawners, recruits)
    
    # Spawning stock for given SPR
    spawn4spr <- (alpha*spr - 1)/beta
    #spawn4spr == EQsp
    #all.equal(EQsp, spawn4spr) # nearly equal, there's a slight rounding error
    p <- ggplot(data = sr, aes(spawners, recruits) ) +
      geom_line(colour = "slateblue") +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1/spr, colour = "red", line_type = 2) +
      geom_point(data = data.frame(x = spawn4spr, y = alpha*spawn4spr/(1+beta*spawn4spr)), aes(x = x, y = y), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(sr$recruits)*1.05)) 
  } 
  if (make_plot == TRUE) {
    print(p)
  }
} # end of plotStockRecruit()

plotStockSpEsc <- function(alpha, beta, EQsp, delta_e, surv1, surv2, surv3, spr, make_plot) {
  
  spawners <- seq(0, max(EQsp)*2, len = 200)
  
  if (length(alpha) > 1 & length(EQsp) > 1) stop("Only alpha or EQsp can be varied in a function call")
  if (make_plot == TRUE) {
    if (length(alpha) > 5) stop("only 1 to 5 values of alpha may be plotted at a time")
    if (length(EQsp) > 5) stop("only 1 to 5 values of EQsp may be plotted at a time")
  }
  #if (is.empty(alpha) | is.empty(beta)) stop("need to enter alpha and beta")
  
  # Set plotting colors
  if (length(alpha)  == 2 | length(EQsp) == 2) cols <- c("blue", "red")
  if (length(alpha)  == 3 | length(EQsp) == 3) cols <- c("blue", "black", "red")
  if (length(alpha)  == 4 | length(EQsp) == 4) cols <- c("blue", "dodgerblue", "plum", "red")
  if (length(alpha)  == 5 | length(EQsp) == 5) cols <- c("blue", "dodgerblue", "black", "plum",  "red")
  if (length(alpha) == 1 & length(EQsp) == 1) cols <- "black"  
  
  # makes a plot of underlying SR relationship for population projections
  # if only one alpha and one beta
  
  if (length(alpha) > 1 & length(EQsp) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(alpha))
    alpha_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(alpha)) {
        recruits[i,j] <- alpha[j]*spawners[i]/(1 + beta[j]*spawners[i])
        alpha_names[j] <- paste("alpha_",toString(signif(alpha[j], 2)), sep = "")
      }
    }
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(alpha)+1)] <- alpha_names
    srm <- melt(sr, id.vars = "spawners")
    srm$spEsc <- with(srm, value*surv1*surv2*delta_e + value*surv1*surv2*surv3*(1-delta_e))
    names(srm)[2:3] <- c("Alpha", "recruits")
    
    p <- ggplot(data = srm, aes(spawners, spEsc, colour = Alpha) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits: Returning Spawners") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
      #geom_point(aes(x = EQsp, y = EQsp), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(srm$spEsc)*1.1))
  }
  
  if (length(EQsp) > 1 & length(alpha) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(EQsp))
    beta_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(beta)) {
        recruits[i,j] <- alpha*spawners[i]/(1 + beta[j]*spawners[i])
        beta_names[j] <- paste("EQsp_",toString(EQsp[j]), sep = "")
      }
    }
    
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(beta)+1)] <- beta_names
    srm <- melt(sr, id.vars = "spawners")
    srm$spEsc <- with(srm, value*surv1*surv2*delta_e + value*surv1*surv2*surv3*(1-delta_e))
    names(srm)[2:3] <- c("Beta", "recruits")
    
    p <- ggplot(data = srm, aes(spawners, spEsc, colour = Beta) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits: Returning Spawners") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
      #geom_point(aes(x = EQsp, y = EQsp), colour = "black", shape = 22, size = 4) +
      xlim(c(0, max(srm$spEsc)*1.1)) + 
      ylim(c(0, max(srm$spEsc)*1.1))
    
  }
  
  if (length(beta) == 1 & length(alpha) == 1) {
    recruits <- alpha*spawners/(1 + beta*spawners)
    # make df for ggplot call
    sr <- data.frame(spawners, recruits)
    sr$spEsc <- with(sr, recruits*surv1*surv2*delta_e + recruits*surv1*surv2*surv3*(1-delta_e))
    # Spawning stock for given SPR
    spawn4spr <- (alpha*spr - 1)/beta
    #spawn4spr == EQsp
    #all.equal(EQsp, spawn4spr) # nearly equal, there's a slight rounding error
    p <- ggplot(data = sr, aes(spawners, spEsc) ) + 
      geom_line(colour = "slateblue") +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits: Returning Spawners") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1, colour = "red", line_type = 2) +
      #geom_point(aes(x = EQsp, y = EQsp), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(sr$spEsc)*1.05)) 
  } 
  if (make_plot == TRUE) {
    print(p)
  }
} # end of plotStockSpEsc()


# age-structured model that holds early ocean survival constant [low, mid, high, very high]
# only time varying parameter is pre-spawning survival
# Vary: alpha or beta f
# Constant: EQsp, wanted_frac, water management scenario, climate scenario, and initial age-structure

# surv1 <- 0.05
# surv2 <- 0.8
# surv3 <- 0.8
# survPS <- "vary"
# wanted_frac <- 0.5 
# EQsp <- 7500; alpha_vec <- c(1.5, 3, 5, 10)
# # EQsp <- c(2500, 7500, 15000); alpha_vec <- 2
# # EQsp <- 7500; alpha_vec <- 2
# wtrMgmt <- "NoD" # "NoD" or "BAU" ""
# climate <- "A2" # "B1" or "A2" ""
# variableSR <- FALSE
# yrs <- 2010:2099

popProjConstMarine <- function (surv1, 
                                surv2 = 0.8, 
                                surv3 = 0.8,
                                survPS = 1,
                                wanted_frac,
                                alpha_vec,
                                EQsp,
                                wtrMgmt = "", # can be "" for test case of survPS = 1
                                climate = "", # can be "" for test case of survPS = 1
                                yrs = 2010:2099,
                                make_plot,
                                variableSR = FALSE) {
  sim_path <- file.path(".", "data", "simulation_results")
  # some checks on model parameterization
  if (length(surv1) != 1) stop("This function assumes constant early ocean survival, please enter a single value")
  
  if (is.numeric(survPS)) {
    print("survPS = 1 and should yield constant number of spawners over time (EQsp), except with beta is 0")
  } 
  
  if (survPS == "vary") {
    print(paste("Investigating", wtrMgmt, "and", climate, "prespawn mortality scenarios"))
  }
  
  if (any(alpha_vec <= 1) | any(alpha_vec > 10)) stop("1/spr multiplier must be greater than 1 or less than or equal to 10") 
  
  
  # start parameterizing model
  # determine the value for delta_e based on the fraction of age-3 spawner to total spawners of interest
  # The delta_e calcs are based on longterm/equilibrium assumptions.
  # CADFG data are based on recovery BC SRCS CWT recoveries. Might not be definitive, but provide 
  # a range to test
  delta_e <- calc_de(surv3, wanted_frac)  
  
  # Calculate SPR, use survPS = 1 (no oversummer prespawning mortality) to get maximum SPR for assumptions
  # about surv1 and delta_e (assume surv2 and surv3 fixed)
  spr <- SPR_srcs(surv1, surv2, surv3, delta_e = delta_e, survPS = 1)
  
  # alpha_vec contains multipliers that are different assumptions 
  # about the productivity of the population at
  # low abundances relatiove to the slope of the replacement
  # line: alpha_vec = 2 indicates --> alpha = 2*(1/SPR)
  if (exists("alpha")) rm(alpha)
  alpha <- 1/spr*alpha_vec
  
  # For a given SPR, alpha and EQsp, we can calculate beta, for Butte Creek 
  # spring run chinook salmon
  # To calculate beta for multiple alpha values, just pass a vecter of alphas
  # to the alpha argument
  
  beta <- calc_beta(alpha = alpha, EQsp = EQsp, spr = spr)
  
  # Plot underlying stock recruitment functions 
  # Spawners vs. Age-1 
  
  plotStockRecruit(alpha = alpha, beta = beta, EQsp = EQsp, delta_e = delta_e, spr = spr, make_plot = make_plot)
  #print("works to plotStockRecruit") 
  # Spawners vs. returning spawners (3 and 4 years later from resultant cohort)
  plotStockSpEsc(alpha = alpha, beta = beta, EQsp = EQsp, delta_e = delta_e, surv1 = surv1, surv2 = surv2, surv3 = surv3, spr = spr, make_plot = make_plot)
  #print("works to plotStockSpEsc") 
  # Select PreSpawn survival scenarios
  
  if (wtrMgmt == "NoD" & climate == "A2") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "NoDiversion", "SALMODspawnFemSurvA2.txt"), sep = ","  )[,1:6])
  } 
  
  if (wtrMgmt == "NoD" & climate == "B1") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "NoDiversion", "SALMODspawnFemSurvB1.txt"), sep = ","  )[,1:6])
  } 
  
  if (wtrMgmt == "BAU" & climate == "A2") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "BAU", "SALMODspawnFemSurvA2.txt"), sep = ","  )[,1:6])     
  }
  
  if (wtrMgmt == "BAU" & climate == "B1") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "BAU", "SALMODspawnFemSurvB1.txt"), sep = ","  )[,1:6])     
  }
  
  if (wtrMgmt == "" & climate == "") survPS <- 1
  
  gcm_names <- dimnames(survPS)[[2]]
  
  # Set up population projection model ####
  
  # Set up an initial age vector
  # start from equilibrium age structure in ocean just before returning to spawn
  
  n_alpha <- length(alpha)
  n_EQsp <- length(EQsp)
  
  if (n_EQsp > 1 & n_alpha > 1) stop("n_alpha and n_EQsp cannot both be greater than 1")
  
  if (n_EQsp > 1 & n_alpha == 1)  n0 <- matrix(0, nrow = 4, ncol = n_EQsp)
  if (n_EQsp == 1 & n_alpha >= 1)  n0 <- matrix(0, nrow = 4, ncol = 1)
  
  # When using multiple EQ spawning levels, need to adjust the intial age distribution
  # based on higher or lower equilibrium abundances
  # There will be an n0 for each EQsp
  # For a single EQsp, there is only a single 
  
  if (n_EQsp > 1) {
    for (i in 1:ncol(n0)) {
      n0[1,i] <- (EQsp[i]*(1-wanted_frac))
      n0[2,i] <- n0[1,i]/(surv3*(1-delta_e))
      n0[3,i] <- n0[2,i]/surv2
      n0[4,i] <- n0[3,i]/surv1
    } 
  } else {
    n0[1] <- (EQsp*(1-wanted_frac))
    n0[2] <- n0[1]/(surv3*(1-delta_e))
    n0[3] <- n0[2]/surv2
    n0[4] <- n0[3]/surv1
  }
  
  n0 <- apply(n0,2,rev) # This flips the matix, so that age-1 is on the first row rather than the 4th row 
  
  # Storage array for population simulations
  # 1-d for each climate model (6)
  # 1-d for each productivity hypothesis (1-5) 
  # 1-d for time (90)
  # 1-d for each age (4)
  # so need a 6 X 1-5 X 90 array
  
  # Create arrays for storing model output
  if (n_EQsp > 1 & n_alpha == 1) {
    pop <- array(0, c(4, length(yrs), 6, n_EQsp))
    spawners <-  array(0, c(length(yrs), 6, n_EQsp))
  }
  
  if (n_alpha >= 1 & n_EQsp == 1) {
    pop <- array(0, c(4, length(yrs), 6, n_alpha))
    spawners <-  array(0, c(length(yrs), 6, n_alpha))
  }
  
  if (n_EQsp == 1 & n_alpha == 1) {
    pop <- array(0, c(4, length(yrs), 6, 1))
    spawners <-  array(0, c(length(yrs), 6, 1))
  }
  
  if (n_EQsp == 1 & n_alpha == 1) {prodH0 <- "oneSR"; j_len <- 1}
  if (n_EQsp > 1 & n_alpha == 1)  {prodH0 <- "EQspVary"; j_len <- n_EQsp}
  if (n_EQsp == 1 & n_alpha > 1)  {prodH0 <- "alphaVary"; j_len <- n_alpha}
  
  for (j in 1:j_len) {  
    for (c in 1:6) {
      for (t in 1:length(yrs)) {
        
        if (t == 1) { 
          # pop[a, t, c, j] # a(ge), t(years), c(limate model), j(hypothesis about EQsp or alpha)
          if (n_EQsp == 1) {
            pop[1, t, c, j] <- n0[1]
            pop[2, t, c, j] <- n0[2]
            pop[3, t, c, j] <- n0[3]
            pop[4, t, c, j] <- n0[4]
          } 
          
          if (n_EQsp > 1) {
            pop[1, t, c, j] <- n0[1,j]
            pop[2, t, c, j] <- n0[2,j]
            pop[3, t, c, j] <- n0[3,j]
            pop[4, t, c, j] <- n0[4,j]
          }
          
        } else { 
          if (length(survPS) > 1 ) {
            spawners[t-1, c, j] <- pop[3, t-1, c, j]*delta_e*survPS[t-1, c] + pop[4,t-1, c, j]*survPS[t-1, c] 
          } else {
            spawners[t-1, c, j] <- pop[3, t-1, c, j]*delta_e + pop[4,t-1, c, j]
          }
          
          if (variableSR == FALSE) {
            if (n_alpha > 1 & n_EQsp == 1) pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha[j])/(1+spawners[t-1, c, j]*beta[j])
            if (n_alpha == 1 & n_EQsp > 1) pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta[j])
            if (n_alpha == 1 & n_EQsp == 1) pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta)
          } else {
            if (n_alpha > 1 & n_EQsp == 1)  pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha[j])/(1+spawners[t-1, c, j]*beta[j])*(exp(sig_r*rnorm(1,mean=0, sd=1)))
            if (n_alpha == 1 & n_EQsp > 1)  pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta[j])*(exp(sig_r*rnorm(1,mean=0, sd=1)))
            if (n_alpha == 1 & n_EQsp == 1)  pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta)*(exp(sig_r*rnorm(1,mean=0, sd=1)))
          }
          pop[2,t, c, j] <- pop[1,t-1, c, j]*surv1
          pop[3,t, c, j] <- pop[2, t-1, c, j]*surv2
          pop[4,t, c, j] <- pop[3, t-1, c, j]*surv3*(1-delta_e)
        }
        
        if (t == length(yrs) & length(survPS) > 1 )  spawners[t] <- pop[3, t, c, j]*delta_e*survPS[t, c] + pop[4, t, c, j]*survPS[t, c] 
        if (t == length(yrs) & length(survPS) == 1 )  spawners[t] <- pop[3,t, c, j]*delta_e + pop[4, t, c, j]
      }
    }
  }
  
  if (n_alpha > 1) {
    dimnames(pop)[4] <- list(alpha)
    dimnames(pop)[3] <- list(gcm_names)
    dimnames(spawners)[3] <- list(alpha)
    dimnames(spawners)[2] <- list(gcm_names)
    popm <- melt(pop, varnames = c("age", "year", "gcm", "alpha"))
    spm <- melt(spawners, varnames = c("year", "gcm", "alpha"))
    spm_sc <- ddply(spm, .(gcm, alpha), transform, value = value/(max(value, na.rm = TRUE)))
  }
  if (n_EQsp > 1) {
    dimnames(pop)[4] <- list(EQsp)
    dimnames(pop)[3] <- list(gcm_names)
    dimnames(spawners)[3] <- list(EQsp)
    dimnames(spawners)[2] <- list(gcm_names)
    popm <- melt(pop, varnames = c("age", "year", "gcm", "EQsp"))
    
    #popm_sc <- ddply(spm, .(gcm, EQsp), transform, scale_abund = value/(max(value, na.rm = TRUE)))
    spm <- melt(spawners, varnames = c("year", "gcm", "EQsp"))
    spm_sc <- ddply(spm, .(gcm, EQsp), transform, value = value/(max(value, na.rm = TRUE)))
  }
  
  survPS <- data.frame(year = 1:90, survPS)
  survPSm <- melt(survPS, id.var = "year")
  names(survPSm) <- c("year", "gcm", "value")
  
  if (n_alpha > 1)   survPSm$alpha = 1 # 'ref'
  if (n_EQsp > 1)   survPSm$EQsp = 1 # 'ref'
  
  survPS <- survPSm[,c(1,2,4,3)]
  spm_sc <- rbind(spm_sc, survPSm)
  return(list(popm = popm, spm = spm, spm_sc = spm_sc, 
              delta_e = delta_e, spr = spr, alpha=alpha, 
              beta=beta))
  
} # end of popProjConstMarine() 


plotPopProjConst <- function(population, spawners, spawners_sc, vary) {
  if (vary == "alpha") {
    p_pop <- ggplot(population, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      #geom_point() + 
      facet_grid(age ~ alpha, scales = "free_y") +
      ylab("Abundance") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    
    print(p_pop)
    
    p_sp <- ggplot(spawners, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(alpha ~ ., scales = "free_y") + 
      ylab("Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_sp)
    
    p_spm_sc <- ggplot(spawners_sc, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(alpha ~ ., scales = "free_y") + 
      ylab("Scaled Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_spm_sc)
  }
  
  if (vary == "EQsp") {
    p_pop <- ggplot(population, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      #geom_point() + 
      facet_grid(age ~ EQsp, scales = "free_y") +
      ylab("Abundance") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_pop)
    
    p_sp <- ggplot(spawners, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(EQsp ~ ., scales = "free_y") + 
      ylab("Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_sp)
    
    p_spm_sc <- ggplot(spawners_sc, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(EQsp ~ ., scales = "free_y") + 
      ylab("Scaled Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_spm_sc)
  }
  
  if (vary == "oneEach") {
    p_pop <- ggplot(popm, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      #geom_point() + 
      ylab("Abundance") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_pop)
    
    p_sp <- ggplot(spm, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      ylab("Spawners") + 
      xlab("Year")
    print(p_sp) +
      theme_few() +
      scale_colour_few()
  }
  
} # end plotPopProjConst()


# this function counts the first four consecutive years of zero (survival or abundance)
# and then records how many years to the first or last of those consecutive year and the year
# climate <- "B1"
# mgmt <- "NoDiversion"
# 
# data <- case1A2[[3]]
calcExtTimeConst <- function(data, climate, mgmt) {
  
  library(ggthemes)
  if (names(data)[3] == "EQsp") {
    sr <- "EQsp" }  else if (names(data[3]) == "alpha") {
      sr <- "alpha" } else {
        print("names data[3] must be EQsp or alpha")
      }
  names(data)[3] <- "SR"
  datac <- dcast(data, year ~ SR + gcm) 
  extInd <- vector("integer", (ncol(datac)-1) )
  for (j in 2:ncol(datac)) {
    zeroInd <- which(datac[,j] == 0) 
    #zeroInd <- c(1,0,0,1,0,2,3,4,3,4,1)
    extRuns <- vector("integer", 0)
    k <- 1
    for (i in 1:(length(zeroInd)-3)) {
      if (length(zeroInd) >= 4) {
        if (zeroInd[i] == (zeroInd[i+1] - 1) & zeroInd[i] == (zeroInd[i+2] - 2) & 
            zeroInd[i] == (zeroInd[i+3] - 3)  ) {
          extRuns[k] <- (zeroInd[i] + 3)
        } else {
          extRuns[k] <- 90
        }
      } else {
        extRuns[k] <- 90
      }
      k <- k + 1 
    }
    if (length(extRuns > 0)) extInd[j-1] <- min(extRuns) else extInd[j-1] <- NA
  }
  
  yrs <- 2010:2099
  extYr <- vector("integer", length(extInd))
  for (i in 1:length(extInd)) {
    extYr[i] <- yrs[extInd[i]]
  }
  
  if (sr == "EQsp") {
    SRinfo <- trunc(as.numeric(gsub("([0-9]{1,5})(_)([a-z0-9]+)", "\\1" , names(datac)[-1])))
    gcm <- gsub("([0-9]{1,5})(_)([a-z0-9]+)", "\\3" , names(datac)[-1])
  } else if (sr == "alpha") {
    SRinfo <- trunc(as.numeric(gsub("([0-9\\.]{1,7})(_)([a-z0-9]+)", "\\1" , names(datac)[-1])))
    gcm <- gsub("([0-9\\.]{1,7})(_)([a-z0-9]+)", "\\3" , names(datac)[-1])
  }
  outData <- data.frame(SRinfo = SRinfo, gcm = gcm, yrs2ext = extInd, extYr = extYr)
  
  #   time2extGCM <- ggplot(data = outData, aes(x = yrs2ext, y = gcm, colour = gcm)) + 
  #     # geom_point(size = 4, alpha = .7, colour = "black") + # position = position_jitter(w = 0, h = 0.5), 
  #     geom_point(size = 3, alpha = .7, position = position_jitter(w = 0, h = 0.5)) + # 
  #     scale_color_brewer(palette = "Dark2") +
  #     theme_bw() +
  #     xlim(c(0,90)) +
  #     xlab("Years to extinction") +
  #     ggtitle(title = paste("Mean time to extinction for", climate, "&", mgmt, "scenario\nvarying", sr, sep = " "))
  #   print(time2extGCM)
  
  xTicks <- unique(outData$SRinfo)
  
  time2extSR <- ggplot(data = outData, aes(x = yrs2ext, y = SRinfo, colour = gcm)) + 
    geom_point(size = 4, colour = "black") +
    geom_point(size = 3) + 
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    xlim(c(0,90)) +
    ylab(paste(sr)) +
    xlab("Years to extinction") +
    ggtitle(paste("Mean time to extinction for", climate, "&", mgmt, "scenario\nvarying", sr, sep = " ")) +
    scale_y_continuous(breaks= xTicks)  
  print(time2extSR)
  
  return(outData)
}


# Lisa says that extinction is a QE of 4 years sans 20 fish

# # Vary alpha, A2 BAU
# data = case1A2#  
# mgmt = "BAU"
# climate = "A2"
# nSpawners = 7500
# working = YES
# getSALMODnums(case1A2) # Works



# Vary beta, A2 BAU
# data <- case2A2
# mgmt = "BAU"
# climate = "A2"
# calcQEtimeConst(case2A2, "A2", "BAU", QEthr = 100)
# working = YES (but should check)
# getSALMODnums(case2A2) #  WORKS

# # Vary alpha, B1 No Diversion
# data = case1B1
# mgmt = "NoDiversion"
# climate = "B1"
# calcQEtimeConst(case1B1, "B1", "No Diversion", QEthr = 100)
# working = NO 
# getSALMODnums(case1B1) # WORKS


# # Vary alpha, B1 No Diversion
# data = case2B1
# mgmt = "NoDiversion"
# climate = "B1"
# calcQEtimeConst(case2B1, "B1", "No Diversion", QEthr = 100)
# working = Yes
# getSALMODnums(case2B1) # WORKS

# Pick a QE threshold
# QEthr <- 2 # 100, 20, 0

getSALMODnums <- function(Sdata) {
  # Function to get the number of females surviving to spawn in the fall calculated by SALMOD
  # to compare SALMOD time to QE with those from the deterministic population projection model
  # Just get the SALMOD prespawn survival output (stored in the [[3]] list from deterministic population projection model)
  names(Sdata[[3]])[3] <- "SR"
  SALMOD_surv  <-  dcast(subset(Sdata[[3]], SR == 1)[,c(1,2,4)], year ~ gcm)
  SALMOD_num <- cbind(SALMOD_surv[,1], SALMOD_surv[,2:ncol(SALMOD_surv)] * 7568) # 7568 is number of spawners used by LCT and CM
  return(SALMOD_num)
} # end getSALMODnums()


customFR2ts <- function(N, # number of time steps
                        reps,
                        r_seed = 1,
                        amp,
                        freqs = 200) {
  # To create white noise (equal variance at all frequencies) 
  # generate a sine waves of frequencies from >0 to up to the maximum frequency 0.5, 
  # then randomly assign each sine wave a phase picked from a 
  # uniform distribution between 0 and 2 pi, 
  # then add them together and normalize or scale to desired time series variance.
  
  # freqs: cap the low frequency and limit the number of sine waves
  
  t <- 1:N # time index
  # f <- (1:(N/2))/N # frequencies from 1/N to 0.5
  f <- seq(0.00001, 0.5, length.out = freqs)# frequencies from 1/N to 0.5
  tsReps <- matrix(NA, nrow = N, ncol = reps)
  set.seed(r_seed)
  for (h in 1:reps) { #for each col in tsReps (replicate of noise time series)
    
    rtheta <- runif(length(f), 0, 2*pi) # the random phases to apply to each frequency
    
    ts <- matrix(NA, nrow = N, ncol = length(f) ) #empty matrix, each col is timeseries at a frequency
    
    for (i in 1:length(f)) { #for each freq (going across columns)
      if (length(amp) == 1) { #if the period 
        ts[,i] <- amp * cos(2*pi*f[i]*t - rtheta[i])
      } else {
        ts[,i] <- amp[i] * cos(2*pi*f[i]*t - rtheta[i])         
      }
    }
    
    noise <- rowSums(ts) # add up the curves
    noise <- (noise - mean(noise, na.rm = TRUE) )/ sd(noise, na.rm = TRUE) # mean = 0, sd of 1
    tsReps[,h] <- noise
  }
  return(tsReps)
  
} # end of customFR2ts

# a wrapper for customFR2ts ####
noiseWrapper <- function(N, reps, r_seed = 1) {
  # create "reps" number of white noise time series of length N using random sine wave approach
  # with selected frequency contents
  white_n <- customFR2ts(N = N, # number of time steps
                         reps = reps,
                         r_seed = r_seed,
                         amp = mk_white(N)) # mean = 0, variance/sd = 1
  
  # Bandpass period 3-4 (frequencies 1/4 to 1/3) in white noise signals.
  
  rsin_34_n <- customFR2ts(N = N, # number of time steps
                           reps = reps,
                           r_seed = r_seed,
                           amp = mk_rsin(N, highF=1/3, lowF=1/4)) # mean = 0, variance/sd = 1
  
  # Bandreject period 3-4 (frequencies 1/4 to 1/3) in white noise signals.
  
  rsin_34_reject_n <- customFR2ts(N = N,
                                  reps = reps,
                                  r_seed = r_seed,
                                  amp = mk_rsin_reject(N, lowF=1/4, highF=1/3) )
  # Band-pass greater than period 4 (frequencies lower than 0.25)
  rsin_gt4_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=0, highF=1/4) )
  
  # Band-pass less than period 3 (frequencies higher than 0.33)
  rsin_lt3_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=1/3, highF=1/2) )
  
  # Band-pass greater than period 3 (frequencies lower than 0.33)
  rsin_gt3_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=0, highF=1/3) )
  
  # Band-pass less than period 4 (frequencies higher than 0.25)
  rsin_lt4_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=1/4, highF=1/2) )
  
  # Band-pass greater than period 10 (frequencies lower than 0.1)
  rsin_gt10_n <- customFR2ts(N = N,
                             reps = reps,
                             r_seed = r_seed,
                             amp = mk_rsin(N, lowF=0, highF=1/10) )  
  
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
                    nose_34_reject = rsin_34_reject_n,
                    noise_gt4 = rsin_gt4_n,
                    noise_lt3 = rsin_lt3_n,
                    noise_gt3 = rsin_gt3_n,
                    noise_lt4 = rsin_lt4_n,
                    noise_gt10 = rsin_gt10_n,
                    noise_34gt10 = rsin_34_gt10_n)
  
  # rm individual sets of noise
  rm(white_n, rsin_34_n, rsin_34_reject_n, rsin_gt4_n, rsin_lt3_n, rsin_gt3_n, rsin_lt4_n, rsin_gt10_n, rsin_34_gt10_n)
  return(noiseList)
} # end of noiseWrapper <- function(N, reps, r_seed = 1)

# plot the mean frequency response of multiple spectra generated by the same noise ####

plotMeanFreqR <- function(dataMat, N) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = ncol(dataMat))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(dataMat)) {
    ifelse(all(dataMat[,i] == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           spcMean[,i] <- spec.pgram(scale(dataMat[,i]), plot = F, c(m,m))$spec )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
  
  plot(freq, mean_spc, type = "n", lwd = 3.7, lty = 1, col = "black", xlab = "Frequency", ylab = "Relative Magnitude", las = 1, 
       ylim = c(0, ceiling(max(spc90, na.rm = TRUE))))
  
  lines(freq, spc90, lwd = 2)
  lines(freq, spc10, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc10, rev(spc90)), col = "lightgrey")
  
  lines(freq, spc75, lwd = 2)
  lines(freq, spc25, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc25, rev(spc75)), col = "darkgrey")
  
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
  lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
  lines(freq, mean_spc, lwd = 2, lty = 2)
  
  box(lwd=2)
}

# plot mean Frequency response for data.tables

plotMeanFreqR_DT <- function(dataTable, N, surv, scale = "CV", yaxis_lim = c(0, ceiling(max(spc90, na.rm = TRUE))) ) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(spcMean)) {
    ifelse(all(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                         j = white] == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                       j = white]/mean(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                                                 j = white]), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                             j = white]), plot = F )$spec ) )
    #            spcMean[,i] <- spec.pgram(scale(dataTable[i = reps_c == i & meanPS_c == surv,
    #                                                      j = white]), plot = F, c(m,m))$spec )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
  
  #if (!exists("yaxis_lim")) yaxis_lim <- c(0, ceiling(max(spc90, na.rm = TRUE)))
  
  plot(freq, mean_spc, type = "n", lwd = 3.7, lty = 1, col = "black", xlab = "Frequency", ylab = "Relative Magnitude", las = 1, 
       ylim = yaxis_lim)
  
  lines(freq, spc90, lwd = 2)
  lines(freq, spc10, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc10, rev(spc90)), col = "lightgrey")
  
  lines(freq, spc75, lwd = 2)
  lines(freq, spc25, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc25, rev(spc75)), col = "darkgrey")
  
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
  lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
  lines(freq, mean_spc, lwd = 2, lty = 2)
  
  box(lwd=2)
}

# Plot only the mean Frequency Response, and allow adding addition Frequency Response plots (lines)

# plotMeanFR_DTmanyAR <- function(dataTable, N, surv, scale = "CV", AR_col, yaxis_lim) {
#   # spectral frequency 
#   if (trunc(sqrt(N)) %% 2 == 0) {
#     m <- trunc(sqrt(N)) + 1 
#   } else {
#     m <- trunc(sqrt(N))
#   }
#   
#   spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
#   
#   freq <- (1:(N/2))/N
#   for (i in 1:ncol(spcMean)) {
#     
#     ts <- as.ts(droplevels(subset(dataTable, reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001, select = AR_col) ))
#     
#     ifelse(all(ts == 0), 
#            spcMean[,i] <- rep(NA, times = N/2),
#            ifelse(scale == "CV", 
#                   spcMean[,i] <- periodogram(ts/mean(ts), plot = F )$spec, 
#                   spcMean[,i] <- periodogram(scale(ts), plot = F )$spec ) )
#   }
#   mean_spc <- rowMeans(spcMean, na.rm = TRUE)
#   
#   spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
#   spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
#   spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
#   spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
#   
#   if ( !exists("yaxis_lim") ) ifelse(scale == "CV", yaxis_lim <- c(0,1.2*max(mean_spc[which(freq > 0.2)])), yaxis_lim <- c(0, 20))
#   
#   plot(freq, mean_spc, type = "n", xlab = "Frequency", ylab = "Relative Magnitude", las = 1, ylim = yaxis_lim)
#   lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
#   lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
#   lines(freq, mean_spc, lwd = 2, lty = 2)
#   
#   box(lwd=2)
# }

plotMeanFR_DTmany <- function(dataTable, N, surv, 
                              scale = "CV", yaxis_lim) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  freq <- (1:(N/2))/N
  
  for (i in 1:ncol(spcMean)) {
    ifelse(all(dataTable[i = reps_c == i & meanPS_c == surv,
                         j = white] == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(dataTable[i = reps_c == i & meanPS_c == surv,
                                                       j = white]/mean(dataTable[i = reps_c == i & meanPS_c == surv,
                                                                                 j = white]), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(dataTable[i = reps_c == i & meanPS_c == surv,
                                                             j = white]), plot = F )$spec ) )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  if ( !exists("yaxis_lim") ) ifelse(scale == "CV", yaxis_lim <- c(0,1.2*max(mean_spc[which(freq > 0.2)])), yaxis_lim <- c(0, 20))
  plot(freq, mean_spc, type = "n", ylab = "Relative Magnitude", 
       las = 1, ylim = yaxis_lim,
       xlab = expression('Frequency' ~ y^-1))
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
  lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
  lines(freq, mean_spc, lwd = 2, lty = 2)
  
  box(lwd=2)
}

linesMeanFR_DTmany <- function(dataTable, N, surv, scale = "CV", line_color) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(spcMean)) {
    ifelse(all(dataTable[i = reps_c == i & meanPS_c == surv,
                         j = white] == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(dataTable[i = reps_c == i & meanPS_c == surv,
                                                       j = white]/mean(dataTable[i = reps_c == i & meanPS_c == surv,
                                                                                 j = white]), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(dataTable[i = reps_c == i & meanPS_c == surv,
                                                             j = white]), plot = F )$spec ) )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
  
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = line_color)
  
  box(lwd=2)
}


calcQEtimeConstPers <- function(data, climate, mgmt, QEthr) {
  
  library(ggthemes)
  
  sr <- "alpha" 
  names(data[[2]])[3] <- "SR"
  datac <- dcast(data[[2]], year ~ SR + gcm) # numbers from deterministic population model 
  salmod_nums <- getSALMODnums(data) # SALMOD Numbers
  
  # Calculate time to Quasi-extinction for the deterministic population model
  extInd <- vector("integer", (ncol(datac)-1) )
  for (j in 2:ncol(datac)) {
    qeInd <- which(datac[,j] <= QEthr) 
    extRuns <- vector("integer", 0)
    k <- 1
    if (length(qeInd) >= 4) {
      for (i in 1:(length(qeInd)-3)) {
        if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
            qeInd[i] == (qeInd[i+3] - 3)  ) {
          extRuns[k] <- (qeInd[i] + 3)
          k <- k + 1
        }  
      }
    }else {
      extRuns[k] <- 90
      k <- k + 1
    }
    if (length(extRuns > 0)) extInd[j-1] <- min(extRuns) else extInd[j-1] <- 90
  }
  # Calculate time to Quasi-extinction for the SALMOD data
  SALMODextInd <- vector("integer", (ncol(salmod_nums)-1) )
  for (j in 2:ncol(salmod_nums)) {
    qeInd <- which(salmod_nums[,j] <= QEthr) 
    SALMODextRuns <- vector("integer", 0)
    k <- 1
    if (length(qeInd) >= 4) {
      for (i in 1:(length(qeInd)-3)) {
        if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
            qeInd[i] == (qeInd[i+3] - 3)  ) {
          SALMODextRuns[k] <- (qeInd[i] + 3)
          k <- k + 1
        }  
      } 
    } else {
      SALMODextRuns[k] <- 90
      k <- k + 1
    }
    if (length(SALMODextRuns > 0)) SALMODextInd[j-1] <- min(SALMODextRuns) else SALMODextInd[j-1] <- 90
  }
  
  yrs <- 2010:2099
  # Vector of extinction year for each climate model/stock recruit (pop dy) scenario
  extYr <- vector("integer", length(extInd))
  for (i in 1:length(extInd)) {
    extYr[i] <- yrs[extInd[i]]
  }
  # Vector of extinction year for each climate model/stock recruit (pop dy) scenario
  SALMODextYr <- vector("integer", length(SALMODextInd))
  for (i in 1:length(SALMODextInd)) {
    SALMODextYr[i] <- yrs[SALMODextInd[i]]
  }
  
  SRinfo <- trunc(as.numeric(gsub("([0-9\\.]{1,7})(_)([a-z0-9]+)", "\\1" , names(datac)[-1])))
  gcm <- gsub("([0-9\\.]{1,9})(_)([a-z0-9]+)", "\\3" , names(datac)[-1])
  
  outData <- data.frame(SRinfo = SRinfo*data[[5]], gcm = gcm, yrs2ext = extInd, extYr = extYr)
  SALMOD <- data.frame(gcm = unique(gcm), yrs2ext = SALMODextInd, extYr = SALMODextYr)
  xTicks <- 1:10#seq(0,1, by = 0.1) # unique(outData$SRinfo) 
  
  time2extSR <- ggplot(data = outData, aes(y = yrs2ext, x = SRinfo, colour = gcm)) + 
    geom_line(linetype = 1, size = 1) +
    geom_hline(aes(yintercept=90), size = 2) +
    scale_color_brewer(palette = "Dark2", name = "") +
    theme_classic() +
    xlab(expression(paste(alpha,  " x SPR"))) +
    ylab("Years to QE") +
    labs(title = paste(climate, " & ", mgmt, ", QE = ", QEthr, sep = "")) +
    scale_x_continuous(breaks = xTicks) +
    ylim(c(0,90)) + 
    geom_hline(data = SALMOD, aes(yintercept = yrs2ext, colour = gcm), linetype = 2,
               size = 1, alpha = 0.8) 
  
  print(time2extSR)
  
  
  extYrSR <- ggplot(data = outData, aes(y = extYr, x = SRinfo, colour = gcm)) + 
    geom_line(linetype = 1, size = 1) +
    scale_color_brewer(palette = "Dark2", name = "") +
    theme_classic() +
    xlab(expression(paste(alpha,  " x SPR"))) +
    ylab("Year of QE") +
    labs(title = paste(climate, " & ", mgmt, ", QE = ", QEthr, sep = "")) +
    scale_x_continuous(breaks= xTicks) +
    ylim(c(2010,2100)) +
    geom_hline(data = SALMOD, aes(yintercept = extYr, colour = gcm), linetype = 2,
               size = 1, alpha = 0.8) 
  print(extYrSR)
  
  return(outData) 
  
} # end calcQEtimeConst 

calcQEtimeConstEQvQE <- function(data, climate, mgmt) {
  
  # data = spawner output for many runs for each climate and mgmt e.g.,  EQsp <- seq(500, 50000, by = 500)
  # for each climate and mgmt loop over each EQsp and QEthr 
  
  QEthr <- seq(0, 1000, by = 10)
  sr <- "EQsp" 
  names(data[[2]])[3] <- "SR"
  datac <- dcast(data[[2]], year ~ SR + gcm) # numbers from deterministic population model 
  
  # Calculate time to Quasi-extinction for the deterministic population model
  extInd <- extYr <- matrix(NA, nrow = length(QEthr), ncol = (ncol(datac)-1))
  
  for (q in 1:length(QEthr)) {
    for (j in 2:ncol(datac)) {
      qeInd <- which(datac[,j] <= QEthr[q]) 
      extRuns <- vector("integer", 0)
      k <- 1
      if (length(qeInd) >= 4) {
        for (i in 1:(length(qeInd)-3)) {
          if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
              qeInd[i] == (qeInd[i+3] - 3)  ) {
            extRuns[k] <- (qeInd[i] + 3)
            k <- k + 1
          }  
        }
      }else {
        extRuns[k] <- 90
        k <- k + 1
      }
      if (length(extRuns > 0)) extInd[q, j-1] <- min(extRuns) else extInd[q, j-1] <- 90
    }
  } # end of q loop
  
  yrs <- 2010:2099
  # Matrix of extinction year for each climate model/stock recruit (pop dy) scenario
  for (i in 1:nrow(extInd)) {
    for (j in 1:ncol(extInd)) {
      extYr[i,j] <- yrs[extInd[i,j]]
    }
  }
  
  # DATA should look like this:
  #   GCM  |  SRinfo  | QEthr | extInd  | extYr
  #   a    |   1000   | 10    |   50    |  ...
  #   a    |   1000   | 20    |   50    |  ...
  #   a    |   1000   | 50    |   50    |  ...
  #   a    |   2000   | 10    |   50    |  ...
  #   a    |   2000   | 20    |   50    |  ...
  #   a    |   2000   | 50    |   50    |  ...
  #   b    |   1000   | 10    |   50    |  ...
  #   b    |   1000   | 20    |   50    |  ...
  #   b    |   1000   | 50    |   50    |  ...
  #   b    |   2000   | 10    |   50    |  ...
  #   b    |   2000   | 20    |   50    |  ...
  #   b    |   2000   | 50    |   50    |  ...
  # and have GCM * SRinfo * QEthr rows
  
  # wrangle extInd data into above format
  extIndDF <- as.data.frame(extInd)
  extYrDF <- as.data.frame(extYr)
  names(extIndDF) <- names(extYrDF) <- names(datac)[-1]
  extIndDF <- cbind(QEthr = QEthr, extIndDF)
  extYrDF <- cbind(QEthr = QEthr, extYrDF)
  extIndDF_m <- melt(extIndDF, id.vars = "QEthr")
  extYrDF_m <- melt(extYrDF, id.vars = "QEthr")
  
  SRinfo <- trunc(as.numeric(gsub("([0-9]{1,16})(_)([a-z0-9]+)", "\\1" ,  extIndDF_m$variable)))#names(datac)[-1])))
  gcm <- gsub("([0-9.]{1,16})(_)([a-z0-9]+)", "\\3" , extIndDF_m$variable)#names(datac)[-1])
  
  outData <- data.frame(gcm = gcm, SRinfo = as.numeric(SRinfo), QEthr = extIndDF_m$QEthr, yrs2ext = extIndDF_m$value, extYr = extYrDF_m$value)
  
  time2extSR <- ggplot(data = outData, aes(x = QEthr, y = SRinfo, z = yrs2ext)) +
    geom_tile(aes(fill = yrs2ext)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "black", midpoint = 45) +
    #stat_contour(aes(colour = ..level..), binwidth = 10) +
    facet_grid(gcm~.) +
    xlab("Quasi-extinction threshold (numbers)") +
    ylab("Equilibrium abundance of spawners") +
    theme_classic() +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) 
  
  print(time2extSR)
  
  extYrSR <- ggplot(data = outData, aes(x = QEthr, y = SRinfo, z = extYr)) + 
    geom_tile(aes(fill = extYr)) +
    scale_fill_gradient2(low="red", mid = "white", high = "black", midpoint = 2055) +
    facet_grid(gcm~.) +
    xlab("Quasi-extinction threshold (numbers)") +
    ylab("Equilibrium abundance of spawners") +
    theme_classic() +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) 
  print(extYrSR)
  
  return(outData)
  
  
} # end calcQEtimeConstEQvQE



QE_counter <- function(data, qeLev = 100, run_length = 4) {
  
  extInd <- vector("integer", ncol(data))
  for (h in 1:ncol(data)) {  
    qeInd <- which(data[,h] <= qeLev)
    extRuns <- vector("integer", 0)
    k <- 1
    if (length(qeInd) >= run_length) {
      for (i in 1:(length(qeInd)-3)) {
        if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
            qeInd[i] == (qeInd[i+3] - 3)  ) {
          extRuns[k] <- (qeInd[i] + 3)
          k <- k + 1
        }  
      }
    } else {
      extRuns[k] <- (nrow(data) + 1) 
      k <- k+1
    }
    if(length(extRuns > 0)) extInd[h] <- min(extRuns) else extInd[h] <- NA
  }
  return(extInd)
} # end QE_count

# A better QE counter ####
# courtesy Jaime Ashander

JA_consec <- function(z, qeLev = 100, run_length = 4) {
  i.low <- which(z <= qeLev)
  
  if (length(i.low) <= run_length) return(as.integer(NA))
  
  i.last <- i.low[1]
  counter <- 1
  for( i in i.low[-1]){
    
    if(counter == run_length)  return(as.integer(i.last)) #  need the double to ensure that data.table gets what it expects
    
    if(i-1 == i.last)  counter <- counter + 1 else counter <- 1
    
    if(i == i.low[length(i.low)] & counter < run_length) i.last <- NA else i.last <- i
  }
  if (is.null(class(i.last))) print("error here")
  return(as.integer(i.last))
}
#if (length(which(z < qeLev)) < run_length)  return(length(z)+1) else i.low <- which(z < 100)


# Stochastic population projection models ####

# create random time series
# white noise
# mk_white <- function(N) {
#   return(rep(1/N, N/2))
# }

mk_white <- function(freq) {
  # return(rep(1/freq, freq/2))
  return(rep(1/freq, freq))
}


mk_1_over_f_beta <- function(N, beta, freq = 200) {
  # fs <- seq(0, 0.5, length.out = (N/2))
  # fs[1] <- 0.00001 # avoid infinity 
  fs <- seq(0.00001, 0.5, length.out = freq)
  one_over_fb <- 1/fs^beta
  return(one_over_fb)
}

# create "band-pass" style time series
mk_rsin <- function(N, lowF, highF) {
  temp <- (1:(N/2))/(N)
  ns <- length(which(temp <= highF & temp >= lowF))
  goodInd <- which(temp <= highF & temp >= lowF)
  temp[-goodInd] <- 0
  temp[goodInd] <- 1/ns
  return(temp)  }

mk_rsin2 <- function(N, lowF, highF, freq) {
  temp <- seq(0.00001, 0.5, length.out = freq)
  ns <- length(which(temp <= highF & temp >= lowF))
  goodInd <- which(temp <= highF & temp >= lowF)
  temp[-goodInd] <- 0
  temp[goodInd] <- 1/ns
  return(temp)  }

# create "band-reject" style time series

mk_rsin_reject <- function(N, lowF, highF) {
  # temp <- (1:(N/2))/(N)
  temp <- seq(0.00001, 0.5, length.out = freq)
  ns <- length(which(temp >= highF & temp <= lowF))
  # print(ns)
  goodInd <- which(temp >= highF & temp <= lowF)
  temp[-goodInd] <- 0
  temp[goodInd] <- 1/ns
  #print(temp)
  return(temp)  
}

# Wavelet filtering of time series to do "band-pass" style filtering
# do band pass filter on the wavelet decomposition of a white noise series.
# This reconstructs the signal using only the pd 3 to 4 components of the
# wavelet power spectrum (T&C 1998 equation 11)
waveFlt34 <- function(white_n) {
  # Wavelet filtering pd 3-4 signal from white noise signal above
  # White noise filtering out pd 3-4 variability from the wavelet power spectrum 
  # vs generating period 3-4 noise using the random sine wave method
  
  # Create a "time" vector
  N <- nrow(white_n)
  x <- 1:N
  
  # I copied/translated Flora's matlab code (Will, Matt, Lauren, T&C) for J1
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1
  # Torrence and Compo 1998 equation 11 constants for reconstruction of time series for
  # Morlet mother wavelets using a delta function 
  # Morlet wavelet constants
  C_delta <- 0.776 #  The factor C_delta comes from the reconstruction of a
  # delta function from its wavelet transform using the function Psi_0(eta)
  Psi_0 <- pi^(-0.25) # removes the energy scaling
  dj <- 0.01
  dt <- 1
  
  reCon <- (dj*dt^(0.5))/(C_delta*Psi_0) # Constant for reconstruction of wavelets
  
  library(biwavelet) # R package for continuous wavelet analysis based on T & C 1998
  
  wave34_n <- matrix(NA, nrow = nrow(white_n), ncol = ncol(white_n))
  
  for (i in 1:ncol(white_n)) {
    white.wt <- wt(cbind(x, white_n[,i]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
    pd34ind <- which(white.wt$period >= 3 & white.wt$period <= 4)
    # new_y <- as.numeric(scale(colMeans(white.wt$power.corr[pd34ind,])))
    # new.wt <- wt(cbind(x,new_y), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
    
    sumWave <- matrix(NA, nrow = nrow(white.wt$wave),
                      ncol = ncol(white.wt$wave))
    
    for (j in 1:nrow(white.wt$wave)) {
      sumWave[j,] <- Re(white.wt$wave[j,])/(white.wt$scale[j]^(0.5))
      # scaling by s^(1/2) converts the wavelet transform to an energy density
    }
    
    # reCon_y <- reCon*colSums(sumWave) # reconstruct the original white noise signal
    sumWave34 <- sumWave[pd34ind,]
    temp <- reCon*colSums(sumWave34)
    wave34_n[,i] <- (temp - mean(temp))/sd(temp)
  }
  return(wave34_n)
} # end of waveFlt34()


popSimPSvary <- function(rand_surv, surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale) {
  # popSimPSvary() is a wrapper for PopProjPSvar().
  # PopProjPSvar() runs a single stochasitic run of the population model with time varying survPS
  # for this specified model parameters: surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale
  # popSimPSvary() just creates output containers (outPop: contains abundance at age for each year;
  # and outSpawn: contains the number of age_3 and age_4 critters surviving to spawn) to collect the output
  # of the replicate runs [for each columon of the rand_surv input]
  outPop <- array(NA, c(reps, 4, nrow(rand_surv))) #abundance at age for each year
  outSpawn <- matrix(NA, ncol = reps, nrow = nrow(rand_surv)) #number of age 3 and 4 surviving to spawn
  
  for (i in 1:ncol(rand_surv)) {
    outList <- PopProjPSvar(survPS = rand_surv[,i], 
                            surv1 = surv1, 
                            surv2 = surv2, 
                            surv3 = surv3,
                            EQsp = EQsp, 
                            wanted_frac = wanted_frac,
                            alpha_scale = alpha_scale)
    outPop[i,,] <- outList[[1]]
    outSpawn[,i] <- outList[[2]]
  }
  return(list(outPop, outSpawn))
} # end of popSimPSvary()

PopProjPSvar <- function(survPS, surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale) {
  # A density dependent (Beverton-Holt), age-structure populaiton model
  # based on previous Botsford lab work (e.g., Worden et al. 2010, ) and modified 
  # based on demographic characteristics of Butte Creek Spring Run Chinook salmon 
  # 
  # survPS = either constant or variable 
  # surv1, surv2 and surv3 = constant. surv1 ~ early ocean survival (0.01-0.05). surv2 == surv3 (0.8 high)
  # EQsp = the number of age-3 and age-4 females returning to spawn in Butte Creek, susceptible to prespawning 
  # mortality
  # wanted_frac = the ratio of age-3 spawners to total spawners (age-3 + age-4), used to determine the fraction 
  # of age-3 fish in the ocean returning to spawn each year
  # alpha_scale = a multiplier of the 1/spr to get the slope of the Beverton-Holt stock recruitment curve
  #
  # Use popSimPSvary() to run multiple simulations and collect output for analysis
  
  # PopProjPSvar() runs a single stochasitic run of the population model with time varying survPS
  # for this specified model parameters: surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale
  
  
  # May 20: check that function arguments are of length == 1
  if ( length(surv1) != 1 | length(surv2) != 1 | length(surv3) != 1 | 
       length(EQsp) != 1 | length(wanted_frac) != 1 | length(alpha_scale) != 1) stop("Function arguments survPS, surv1, surv2, 
                                                                                     surv3, EQsp, wanted_frac, 
                                                                                     be length 1! Doofus...")
  
  delta_e <- calc_de(wanted_frac, surv3)
  # assuming survPS = 1 produces maximum spr, 
  # hence least steep replacement line
  if (length(survPS) > 1) { 
    spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS=1) 
  } else {
    spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS=survPS) 
  } 
  alpha <- 1/spr*alpha_scale 
  beta <- calc_beta(alpha, EQsp, spr)  
  
  # Set up an initial age vector
  # start from equilibrium age structure in ocean just before returning to spawn
  
  n0 <- numeric(4) #matrix(0, nrow = 4, ncol = 1)
  
  # There will be an n0 for a given EQsp, 
  # but have to calculate from oldest age
  # to youngest age
  
  #   Test code to get initial abundance at age vector from reduced EQ based on reduced meanPS level
  #   ---
  #   spr_loweq <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS = PSmean) 
  #   
  #   init_eq <- (alpha*spr_loweq - 1)/beta # equilibrium number of spawners
  #   n0[1] <- (init_eq*(1-wanted_frac)) # age-4
  # ---
  n0[1] <- (EQsp*(1-wanted_frac)) # age-4
  n0[2] <- n0[1]/(surv3*(1-delta_e))
  n0[3] <- n0[2]/surv2
  n0[4] <- n0[3]/surv1
  
  n0 <- rev(n0) # This flips the n0 vector, so that age-1 is now in the 
  # first position rather than the 4th position 
  # Storage array for population simulations
  # 1-d for each age (4)
  # 1-d for time (nyrs)
  # so need a 4 X nyrs dimension MATRIX for pop and 
  # a VECTOR nyrs long for spawners
  
  # Create arrays for storing model output
  #burn_in <- trunc(length(envF)/3)
  nyrs <- length(survPS) #burn_in + length(envF)
  #envF <- c(sample(envF, size = trunc(length(envF)/3), replace = TRUE), envF)
  pop <- matrix(0, nrow = length(n0), ncol = nyrs) #array(0, c(4, n,  reps))
  sp <-  numeric(nyrs) # array(0, c(n, reps))
  
  for (t in 1:nyrs) {
    if (t == 1) {   
      pop[1, t] <- n0[1]
      pop[2, t] <- n0[2]
      pop[3, t] <- n0[3]
      pop[4, t] <- n0[4]
    } else {
      sp[t-1] <- pop[3, t-1]*delta_e*(survPS[t-1]) + pop[4,t-1]*(survPS[t-1]) 
      pop[1,t] <- (sp[t-1]*alpha)/(1+sp[t-1]*beta)  
      pop[2,t] <- pop[1, t-1]*surv1
      pop[3,t] <- pop[2, t-1]*surv2
      pop[4,t] <- pop[3, t-1]*surv3*(1-delta_e)
    }
    sp[t] <- pop[3, t]*delta_e*(survPS[t]) + pop[4, t]*(survPS[t])   
  }  
  
  return(list(pop = pop, sp = sp))
  
} # end PopProjPSvar()

# scale_surv_dat <- function(dat,
#                            surv_mean,
#                            surv_range) {
#   # scale noise series to specified mean and range
#   # shift to min of 0
#   min_dat <- apply(dat, 2, min)
#   shift_dat <- sapply( 1:ncol(dat), function(x) dat[,x] + (-1*min_dat[x]) )
#   # scale to max of 1
#   max_dat <- apply(shift_dat, 2, max)
#   dat_01 <- sapply(1:ncol(shift_dat), function(x) shift_dat[,x]/max_dat[x])
#   # scale and shift to specified mean and range
#   out_dat <- (surv_mean/2) + (dat_01 * surv_range)
#   return(out_dat)
# }

# make_surv_mat <- function(noise_dat,
#                           mean_surv,
#                           sd_surv,
#                           sim_len,
#                           burn_in) {
#   surv <- matrix(NA, nrow = nrow(noise_dat), ncol = ncol(noise_dat))
#   surv[1:burn_in, ] <- noise_dat[1:burn_in, ] * 0.01 + mean_surv
#   surv[(burn_in + 1):(sim_len+burn_in), ] <- noise_dat[(burn_in + 1):(sim_len+burn_in), ] * sd_surv + mean_surv
#   surv[surv > 1] <- 1
#   surv[surv < 0] <- 0
#   return(surv)
# }

# make_surv_mat() generates matrix that is ts length + burn in time rows long
# and reps columns wide (1124 x 1000). Could be white noise, period 3-4 etc.
# 
make_surv_mat <- function(noise_dat,
                          mean_surv,
                          sd_surv,
                          sim_len,
                          phasein_len,
                          burn_in) {
  surv <- matrix(NA, nrow = nrow(noise_dat), ncol = ncol(noise_dat))
  surv[1:burn_in, ] <- 1 # 0.8 # noise_dat[1:burn_in, ] * sd_surv + 0.8
  surv[(burn_in + 1):(phasein_len + burn_in), ] <-
    noise_dat[(burn_in + 1):(phasein_len + burn_in), ] * sd_surv + mean_surv # seq(0.01, sd_surv, length.out = phasein_len) + mean_surv
  surv[(burn_in + phasein_len + 1):(sim_len + phasein_len + burn_in), ] <-
    noise_dat[(burn_in + phasein_len + 1):(sim_len + burn_in + phasein_len), ] * sd_surv + mean_surv
  surv[surv > 1] <- 1
  surv[surv < 0] <- 0
  # surv <- apply(surv, 2, detrend)
  return(surv)
}

# make_surv_mat_range <- function(noise_dat, 
#                           mean_surv, 
#                           range_surv,
#                           sim_len = 1024) {
#   # n_noise <- length(noise_dat)
#   # surv_list <- vector("list", length = n_noise)
#   # for (i in 1:n_noise) {
#     surv <- scale_surv_dat(dat = noise_dat,
#                            surv_mean = mean_surv,
#                            surv_range = range_surv)
#     surv[surv > 1] <- 1
#     surv[surv < 0] <- 0
#   #   surv_list[[i]] <- surv
#   # }
#   return(surv)
# }

##### beta distribution approach #####

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# a fucntion for generation correlated RV from two normal distributions
# http://www.sitmo.com/article/generating-correlated-random-numbers/
# http://r.789695.n4.nabble.com/generate-two-sets-of-random-numbers-that-are-correlated-td3736161.html
corrNorm <- function(n, rho, seed) { 
  
  if ( rho < 0 | rho > 1 ) stop("rho must be between 0 and 1")
  set.seed(seed)
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1) 
  Z <- rho*X1+sqrt(1-rho^2)*X2
  return(list(X1, X2, Z) )
  
} 

# function for ar-1 noise translated from Matt Holland's Matlab code

ar1rand <- function(n, mu, sigma, phi, seed) {
  # R version of Matt Hollands Matlab function for making
  # AR(1) noise
  # calculate AR(1) model parameters
  c <- mu*(1-phi)
  sigma_xi <- sqrt(sigma^2*(1-phi^2))
  set.seed(seed)
  xi <- sigma_xi * rnorm(n, 0, 1)
  # initialize y and define first element
  y <- numeric(n)
  y[1] <- c + phi*mu + xi[1]
  # iterate over the remaining samples of xi, assigning subsequent
  # values of y
  for (i in 2:n) {
    y[i] <- c + phi * y[i-1] + xi[i]
  }
  return(y)
}

ar1_redden <- function(ts, mu, sigma, phi) {
  # Variant of R version of Matt Hollands Matlab function for making
  # AR(1) noise
  # This one reddens an existing time series
  
  # calculate AR(1) model parameters
  n <- length(ts)
  n_phi <- length(phi)
  c <- mu*(1-phi)
  sigma_xi <- sqrt(sigma^2*(1-phi^2))
  
  # initialize y and define first element
  y <- xi <- matrix(0, nrow = n, ncol = n_phi)
  
  for (k in 1:n_phi) {
    xi[,k] <- sigma_xi[k] * ts#rnorm(n, 0, 1)
    y[1,k] <- c[k] + phi[k]*mu + xi[1,k]
    # iterate over the remaining samples of xi, assigning subsequent
    # values of y
    for (i in 2:n) {
      y[i,k] <- c[k] + phi[k] * y[i-1,k] + xi[i,k]
    }  
  }
  colnames(y) <- paste("R=", phi, sep="")
  return(y)
}  

# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
SPR_srcs <- function(surv1, surv2, surv3, delta_e, survPS = 1) {
  SPR <- surv1*surv2*delta_e*survPS + surv1*surv2*surv3*(1-delta_e)*survPS
  return(SPR)
} # end bracket for SPR_scrs

# plot how SPR changes with early ocean survival and fraction of age-3 spawners
plotSPR_de_surv1 <- function (sprDF) {
  ggplot(sprDF, aes(x = delta_e, y = SPR, colour = factor(surv1))) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("black","grey20", "grey40", "grey70"), name = "Early Ocean \nSurvival") +
    ylab("Spawners Per Recruit") +
    xlab("Fraction of spawners returning at age-3") +
    xlim(c(0, 1)) +
    ylim(c(0, 0.08)) +
    theme_bw() +
    theme(panel.grid.major = theme_line(colour = NA), panel.grid.minor = theme_line(colour = NA)) 
}

Zuur_acm_univar <- function(vector, path4file) {
  
  df <- data.frame(years = 2010:2099, surv = vector, period=gl(n = 3, k = 30, labels = c("Early", "Middle", "Late")))
  
  #Makes 9 diagnostic plots of SALMOD output data 
  pdf(file = path4file)
  
  # Outliers
  # 1a Boxplot 
  
  df_sd <- aggregate(df$surv, by = list(df$period), sd)
  df_mn <- aggregate(df$surv, by = list(df$period), mean)
  boxplot <- ggplot(data = df, aes(x = period, y = surv)) +
    geom_jitter(colour = "slateblue", alpha = .8, size = 2.5) +
    geom_boxplot(outlier.colour = "red", alpha=.8, fill = "grey", outlier.size = 3, outlier.shape = 3 ) + 
    theme_classic() +
    ylab("Survival") +
    xlab("Time period\n(2010-2039, 2040-2069, 2070-2099)") +
    annotate("text", x = df_sd$Group.1, y =0.9, label = paste("sd =",round(df_sd$x, 2)), size = 6,
             colour = "red") +
    annotate("text", x = df_mn$Group.1, y =0.1, label = paste("mean =",round(df_mn$x, 2)), size = 6,
             colour = "red")
  
  print(boxplot)
  
  
  # 1b Clevelend dotplot 
  
  dotchart(df$surv)
  # conditioned on "period"  
  dotchart(df$surv, groups = df$period)
  
  # Homegeneity
  # 2. conditional boxplot
  
  # See above
  
  # Normality
  # 3a. Histogram geom_dotplot
  
  histogram <- ggplot(data = df, aes(x = surv)) +
    geom_histogram(aes(y = ..density..), fill = "grey30", binwidth = 0.1) +
    geom_density(size = 1.5, fill = "red", alpha = .5) +
    theme_few() +
    ylab("Frequency") +
    xlab("Survival")
  
  print(histogram)
  
  # 3b. QQ-plot
  
  qqnorm(vector)
  
  # Zero-trouble
  # 4. frequency plot
  
  dotplot <-  ggplot(data = df, aes(x = surv)) + 
    geom_dotplot(binwidth = 0.05) +
    theme_bw() +
    ylab("Frequency") +
    xlab("Survival")
  print(dotplot)    
  
  cond_dotplot <-  ggplot(data = df, aes(x = surv, fill = period)) + 
    geom_dotplot(binwidth = 0.05, stackgroups = TRUE, method = "histodot") +
    scale_fill_tableau("colorblind10") + 
    theme_bw() +
    ylab("Frequency") +
    xlab("Survival") 
  print(cond_dotplot)
  
  # Interactions
  
  coplots <- ggplot(data = df, aes(x = years, y = surv, colour = period)) +
    geom_line(size = 0.4, linetype = 2) +
    scale_colour_tableau("colorblind10") +
    geom_point(size = 3) +
    theme_few() +
    ylab("Survival") +
    xlab("Year")
  
  print(coplots)
  
  # Independence
  
  acf(df$surv)
  
  dev.off()  
} # end of Zuur_acm_univar()

salmod_surv_dp <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  mgmtScenarios <- c("BAU", "NoDiversion") 
  #"ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
  #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  stage <- c("prespawn", "egg", "fry")
  
  for (m in 1:length(mgmtScenarios))  {
    for (c in 1:length(climateScenarios)) {
      
      for (k in 1:length(gcms)) {
        salmod_out <- read.delim(file.path(sim_path, mgmtScenarios[m], paste0(climateScenarios[c],"_", gcms[k]), "SALMODsumOutMerge.txt"), sep = ",")
        mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
        rm(salmod_out)
        mortFW$gcms <- gcms[k]
        
        if(k == 1) mortFWgcms <- mortFW else mortFWgcms <- rbind(mortFWgcms, mortFW)
      } # end of k loop
      
      prespSurv <-  with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, PreSpawn = allMortSF/AFem))
      eggSurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Egg = FryGrad/Eggs))
      frySurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Fry = FryExit/FryGrad))
      prespSurvC <- dcast(prespSurv, year ~ gcms, value.var = "PreSpawn")
      eggSurv$Egg[which(is.na(eggSurv$Egg))] <- 0
      eggSurvC <- dcast(eggSurv, year ~ gcms, value.var = "Egg")
      frySurv$Fry[which(is.na(frySurv$Fry))] <- 0
      frySurvC <- dcast(frySurv, year ~ gcms, value.var = "Fry")
      
      for ( g in 2:7 ) {
        for ( s in 1:3 ) {    
          if (!dir.exists(file.path(".", "output", "diagnostic_plots", mgmtScenarios[m], climateScenarios[c]))) { 
            dir.create(file.path(".", "output", "diagnostic_plots", mgmtScenarios[m], climateScenarios[c]), recursive = TRUE) }
          dp_path <- file.path(".", "output", "diagnostic_plots", mgmtScenarios[m], climateScenarios[c])
          path4file <- file.path(dp_path, paste0(gcms[g-1], "_", stage[s], ".pdf")) 
          ifelse(stage[s] == "prespawn", vector <- prespSurvC[,g],
                 ifelse(stage[s] == "egg", vector <- eggSurvC[,g],
                        ifelse(stage[s] == "fry", vector <- frySurvC[,g])) )
          
          Zuur_acm_univar(vector, path4file)
          
        }
      }
    }
  }
  
} # end of salmod_surv_dp


wavelet_SALMOD_series <- function(management = c("BAU", "NoDiversion"),
                                  emissions = c("A2","B1"),
                                  global_mods = c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1") 
) {
  
  # Function that extracts a specific survival time series (Prespawn, egg and fry) 
  # from SALMOD output files for each climate and water management scenario.
  # Then it prints a pdf with plots of the ts, WPS, global WPS and scale-averaged WPS
  
  sim_path <- file.path(".", "data", "simulation_results")
  
  mgmtScenarios <- management #"ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
  #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- emissions
  gcms <- global_mods 
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1; I translated from your matlab code
  
  for (m in 1:length(mgmtScenarios))  {
    for (c in 1:length(climateScenarios)) {
      
      for (k in 1:length(gcms)) {
        salmod_out <- read.delim(file.path(sim_path, mgmtScenarios[m], 
                                           paste0(climateScenarios[c],"_", gcms[k]), "SALMODsumOutMerge.txt"), 
                                 sep = ",")
        mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
        mortFW$gcms <- gcms[k]
        
        if(k == 1) mortFWgcms <- mortFW else mortFWgcms <- rbind(mortFWgcms, mortFW)
      } # end of k loop
      
      prespSurv <-  with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, PreSpawn = allMortSF/AFem))
      eggSurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Egg = FryGrad/Eggs))
      frySurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Fry = FryExit/FryGrad))
      
      prespSurvC <- dcast(prespSurv, year ~ gcms, value.var = "PreSpawn")
      eggSurv$Egg[which(is.na(eggSurv$Egg))] <- 0
      eggSurvC <- dcast(eggSurv, year ~ gcms, value.var = "Egg")
      frySurv$Fry[which(is.na(frySurv$Fry))] <- 0
      frySurvC <- dcast(frySurv, year ~ gcms, value.var = "Fry")
      
      # For each cliamte scenario and water management alternative, plot the
      # the time series and wavelet spectra for each gcm for
      # prespawn survival, egg survival and fry survival
      # Calculate wavelet spectra
      # Visualize spectra
      dims <- dim(prespSurvC)
      out_wavelet <- file.path(".", "output", "wavelet_plots")
      if (!dir.exists(out_wavelet)) { dir.create(out_wavelet, recursive = TRUE) }
      
      path4file <- file.path(out_wavelet, paste0(mgmtScenarios[m], "_", climateScenarios[c], ".pdf"))
      pdf(file = path4file)
      
      for (i in 2:dims[2]) {
        #i <- 2
        # Time series
        par(fig= c(0,0.65,0.6,1))
        old <- par(mar = c(3, 3, 3, 1) )
        plot(2010:2099, prespSurvC[,i], type = "l", lwd = 2, col = "darkgrey",
             axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             las = 2, cex.axis = 1, tck = 0.02)
        box(lwd = 2)
        title(paste("GCM: ", gcms[i-1], " Prespawn survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))
        par(old)
        
        # Wavelet Power Spectrum
        par(fig= c(0,0.65,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(prespSurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
        plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        par(old)
        # Scale-averaged Power Spectrum - periods of greatest variability
        par(fig= c(0,0.65,0,0.3), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
        par(old)
        # Global Wavelet Power Spectrum
        par(fig= c(.65,1,0.3,.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        yrange <- NULL 
        y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
        FR <- log(rowSums(test_cwt$power.corr))
        plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
        axis(1)
        axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
        box()
        par(old)
      }
      for (i in 2:dims[2]) {
        #i <- 2
        par(fig= c(0,0.65,0.6,1))
        old <- par(mar = c(3, 3, 3, 1) )
        # plot egg surv ts
        plot(2010:2099, eggSurvC[,i], type = "l", lwd = 2, col = "darkgrey",
             axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             las = 2, cex.axis = 1, tck = 0.02)
        box(lwd = 2)
        title(paste("GCM: ", gcms[i-1], " Eggs survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))
        par(old)
        # plot wavelet power spectrum of egg survival
        par(fig= c(0,0.65,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(eggSurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
        plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        par(old)
        # Scale-averaged Power Spectrum - periods of greatest variability
        par(fig= c(0,0.65,0,0.3), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
        par(old)
        # Global Wavelet Power Spectrum
        par(fig= c(.65,1,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        yrange <- NULL 
        y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
        FR <- log(rowSums(test_cwt$power.corr))
        plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
        axis(1)
        axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
        box()
        par(old)
      }
      
      for (i in 2:dims[2]) {
        #i <- 2
        # plot fry survival time series 
        par(fig= c(0,0.65,0.6,1))
        old <- par(mar = c(3, 3, 3, 1) )
        plot(2010:2099, frySurvC[,i], type = "l", lwd = 2, col = "darkgrey",
             axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             las = 2, cex.axis = 1, tck = 0.02)
        box(lwd = 2)
        title(paste("GCM: ", gcms[i-1], " Fry survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))        
        par(old)
        # plot Wavelet power spectrum
        par(fig= c(0,0.65,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(frySurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
        plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        par(old)
        # Scale-averaged Power Spectrum - periods of greatest variability
        par(fig= c(0,0.65,0,0.3), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
        par(old)
        # Global Wavelet Power Spectrum
        par(fig= c(.65,1,0.3,.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        yrange <- NULL 
        y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
        FR <- log(rowSums(test_cwt$power.corr))
        plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
        axis(1)
        axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
        box()
        par(old)
      }
      
      
      dev.off()
      
    } # end of c climate emissions loop
  } # end of m water management loop
  
} # end of wavelet_SALMOD_series()

PDO_wvlt_test <- function() {
  # sim_path <- file.path(".", "data", "simulation_results")
  # setwd("~/NEP_salmon/chap_3/data/PDO/")
  # 
  # Testing the biwavelet function to show Loo 
  # Use PDO time series 
  pdo <- read.table(url("http://research.jisao.washington.edu/pdo/PDO.latest.txt"))
  # row.names = FALSE,
  # stringsAsFactors = FALSE,
  # sep = "\t",
  # skip = 34)
  pdo <- readLines(url("http://research.jisao.washington.edu/pdo/PDO.latest.txt"))
  lpdo <- length(pdo)
  pdo <- (pdo[-c(1:34, (lpdo-20):lpdo)])
  columns <- c("YEAR", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
  pdo_df <- data.frame(matrix(ncol = length(columns), nrow = 0))
  
  for (i in 1:length(pdo) ) {
    df <- read.table(textConnection(pdo[[i]]))
    pdo_df <- rbind(pdo_df, df)
  }
  colnames(pdo_df) <- columns
  
  # pdo <- read.table("PDO_2012SEP.txt", header = TRUE, sep = "\t", fill = TRUE)
  rowMeans(pdo[,2:13], na.rm = TRUE)
  # wavelet analysis of the PDO annual signal
  ann_pdo <- scale(rowMeans(pdo[,2:13], na.rm = TRUE)) # mean of 0, sd of 1, drop year 2012 - since incomplete
  years <- pdo[,1] # drop 2012 since incomplete
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1; I translated from your matlab code
  # Calculate the wavelet, using parameters as in your matlab code
  pdo.wt <- wt(cbind(years, ann_pdo), dj = 0.01, J1 = J1, max.scale = 32, 
               mother = "morlet", sig.test = 0, sig.level = 0.95) # acf(ann_pdo)$acf[2]
  
  # Compare to Flora's matlab: April 24
  # outputs from wavelet function call on PDO through 2012. Inputs in R and Matlab are the same
  # Wavelet coefficients: look the same
  # Period: same
  # Scale: same
  # Wavelet Power: different. matlab code rescales wavelet power spectrum to have units of variance.
  # R {biwavelet} also calculates bias corrected power
  # COI: same
  # signif: different. to get signif R divides wavelet power spectrum by product of var(t.s.) and signif (Pk), 
  # while matlab only divides by signif.
  # Plot time series and the wavelet power spectrum
  old <- par(mfrow = c(2,1), mar = c(3,4,0.5,0.5))
  plot(years, ann_pdo, type = "l", lwd = 3.7, lty = 1, col = "black", xlab = "Year", ylab = "PDO")
  lines(years, ann_pdo, type = "l", lwd = 2, col = "white")
  lines(years, ann_pdo, lwd = 2, lty = 2)
  box(lwd = 2)
  plot(pdo.wt, type = "power.corr.norm", tol = 0.95)
  box(lwd = 2)
  par(old)
  
} # end of PDO_wvlt_test()

# # Plot global wavelt power spectrum
# 
# wt_plot_GWPS <- function(biwv_obj, labs = TRUE, small = TRUE) {
#   library(biwavelet)
#   if (small == TRUE) vold <- par(mar = c(1,2,1,1), cex = .7) else vold <- par(mar = c(4,4,0.5,0.5))
#   yrange <- NULL 
#   y_ticks <- 2^(floor(log2(min(biwv_obj$period, yrange))):(floor(log2(max(biwv_obj$period, yrange)))+1))
#   FR <- log(rowMeans(biwv_obj$power.corr))
#   if (labs == TRUE) {
#     xlabel <- "Relative Magnitude"
#     ylabel <- "Period"
#   } else {
#     xlabel <- ""
#     ylabel <- ""
#   }
#   plot(FR, log2(biwv_obj$period), type = "n", ylim = rev(range(log2(biwv_obj$period))), 
#        axes = FALSE, ylab = ylabel, xlab = ylabel)
#   lines(FR, log2(biwv_obj$period), lty = 1, lwd = 4)
#   lines(FR, log2(biwv_obj$period), lty = 1, lwd = 3, col = "white")
#   lines(FR, log2(biwv_obj$period), lty = 2, lwd = 2.5)
#   axis(1)
#   axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
#   box(lwd = 2)
#   par(vold)
# }

plot_gen_freq_wvlt <- function(noise = noiseList, 
                               burn_in_pd,
                               num_rows2plt = 100,
                               n = 1, 
                               J1 = trunc((log(32/(2 * 1))/log(2))/0.01)) {
  
  line_d <- 2
  ylim <- c(0,2)
  lwd_ts <- 1.5
  rows2plot <- (burn_in_pd+1):(burn_in_pd+num_rows2plt)
  old <- par(mar = c(3,2,2,1), cex = .7)
  # white noise
  # Label
  print("white")
  par(fig=c(0, 0.2, 0.80, 1))
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "White noise", cex = 0.8)
  
  # Generating spectrum
  par(fig=c(0.2, 0.55, 0.80, 1), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 1), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i", yaxt='n')
  axis(2, at = c(0,2), labels = c(0,1))
  rect(xleft = 0, xright  = 0.5, ybottom = 0,  ytop = 0.2, col="gray")
  mtext("a", side = 2, las = 1, at = 2, line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  par(fig=c(0.55, 0.9, 0.80, 1), new = TRUE)
  white.wt <- wt(cbind(1:num_rows2plt, noiseList[[1]][rows2plot,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  mtext("b", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  # Bandpass period 3-4
  # label
  print("cohort")
  par(fig=c(0, 0.20, 0.60, 0.80), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "Period 3-4", cex = 0.8)
  
  # Generating spectrum
  par(fig=c(0.20, 0.55, 0.60, 0.80), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i", yaxt='n')
  axis(2, at = c(0,2), labels = c(0,1))
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  mtext("c", side = 2, las = 1, at = 2, line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  par(fig=c(0.55, 0.9, 0.60, 0.80), new = TRUE)
  p34.wt <- wt(cbind(1:num_rows2plt, noiseList[[2]][rows2plot,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34.wt)
  mtext("d", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  # Bandpass greater than period 10
  # label
  print("low")
  par(fig=c(0, 0.2, 0.40, 0.60), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "Low Frequency", cex = 0.8)
  
  # Generating spectrum
  par(fig=c(0.2, 0.55, 0.40, 0.60), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i", yaxt='n')
  axis(2, at = c(0,2), labels = c(0,1))
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
  mtext("e", side = 2, las = 1, at = 2, line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  par(fig=c(0.55, 0.9, 0.40, 0.60), new = TRUE)
  pgt10.wt <- wt(cbind(1:num_rows2plt, noiseList[[3]][rows2plot,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt10.wt)
  mtext("f", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  # Bandpass greater than period 10 and period 3-4
  # label
  print("both")
  par(fig=c(0, 0.2, 0.2, 0.4), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "Both Cohort and\nLow Frequency", cex = 0.8)
  
  # Generating spectrum
  par(fig=c(0.2, 0.55, 0.2, 0.4), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i", yaxt='n')
  axis(2, at = c(0,2), labels = c(0,1))
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 0.5, col="gray")
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 0.5, col="gray")
  mtext("g", side = 2, las = 1, at = 2, line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  par(fig=c(0.55, 0.9, 0.2, 0.4), new = TRUE)
  p34gt10.wt <- wt(cbind(1:num_rows2plt, noiseList[[4]][rows2plot,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34gt10.wt)
  mtext("h", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  par(old)
  # one over f noise (beta = 1)
  # label
  old <- par(mar = c(4,2,2,1), cex = .7)
  print("1/f^b")
  par(fig=c(0, 0.2, 0, 0.2), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, expression(frac(1, sqrt(f))), cex = 1)
  # text(5,5, expression(1/f^b ~ 'b = 0.5'), cex = 0.8)
  
  # Generating spectrum
  par(fig=c(0.2, 0.55, 0, 0.2), new = TRUE)
  Ns <- 50
  xs <- seq(0, 0.5, length = Ns)
  # ys <- mk_1_over_f_beta(N = Ns*2, beta = 1)
  ys <- mk_1_over_f_beta(N = Ns, beta = 0.5, freq = Ns)
  scaled_ys <- 2 * (ys / max(ys[-1]))
  int_xy <- approx(xs[-1], scaled_ys[-1], xout = seq(0, 0.5, by = 0.005), rule = 2:2)
  plot(x = int_xy$x, y = int_xy$y, type = "l",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i", yaxt='n')
  axis(2, at = c(0,2), labels = c(0,1))
  polygon(c(int_xy$x,rev(int_xy$x)), 
          c(rep(0, length = length(int_xy$x)), rev(int_xy$y)),
          col="gray")
  mtext("i", side = 2, las = 1, at = 2, line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  par(fig=c(0.55, 0.9, 0, 0.2), new = TRUE)
  p_one_over_f.wt <- wt(cbind(1:num_rows2plt, noiseList[[5]][rows2plot,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p_one_over_f.wt, xlab = "Year")
  mtext("j", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  par(fig=c(0.4, 1, 0, 1), new = TRUE)
  image.plot(legend.only = TRUE, zlim = c(0,64), add = FALSE) 
  par(old)
  
} # end plot_gen_freq_wvlt

# grab_zlims <- function(x) {
#   stopifnot(class(x) == "biwavelet")
#   x$power <- x$power.corr
#   zvals <- log2(abs(x$power/x$sigma2))
#   zlim <- range(c(-1, 1) * max(zvals))
#   zvals[zvals < zlim[1]] <- zlim[1]
#   return(zvals)
# }


plot_surv_spawn_ts <- function(spawners = storage, 
                               noise = noiseList,
                               burn_in_pd,
                               phasein_len,
                               num_rows2plt = 100,
                               sim_len,
                               meanSurv = "0.5",
                               sigPSmult = "0.2", 
                               n = 1, 
                               J1 = trunc((log(32/(2 * 1))/log(2))/0.01)) {
  
  spPlot <- copy(spawners[ i = N > burn_in_pd & reps_c == n & sigPSmult_c == sigPSmult & meanPS_c == meanSurv ])
  
  out <- vector("list", length = length(noise))
  for (i in 1:length(noise)) {
    out[[i]] <- make_surv_mat(noise_dat = noise[[i]], 
                              mean_surv = as.numeric(meanSurv), 
                              sd_surv = as.numeric(sigPSmult),
                              sim_len = sim_len,
                              burn_in = burn_in,
                              phasein_len = phasein_len) 
  }
  
  line_d <- 2.2
  sp_ylim <- c(0, 6000)
  s_yim <- c(0, 1)
  lwd_ts <- 1.5
  n_rows <- 1:num_rows2plt
  rows2plot <- (burn_in_pd+1):(burn_in_pd+num_rows2plt)
  qeT <- 100
  
  old <- par(mar = c(3,2,1,2), cex = .7)
  
  # white noise
  # Label
  print("white")
  par(fig=c(0, 0.20, 0.80, 1))
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "White noise", cex = 1)
  
  # time series plot: survival 
  par(fig=c(0.20, 0.60, 0.80, 1), new = TRUE)
  plot(n_rows, out[[1]][rows2plot, n], type = "l", col = "grey20", lwd = lwd_ts, 
       ylim = s_yim, xlab = "", ylab = "")
  mtext("a", side = 2, las = 1, at = 1, line = line_d, cex = 1.2)
  
  # time series plot: spawners 
  par(fig=c(0.60, 1, 0.80, 1), new = TRUE)
  plot(n_rows, spPlot[i = n_rows, j = white], type = "l", col = "grey20", 
       lwd = lwd_ts, xlab = "", yaxt = "n", ylim = sp_ylim)
  abline(h = qeT, col = "black", lty = 2)
  axis(side = 2, at = seq(0, sp_ylim[2], by = sp_ylim[2]/2), 
       labels = seq(0, sp_ylim[2], by = sp_ylim[2]/2))
  mtext("b", side = 2, las = 1, at = sp_ylim[2], line = line_d, cex = 1.2)
  # Bandpass period 3-4
  # label
  print("cohort")
  par(fig=c(0, 0.20, 0.6, 0.8), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "Period 3-4", cex = 1)
  
  # time series plot: survival 
  par(fig=c(0.20, 0.60, 0.6, 0.8), new = TRUE)
  plot(n_rows, out[[2]][rows2plot, n], type = "l", col = "grey20", lwd = lwd_ts, 
       ylim = s_yim, xlab = "", ylab = "")
  mtext("c", side = 2, las = 1, at = 1, line = line_d, cex = 1.2)
  
  # time series plot: spawners
  par(fig=c(0.6, 1, 0.6, 0.8), new = TRUE)
  plot(n_rows, spPlot[i = n_rows, j = p34], type = "l", col = "grey20", 
       lwd = lwd_ts, xlab = "", yaxt="n", ylim = sp_ylim)
  abline(h = qeT, col = "black", lty = 2)
  axis(side = 2, at = seq(0, sp_ylim[2], by = sp_ylim[2]/2), 
       labels = seq(0, sp_ylim[2], by = sp_ylim[2]/2))
  mtext("d", side = 2, las = 1, at = sp_ylim[2], line = line_d, cex = 1.2)
  
  # Bandpass greater than period 10
  # label
  print("low")
  par(fig=c(0, 0.2, 0.4, 0.6), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "Low Frequency", cex = 1)
  
  # time series plot: survival 
  par(fig=c(0.20, 0.60, 0.4, 0.6), new = TRUE)
  plot(n_rows, out[[3]][rows2plot, n], type = "l", col = "grey20", lwd = lwd_ts, 
       ylim = s_yim, xlab = "", ylab = "")
  mtext("e", side = 2, las = 1, at = 1, line = line_d, cex = 1.2)
  
  # time series plot: spawners
  par(fig=c(0.6, 1, 0.4, 0.6), new = TRUE)
  plot(n_rows, spPlot[i = n_rows, j = pgt10], type = "l", col = "grey20", 
       lwd = lwd_ts, xlab = "", yaxt="n", ylim = sp_ylim)
  abline(h = qeT, col = "black", lty = 2)
  axis(side = 2, at = seq(0, sp_ylim[2], by = sp_ylim[2]/2), 
       labels = seq(0, sp_ylim[2], by = sp_ylim[2]/2))
  mtext("f", side = 2, las = 1, at = sp_ylim[2], line = line_d, cex = 1.2)
  
  # Bandpass greater than period 10 and period 3-4
  # label
  print("both")
  par(fig=c(0, 0.2, 0.2, 0.4), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, "Both Cohort and\nLow Frequency", cex = 1)
  
  # time series plot: survival 
  par(fig=c(0.20, 0.60, 0.2, 0.4), new = TRUE)
  plot(n_rows, out[[4]][rows2plot, n], type = "l", col = "grey20", lwd = lwd_ts, 
       ylim = s_yim, xlab = "", ylab = "")
  mtext("g", side = 2, las = 1, at = 1, line = line_d, cex = 1.2)
  
  # time series plot: spawners
  par(fig=c(0.6, 1, 0.2, 0.4), new = TRUE)
  plot(n_rows, spPlot[i = n_rows, j = p34gt10], type = "l", col = "grey20", 
       lwd = lwd_ts, xlab = "", yaxt="n", ylim = sp_ylim)
  abline(h = qeT, col = "black", lty = 2)
  axis(side = 2, at = seq(0, sp_ylim[2], by = sp_ylim[2]/2), 
       labels = seq(0, sp_ylim[2], by = sp_ylim[2]/2))
  mtext("h", side = 2, las = 1, at = sp_ylim[2], line = line_d, cex = 1.2)
  par(old)
  # one over f noise (beta = 1)
  # label
  old <- par(mar = c(4,2,1,2), cex = .7)
  print("1/f^b")
  par(fig=c(0, 0.2, 0, 0.2), new = TRUE)
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  # text(5,5, expression(1/f^b ~ 'b = 0.5'), cex = 1)
  text(5,5, expression(frac(1, sqrt(f))), cex = 1)
  
  # time series plot: survival 
  par(fig=c(0.20, 0.60, 0, 0.2), new = TRUE)
  plot(n_rows, out[[5]][rows2plot, n], type = "l", col = "grey20", lwd = lwd_ts, 
       ylim = s_yim, xlab = "Years", ylab = "")
  mtext("i", side = 2, las = 1, at = 1, line = line_d, cex = 1.2)
  
  # time series plot: spawners
  par(fig=c(0.6, 1, 0, 0.2), new = TRUE)
  plot(n_rows, spPlot[i = n_rows, j = one_over_f], type = "l", col = "grey20", 
       lwd = lwd_ts, xlab = "Years", yaxt="n", ylim = sp_ylim)
  abline(h = qeT, col = "black", lty = 2)
  axis(side = 2, at = seq(0, sp_ylim[2], by = sp_ylim[2]/2), 
       labels = seq(0, sp_ylim[2], by = sp_ylim[2]/2))
  mtext("j", side = 2, las = 1, at = sp_ylim[2], line = line_d, cex = 1.2)
  
  par(old)
  
} # end plot_surv_spawn_ts()



# summaryCh3tsPlotNoise <- function(noise = noiseList, n = 1, J1 = trunc((log(32/(2 * 1))/log(2))/0.01)) {
#   
#   ylim <- c(0,2)
#   lwd_ts <- 1.5
#   n_rows <- 1:200
#   
#   old <- par(mar = c(1,2,1,1), cex = .7)
#   # white noise
#   # Label
#   print("white")
#   par(fig=c(0, 0.16, 0.85, 1))
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "White noise", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.85, 1), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.85, 1), new = TRUE)
#   plot(n_rows, noiseList[[1]][n_rows,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.85, 1), new = TRUE)
#   white.wt <- wt(cbind(n_rows, noiseList[[1]][n_rows,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(white.wt)
#   
#   # Bandpass period 3-4
#   # label
#   print("cohort")
#   par(fig=c(0, 0.16, 0.7, 0.85), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "Period 3-4", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.7, 0.85), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.7, 0.85), new = TRUE)
#   plot(n_rows, noiseList[[2]][n_rows,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.7, 0.85), new = TRUE)
#   p34.wt <- wt(cbind(n_rows, noiseList[[2]][n_rows,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p34.wt)
#   
#   # Bandpass greater than period 10
#   # label
#   print("low")
#   par(fig=c(0, 0.16, 0.55, 0.7), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "Low Frequency", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.55, 0.7), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.55, 0.7), new = TRUE)
#   plot(n_rows, noiseList[[3]][n_rows,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.55, 0.7), new = TRUE)
#   pgt10.wt <- wt(cbind(n_rows, noiseList[[3]][n_rows,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(pgt10.wt)
#   
#   
#   # Bandpass greater than period 10 and period 3-4
#   # label
#   print("both")
#   par(fig=c(0, 0.16, 0.4, 0.55), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "Both Cohort and\nLow Frequency", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.4, 0.55), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
#   rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.4, 0.55), new = TRUE)
#   plot(n_rows, noiseList[[4]][n_rows,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.4, 0.55), new = TRUE)
#   p34gt10.wt <- wt(cbind(n_rows, noiseList[[4]][n_rows,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p34gt10.wt)
#   
#   # one over f noise (beta = 1)
#   # label
#   print("1/f^b")
#   par(fig=c(0, 0.16, 0.25, 0.4), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, expression(1/f^b ~ 'b = 0.5'), cex = 1)
# 
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.25, 0.4), new = TRUE)
#   Ns <- 50
#   xs <- seq(0,0.5, length = Ns)
#   ys <- mk_1_over_f_beta(N = Ns*2, beta = 1)
#   scaled_ys <- 2 * (ys / max(ys[-1]))
#   plot(x = xs[-1], y = scaled_ys[-1], type = "l",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   polygon(c(xs[-1],rev(xs[-1])),
#           c(rep(0, length = Ns-1), rev(scaled_ys[-1])),
#           col="gray")
# 
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.25, 0.4), new = TRUE)
#   plot(n_rows, noiseList[[5]][n_rows,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.25, 0.4), new = TRUE)
#   p_one_over_f.wt <- wt(cbind(n_rows, noiseList[[5]][n_rows,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p_one_over_f.wt)
#   
#   # ar-noise phi = 0.5   
#   # label
#   print("AR1")
#   par(fig=c(0, 0.16, 0.1, 0.25), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, expression(AR1 ~ phi: 0.5), cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.1, 0.25), new = TRUE)
#   plot(x = 0:2, y = 0:2, type = "n",
#        axes = F, xlab = "", ylab = "")
#   text(1,1, "Not applicable", cex = 1)
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.1, 0.25), new = TRUE)
#   plot(n_rows, noiseList[[7]][n_rows,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.1, 0.25), new = TRUE)
#   p_ar1.wt <- wt(cbind(n_rows, noiseList[[7]][n_rows,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p_ar1.wt)
#   
#   par(old)
#   
# } # end summaryCh3tsPlotNoise


# summaryCh3tsPlotSpawners <- function(spawners = storage, 
#                                      meanSurv = 0.5,
#                                      sigma = 0.2, 
#                                      n = 1, 
#                                      J1 = trunc((log(32/(2 * 1))/log(2))/0.01)) {
#   
#   spPlot <- copy(spawners[ i = N > 400 & reps_c == 1 & sigPSmult_c == sigma & meanPS_c == meanSurv ])
#   
#   ylim <- c(0,2)
#   lwd_ts <- 1.5
#   n_rows <- 1:200
#   
#   old <- par(mar = c(1,2,1,1), cex = .7)
#   
#   # white noise
#   # Label
#   print("white")
#   par(fig=c(0, 0.16, 0.85, 1))
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "White noise", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.85, 1), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.85, 1), new = TRUE)
#   plot(n_rows, spPlot[i = n_rows, j = white], type = "l", col = "slateblue", 
#        lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.85, 1), new = TRUE)
#   white.wt <- wt(cbind(n_rows, spPlot[i = n_rows, j = white]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(white.wt)
#   
#   # Bandpass period 3-4
#   # label
#   print("cohort")
#   par(fig=c(0, 0.16, 0.7, 0.85), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "Period 3-4", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.7, 0.85), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.7, 0.85), new = TRUE)
#   plot(n_rows, spPlot[i = n_rows, j = p34], type = "l", col = "slateblue", 
#        lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.7, 0.85), new = TRUE)
#   p34.wt <- wt(cbind(n_rows, spPlot[i = n_rows, j = p34]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p34.wt)
# 
#   # Bandpass greater than period 10
#   # label
#   print("low")
#   par(fig=c(0, 0.16, 0.55, 0.7), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "Low Frequency", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.55, 0.7), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.55, 0.7), new = TRUE)
#   plot(n_rows, spPlot[i = n_rows, j = pgt10], type = "l", col = "slateblue", 
#        lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.55, 0.7), new = TRUE)
#   pgt10.wt <- wt(cbind(n_rows, spPlot[i = n_rows, j = pgt10]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(pgt10.wt)
#   
#   ########### REF
#   # Bandpass greater than period 10 and period 3-4
#   # label
#   print("both")
#   par(fig=c(0, 0.16, 0.4, 0.55), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, "Both Cohort and\nLow Frequency", cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.4, 0.55), new = TRUE)
#   plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
#   rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.4, 0.55), new = TRUE)
#   plot(n_rows, spPlot[i = n_rows, j = p34gt10], type = "l", col = "slateblue", 
#        lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.4, 0.55), new = TRUE)
#   p34gt10.wt <- wt(cbind(n_rows, spPlot[i = n_rows, j = p34gt10]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p34gt10.wt)
#   #### DOWN
#   # one over f noise (beta = 1)
#   # label
#   print("1/f^b")
#   par(fig=c(0, 0.16, 0.25, 0.4), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, expression(1/f^b ~ 'b = 0.5'), cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.25, 0.4), new = TRUE)
#   Ns <- 50
#   xs <- seq(0,0.5, length = Ns)
#   ys <- mk_1_over_f_beta(N = Ns*2, beta = 1)
#   scaled_ys <- 2 * (ys / max(ys[-1]))
#   plot(x = xs[-1], y = scaled_ys[-1], type = "l",
#        xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
#        xaxs = "i", yaxs = "i")
#   polygon(c(xs[-1],rev(xs[-1])),
#           c(rep(0, length = Ns-1), rev(scaled_ys[-1])),
#           col="gray")
#   
#   # time series plot
#   par(fig=c(0.44, 0.72, 0.25, 0.4), new = TRUE)
#   plot(n_rows, spPlot[i = n_rows, j = one_over_f], type = "l", col = "slateblue", 
#        lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.25, 0.4), new = TRUE)
#   p_one_over_f.wt <- wt(cbind(n_rows, spPlot[i = n_rows, j = one_over_f]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p_one_over_f.wt)
#   
#   # ar-noise phi = 0.5   
#   # label
#   print("AR1")
#   par(fig=c(0, 0.16, 0.1, 0.25), new = TRUE)
#   plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
#   text(5,5, expression(AR1 ~ phi: 0.5), cex = 1)
#   
#   # Generating spectrum
#   par(fig=c(0.16, 0.44, 0.1, 0.25), new = TRUE)
#   plot(x = 0:2, y = 0:2, type = "n",
#        axes = F, xlab = "", ylab = "")
#   text(1,1, "Not applicable", cex = 1)
#   
#   # time series plot:
#   par(fig=c(0.44, 0.72, 0.1, 0.25), new = TRUE)
#   plot(n_rows, spPlot[i = n_rows, j = ar1], type = "l", col = "slateblue", 
#        lwd = lwd_ts, xlab = "Time (years)", ylab = "")
#   
#   # plot wavelet power spectrum
#   par(fig=c(0.72, 1, 0.1, 0.25), new = TRUE)
#   p_ar1.wt <- wt(cbind(n_rows, spPlot[i = n_rows, j = ar1]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
#   plot(p_ar1.wt)
#   
#   par(old)
#   
# } # end summaryCh3tsPlotSpawners
# 
# corrCoefDist <- function(white, wave34, rsw34, hi  = 1, lo = -(hi)) {
#   
#   w <- cor(white)
#   wtCorrs <- w[col(w) < row(w)]
#   print(range(wtCorrs))
#   wt <- cor(wave34)
#   wv34Corrs <- wt[col(wt) < row(wt)]
#   rsw <- cor(rsw34)
#   rsw34Corrs <- rsw[col(rsw) < row(rsw)]
#   
#   brk <- seq(lo, hi, by = 0.05)
#   
#   empCIwt  <- quantile(wtCorrs, c(0.05, 0.95))
#   empCIwv34  <- quantile(wv34Corrs, c(0.05, 0.95))
#   empCIrsw34  <- quantile(rsw34Corrs, c(0.05, 0.95))
#   
#   hist(wtCorrs, breaks = brk, col = rgb(0,0,0,.8), xlim = c(lo, hi), main = "", 
#        xlab = expression(paste(rho)))
#   hist(wv34Corrs, breaks = brk, col = rgb(1,0,0,.6), add = TRUE)
#   hist(rsw34Corrs, breaks = brk, col = rgb(0,0,1,.4), add = TRUE)
#   
#   abline(v = empCIwt[1], col = rgb(0,0,0,1), lty = 2, lwd = 2)
#   abline(v = empCIwt[2], col = rgb(0,0,0,1), lty = 2, lwd = 2)
#   abline(v = empCIwv34[1], col = rgb(1,0,0,1), lty = 2, lwd = 2)
#   abline(v = empCIwv34[2], col = rgb(1,0,0,1), lty = 2, lwd = 2)
#   abline(v = empCIrsw34[1], col = rgb(0,0,1,1), lty = 2, lwd = 2)
#   abline(v = empCIrsw34[2], col = rgb(0,0,1,1), lty = 2, lwd = 2)
#   legend("topright", c("white", "wavelet34", "rsw34"), col = c("black", "red", "blue"), lty = 1, lwd = 3, cex = .8)
# }
# 

# linesMeanFR_DTmanyAR <- function(dataTable, N, surv, scale = "CV", line_color, AR_col) {
#   # spectral frequency 
#   if (trunc(sqrt(N)) %% 2 == 0) {
#     m <- trunc(sqrt(N)) + 1 
#   } else {
#     m <- trunc(sqrt(N))
#   }
#   
#   spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
#   
#   freq <- (1:(N/2))/N
#   for (i in 1:ncol(spcMean)) {
#     
#     ts <- as.ts(droplevels(subset(dataTable, reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001, select = AR_col) ))
#     
#     ifelse(all(ts == 0), 
#            spcMean[,i] <- rep(NA, times = N/2),
#            ifelse(scale == "CV", 
#                   spcMean[,i] <- periodogram(ts/mean(ts), plot = F )$spec, 
#                   spcMean[,i] <- periodogram(scale(ts), plot = F )$spec ) )
#   }
#   mean_spc <- rowMeans(spcMean, na.rm = TRUE)
#   
#   spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
#   spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
#   spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
#   spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
#   
#   lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = line_color)
#   
#   box(lwd=2)
# }
