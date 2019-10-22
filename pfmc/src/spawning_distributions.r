# Spawning Biomass Distributions

# Plot spawning biomass distribution -- treat as probability distribution
# y axis = probability of spawning
# x axis = age

# run life table information and load functions
source("C:/Users/Mikaela/Documents/GitHub/pfmc/src/life_table_calculations.r")

# load parms df
parms = read.csv("C:/Users/Mikaela/Documents/GitHub/pfmc/parms_data/pfmc_parms.csv",
                 header=TRUE,stringsAsFactors = FALSE)
parms = as.data.frame(parms)
# remove problem species:
#parms <- parms[parms$problem_spp =="no",]

# store distribution CV, mode, sd here:
cvs_mode = rep(NA, length=length(parms$spp))
mode_age = rep(NA, length=length(parms$spp))
sd_mode = rep(NA, length=length(parms$spp))

p <- as.list(parms$spp)
names(p) <- parms$spp

spawning_dist_plotting_dataL <- as.list(parms$spp)
names(spawning_dist_plotting_dataL) <- parms$spp

for (i in 1:length(parms$spp)) { # step through each spp
  
  if(parms[parms$spp == parms$spp[i],]$problem_spp == "yes") {#if a problem, then
    p[[i]] <- parms$spp[i]} 
  else {
  
  # define parms for LEP calculations:
  maxage <- as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  L_a    <- La_list[[i]]
  B_a    <-  Ba_list[[i]]
  propmat_a <-  prop_list[[i]]
  eggs_a <- eggs_list[[i]]
  vul_a  <- vul_list[[i]]
  M      <- as.numeric(parms[parms$spp == parms$spp[i],]$M)
  
  # calculate LEP at each age
  out <- assemble_Leslie(maxage=maxage,
                                  L_a=L_a,
                                  B_a=B_a,
                                  propmat_a=propmat_a,
                                  eggs_a=eggs_a,
                                  vul_a=vul_a,
                                  M=M,
                                  Fval=0) #for now, set F=0
  
  # calculate probability of spawning at age (LEP-at-age/sumLEP)
  Age <- seq(from=1,to=maxage,by=1)
  p_spawn = as.data.frame(out$LTABLE$LEP_a/sum(out$LTABLE$LEP_a))
  keep= cbind(p_spawn,Age)
  colnames(keep) <- c("p_spawn","Age")
  keep$spp <- rep(parms$spp[i],length=length(keep[,1]))
  
  # using mode in sd
  mode_age[i] = keep$Age[which.max(keep$p_spawn)] # what is the age with highest probability?
  sd_mode[i] = round(sqrt(sum(keep$p_spawn*(keep$Age-mode_age[i])^2) ),digits=2) # stdev
  cvs_mode[i] = round(sd_mode[i]/mode_age[i],digits=2) # coefficient of variation 
  
  # Plot spawning distribution for each population:
  p[[i]] <- ggplot(keep,aes(x=Age,y=p_spawn)) +
    geom_line() + theme_classic() + xlab("") + 
    ylab("") + 
    ggtitle(paste(parms$spp[i]," sd=",sd_mode[i]," cv=",cvs_mode[i],sep="")) +
    scale_y_continuous(limits = c(0,0.35)) + #y axis for not adjusted
    xlim(0,110) +
    geom_vline(xintercept=mode_age[i],linetype="dashed") +
    geom_text(x=(mode_age[i]+2), y=0.3, label=mode_age[i], size=4) +
    theme(text = element_text(size = 10))
  
  # Save keep df to plot (for plotting later)
  spawning_dist_plotting_dataL[[i]] <- keep
  }

}
spp <- parms$spp
problem_spp <- parms$problem_spp
maxage <- parms$maxage
spawndistmetrics <- as.data.frame(cbind(spp,problem_spp,mode_age,sd_mode,cvs_mode,maxage))
spawning_dist_data <- bind_rows(spawning_dist_plotting_dataL[parms[parms$problem_spp == "no",]$spp])
rm(spp,problem_spp,i,Age,p_spawn,keep,out,maxage,L_a,B_a,propmat_a,vul_a,M,spawning_dist_plotting_dataL)

# only run code below if I want to generate new plots:
#subset p to only spp that are not problems
#rm(p)
psub <- p[parms[parms$problem_spp == "no",]$spp]

#tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/spawning_distribution_v2.tiff', units="in", width=12, height=9, res=300)
#do.call(grid.arrange,c(psub,ncol=4,top="LEP-at-age/total LEP, where LEP-at-age = eggs-at-age*p(mature-at-age)*survival-to-that-age",left="Pr(spawning)"))
#dev.off()




