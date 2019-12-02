# *** no longer abundance ts
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

plot_wvlt_CVnotstd_AUCnot1 <- function(burn_in_pd = rm_first_timesteps, #remove from beginning of ts
                                       num_rows2plt = 1000, #number of years to plot
                                       timesteps = timesteps,
                                       span.multiplier = 1,
                                       noise = noisetypes[[n]],
                                       noisename = noisenames[n],
                                       ts.data = ts.data, #timeseries dataset
                                       ts.spec = ts.spec, #spectrum dataset
                                       spawning_dist_data=spawning_dist_data) {
  
  # ---
  rows2plot <- (burn_in_pd+1):(burn_in_pd+num_rows2plt)
  n = 1 #for indexing column n in noise df?
  J1 = trunc((log(32/(2 * 1))/log(2))/0.01)
  test = subset(parms[parms$problem_spp=="no",],select=c("spp","maxage"))
  test$maxage <- as.numeric(test$maxage)
  test <- test[order(test$maxage),]
  spp <- test$spp #put spp in order of max age
  #df.list <- df.list[spp] #set order
  line_d <- 2 #which margin to start counting at zero, 2=left
  tsylim <- c(min(ts.data[ts.data$year %in% rows2plot &
                            ts.data$EI==0 &
                            ts.data$variable=="Nsize",]$value ),
              max(ts.data[ts.data$year %in% rows2plot &
                            ts.data$EI==0 &
                            ts.data$variable=="Nsize",]$value))
  
  length(ts.data[ts.data$year %in% rows2plot &
                   ts.data$EI==0,]$value)
  spwn_ylim <- c(0,0.35) #ylimits for spawning biomass dist
  spwn_xlim <- c(0,110)  #xlimits for spawning biomass dist
  
  n_rows <- seq(from=1,to=length(rows2plot))
  old <- par(mar = c(3,2,2,1), cex = .7)
  specminval <- as.numeric(min(ts.spec$spec))
  specmaxval <- as.numeric(max(ts.spec$spec))
  # # set Daniell smoother for frequency response plot
  # tmp <- ceiling(sqrt(length(1:(timesteps-burn_in_pd-1)))) #sq root of timeseries lgth, rounded
  # if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp} #make it odd, if the square root is even
  # m = m * span.multiplier
  
  # mtext()
  #side 2 is left
  #las 1:horizontal 2:perpendicular to axis 3:vertical
  #at position along y axis
  #line on which MARgin line, starting at 0 counting outwards
  # ---
  
  tiff(file=paste('C:/Users/Mikaela/Documents/GitHub/pfmc/results/wavelet_plots/3_A_CVstd_AUC1_',noisename,'.tiff',sep=""), 
       units="in", width=8, height=11, res=300) 
  
  par(mfrow=c(7,4),mai=c(0.2,0.3,0.2,0.2))
  
  # --------------
  # 1) White noise
  # Label 
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  text(5,5, noisename, cex = 1.3)
  # Plot recruit time series
  plot(x = n_rows,
       y = noise[rows2plot], cex=0.8, ylab="recruit",
       type = "l", xaxs = "i", yaxs = "i")
  # plot wavelet power spectrum
  noise.wt <- wt(cbind(1:num_rows2plt, noise[rows2plot]), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(noise.wt)
  # Plot frequency response
  p_spec <- spec.pgram(x=noise,
                       spans=c(m,m),taper=0.1,plot = FALSE,demean = TRUE)
  plot(x=p_spec$freq,y=p_spec$spec,xaxs = "i", yaxs = "i",
       type = "l",log="y",yaxt="n",ylab="")
  rm(p_spec,noise.wt)
  
  # --------------
  i=1
  spp[i]
  
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  # --------------
  i=2
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  i=3
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  i=4
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  i=5
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  i=6
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  dev.off()
  
  
  # --------
  tiff(file=paste('C:/Users/Mikaela/Documents/GitHub/pfmc/results/wavelet_plots/3_B_CVstd_AUC1_',noisename,'.tiff',sep=""), 
       units="in", width=8, height=11, res=300) 
  par(mfrow=c(7,4),mai=c(0.2,0.3,0.2,0.2))
  
  i=7
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  i=8
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  
  i=9
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  
  i=10
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  
  i=11
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  
  
  i=12
  spp[i]
  # Label 
  #par(fig=c(0, 0.15, 0.80, 1))
  #plot(1:10, type = "n", axes = F, xlab = "", ylab = "")
  #text(5,5, spp[i], cex = 1.3)
  #text(5,5, "max age = ", cex = 1.3)
  # Plot spawning distribution
  plot(x = spawning_dist_data[spawning_dist_data$spp == spp[i],]$Age,
       y = spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn,
       type = "l", ylim=spwn_ylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="Age",xlab="p(spawning)",xlim=spwn_xlim)
  title(ylab="pr(spawn)",line=0)
  text(x=100,y=0.2,label=paste(spp[i],"\nmax age=",parms[parms$spp==spp[i],]$maxage,
                               "\nCV=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$cvs_mode,
                               "\npeak=",spawndistmetrics[spawndistmetrics$spp==spp[i],]$mode_age,sep=""),
       cex=1.2,pos=2)
  segments(x0=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           x1=which.max(spawning_dist_data[spawning_dist_data$spp == spp[i],]$p_spawn),
           y0=0,y1=0.33,lty=2,col="red",lwd=2)
  # Plot recruit time series
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                     ts.data$EI==0 &
                     ts.data$noise==noisename &
                     ts.data$variable=="Nsize" &
                     ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(SSB)",line=0)
  # plot wavelet power spectrum
  ts <- ts.data[ts.data$spp==spp[i] & 
                  ts.data$EI==0 & 
                  ts.data$noise==noisename &
                  ts.data$variable=="Nsize" &
                  ts.data$year %in% rows2plot,]$value
  wavept <- wt(cbind(1:num_rows2plt,(ts-mean(ts))/sd(ts)), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  rm(ts)
  plot(wavept)
  # Plot frequency response
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="Nsize" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
  dev.off()
  
  
} # end plot_gen_freq_wvlt
