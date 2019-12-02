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

plot_gen_freq_wvlt <- function( ) {
  burn_in_pd = rm_first_timesteps #remove from beginning of ts
  num_rows2plt = 1000 #number of years to plot
  rows2plot <- (burn_in_pd+1):(burn_in_pd+num_rows2plt)
  n = 1 #for indexing column n in noise df?
  J1 = trunc((log(32/(2 * 1))/log(2))/0.01)
  timesteps = timesteps
  span.multiplier = 1
  noise = noisetypes[[1]]
  noisename = noisenames[1]
  # ---
  test = subset(parms[parms$problem_spp=="no",],select=c("spp","maxage"))
  test$maxage <- as.numeric(test$maxage)
  test <- test[order(test$maxage),]
  spp <- test$spp #put spp in order of max age
  #df.list <- df.list[spp] #set order
  line_d <- 2 #which margin to start counting at zero, 2=left
  tsylim <- c(min(ts.data[ts.data$year %in% rows2plot,]$value),
              max(ts.data[ts.data$year %in% rows2plot,]$value))
  spwn_ylim <- c(0,0.35) #ylimits for spawning biomass dist
  spwn_xlim <- c(0,110)  #xlimits for spawning biomass dist
  
  n_rows <- seq(from=1,to=length(rows2plot))
  old <- par(mar = c(3,2,2,1), cex = .7)
  specminval <- as.numeric(min(ts.spec$spec))
  specmaxval <- as.numeric(max(ts.spec$spec))
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
  
  tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/wvlet_panels_1-7_ENSO_v2.tiff', units="in", width=8, height=11, res=300) 
  
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
  #mtext("a", side = 2, las = 1, at = max(whitenoise), line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  noise.wt <- wt(cbind(1:num_rows2plt, noise[rows2plot]), 
                 dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
                 sig.test = 0, sig.level = 0.95)
  plot(noise.wt)
  #mtext("b", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  # Plot frequency response
  #m = 55
  p_spec <- spec.pgram(x=noise,
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
  #par(fig=c(0.1, 0.4, 0.80, 1), new = TRUE)
  plot(x = n_rows,
       y = ts.data[ts.data$spp == spp[i] &
                   ts.data$EI==0 &
                   ts.data$noise==noisename &
                   ts.data$variable=="recruits" &
                   ts.data$year %in% rows2plot,]$value,
       type = "l", ylim=tsylim, xaxs = "i", yaxt="n", cex=0.8, 
       ylab="",log="y")
  title(ylab="log(recruit)",line=0)
  
  #mtext("d", side = 2, las = 1, at = max(tsylim), line = line_d, cex = 1.2)
  
  # plot wavelet power spectrum
  #par(fig=c(0.4, 0.7, 0.80, 1), new = TRUE)
  wavept <- wt(cbind(1:num_rows2plt,ts.data[ts.data$spp==spp[i] & 
                         ts.data$EI==0 & 
                         ts.data$noise==noisename &
                         ts.data$variable=="recruits" &
                         ts.data$year %in% rows2plot,]$value), 
               dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", 
               sig.test = 0, sig.level = 0.95)
  plot(wavept)
  #mtext("e", side = 2, las = 1, at = 1, line = line_d-1, cex = 1.2)
  
  # Plot frequency response
  #par(fig=c(0.7, 1, 0.8, 1), new = TRUE)
  #m = as.numeric(parms[parms$spp == parms$spp[i],]$maxage)
  plot(x=ts.spec[ts.spec$spp==spp[i] &
                 ts.spec$noise==noisename &
                 ts.spec$variable=="recruits" &
                 ts.spec$EI==0,]$freq,
       y=ts.spec[ts.spec$spp==spp[i] &
                   ts.spec$noise==noisename &
                   ts.spec$variable=="recruits" &
                   ts.spec$EI==0,]$spec,
       type = "l",ylab="",xlab="frequency",log="y",yaxt="n")
  title(ylab="log(power)",line=0)
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
  axis(2,at=c(10,100,1000,10000))
  #mtext("u", side = 2, las = 1, at = 0.5, line = line_d-1, cex = 1.2)
  
  dev.off()
  
  # --------
  tiff(file='C:/Users/Mikaela/Documents/GitHub/pfmc/results/wvlet_panels_8-14_ENSO_v1.tiff', units="in", width=8, height=11, res=300) 
  
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
       type = "l",log="y",yaxt="n",ylab="",ylim(c(specminval,specmaxval)),xlab="frequency")
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
