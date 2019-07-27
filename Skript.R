library('rjags')
library('ggplot2')

dosing_events <- data.frame(time=c(0,12),
                            amt=c(100,100))

tdm_data <- data.frame(conc=c(2,3, 1.5, 0.72),
                       time=c(4,6, 30, 50))

thetas  <- c(0.5,
             60,
             0.09,
             0.5)

omegas <- c(0.51,0.54,0.53)

TIME <- seq(0, 72, by=0.5)

pk_1cmt_oral <- function(theta, eta, dosing_events, times){
  
  dosing_time <- dosing_events[,1]
  amt <- dosing_events[,2]
  
    k <- theta[3]*exp(eta[3])
    v_F <- theta[2]*exp(eta[2])
    ka <- theta[1]*exp(eta[1])
    t_lag <- theta[4]

    IPRED <- vector()
    
    for(t in 1:length(times)){
      temp_conc <- c(nrow = length(amt))
      for(i in 1:length(amt)) {
        
        temp_conc[i] <- ifelse((times[t]-dosing_time[i]) < t_lag, 0, amt[i]*ka/(v_F*(ka-k))*(exp(-k*(times[t]-dosing_time[i]-t_lag))-exp(-ka*(times[t]-dosing_time[i]-t_lag))))
        
      }
      IPRED[t] <- sum(temp_conc)
    }
    
    return(IPRED)
}
  
IPRED <- pk_1cmt_oral(theta=thetas,
                       eta=c(0,0,0),
                       dosing_events = dosing_events,
                       times=TIME)



pk_data <- data.frame(TIME, IPRED)

ggplot(pk_data, aes(x=TIME, y=IPRED)) + geom_line() + geom_point(data=tdm_data,aes(x=time, y=conc))

jags <- jags.model('1cmt_multiple_dose.bug',
                   data = list('c' = tdm_data$conc,
                               'amt' = dosing_events$amt, 
                               'dosing_time' = dosing_events$time,
                               't_lag' = thetas[4],
                               'ts'= tdm_data$time,
                               "theta"=thetas,
                               "omega"=omegas),
                   n.chains = 4,
                   n.adapt = 5000)
d <- coda.samples(jags,
                  c('eta1', 'eta2', 'eta3'),
                  5000, thin=1)

plot(do_plot(d))


mcmc_plots_1 <- mcmc_diagnosticplots(1, d, nburn=500, omega=omegas, "red")
mcmc_plots_2 <- mcmc_diagnosticplots(2, d, nburn=500, omega=omegas, "orange")
mcmc_plots_3 <- mcmc_diagnosticplots(3, d, nburn=500, omega=omegas, "yellow")
mcmc_plots_4 <- mcmc_diagnosticplots(4, d, nburn=500, omega=omegas)


gridExtra::grid.arrange(mcmc_plots_1$p_iter_ETA1, mcmc_plots_1$p_dens_ETA1,   
                        mcmc_plots_1$p_iter_ETA2, mcmc_plots_1$p_dens_ETA2,
                        mcmc_plots_1$p_iter_ETA3, mcmc_plots_1$p_dens_ETA3, nrow=3, ncol=2)

gridExtra::grid.arrange(mcmc_plots_2$p_iter_ETA1, mcmc_plots_2$p_dens_ETA1,   
                        mcmc_plots_2$p_iter_ETA2, mcmc_plots_2$p_dens_ETA2,
                        mcmc_plots_2$p_iter_ETA3, mcmc_plots_2$p_dens_ETA3, nrow=3, ncol=2)

gridExtra::grid.arrange(mcmc_plots_3$p_iter_ETA1, mcmc_plots_3$p_dens_ETA1,   
                        mcmc_plots_3$p_iter_ETA2, mcmc_plots_3$p_dens_ETA2,
                        mcmc_plots_3$p_iter_ETA3, mcmc_plots_3$p_dens_ETA3, nrow=3, ncol=2)

gridExtra::grid.arrange(mcmc_plots_4$p_iter_ETA1, mcmc_plots_4$p_dens_ETA1,   
                        mcmc_plots_4$p_iter_ETA2, mcmc_plots_4$p_dens_ETA2,
                        mcmc_plots_4$p_iter_ETA3, mcmc_plots_4$p_dens_ETA3, nrow=3, ncol=2)


do_plot <- function(jags_result, nburn=500){

    df <- as.data.frame(jags_result[[4]])

    df <- df[-(1:nburn),]
    
    mcmc_se <- list()

      for (i in 1:nrow(df)) {
        mcmc_se[[i]] <- pk_1cmt_oral(theta=thetas,
                                     eta=c(df$eta1[i],df$eta2[i],df$eta3[i]),
                                     dosing_events = dosing_events,
                                     times=TIME)
        
          
      }
    
    df_temp <- NULL
    
    for (k in 1:1000) {
      df_temp <- rbind(df_temp, mcmc_se[[k]])
    
    }
    
    dens <- density(df$eta1, adjust=2)
    max_eta1 <- dens$x[which.max(dens$y)]
    
    dens <- density(df$eta2, adjust=2)
    max_eta2 <- dens$x[which.max(dens$y)]
    
    dens <- density(df$eta3, adjust=2)
    max_eta3 <- dens$x[which.max(dens$y)]
    
    
    s <- apply(df_temp,2,function(x) quantile(x,probs=c(0.05, 0.10, 0.15, 0.20, 0.8, 0.85, 0.9, 0.95)))
    

    
    pk_data <- data.frame(time=TIME,
                          s1=s[1,],s2=s[8,], 
                          s3=s[2,],s4=s[7,],
                          s5=s[3,],s6=s[6,],
                          s7=s[4,],s8=s[5,],
                          max=pk_1cmt_oral(theta=thetas,
                                           eta=c(max_eta1,max_eta2,max_eta3),
                                           dosing_events = dosing_events,
                                           times=TIME),
                          avg=pk_1cmt_oral(theta=thetas,
                                           eta=c(0,0,0),
                                           dosing_events = dosing_events,
                                           times=TIME))
                            
    print(pk_data)
    
    p <- ggplot(pk_data) + 
      geom_ribbon(aes(ymin=s1, ymax=s2, x=time), fill="red", alpha=0.25) + 
      geom_ribbon(aes(ymin=s3, ymax=s4, x=time), fill="red", alpha=0.25) + 
      geom_ribbon(aes(ymin=s5, ymax=s6, x=time), fill="red", alpha=0.25) + 
      geom_ribbon(aes(ymin=s7, ymax=s8, x=time), fill="red", alpha=0.25) + 
      geom_line(aes(y=avg, x=time), linetype=2) +
          geom_line(aes(y=max, x=time))+
          geom_point(data=tdm_data, aes(x=time, y=conc))
    
    return(p)
}



## todo: detect number of etas automatically
mcmc_diagnosticplots <- function(chain=1, jags_result, nburn = 500, omega, colour="blue") {
  
  df <- as.data.frame(jags_result[[chain]]) ## Dataframe with etas
  
  
  niter <- nrow(df)
  
  df$iteration <- (1:niter) 
  
  df <- df[-(1:nburn),]
  
  dk_ETA1 <- density(df$eta1, adjust=2)
  dk_ETA1 <- data.frame(ETA1=dk_ETA1$x,freq=dk_ETA1$y)
  xk_ETA1 <- seq(min(dk_ETA1$ETA1),max(dk_ETA1$ETA1), length.out = 512)
  
  dk_ETA2 <- density(df$eta2, adjust=2)
  dk_ETA2 <- data.frame(ETA2=dk_ETA2$x,freq=dk_ETA2$y)
  xk_ETA2 <- seq(min(dk_ETA2$ETA2),max(dk_ETA2$ETA2), length.out = 512)
  
  dk_ETA3 <- density(df$eta3, adjust=2)
  dk_ETA3 <- data.frame(ETA3=dk_ETA3$x,freq=dk_ETA3$y)
  xk_ETA3 <- seq(min(dk_ETA3$ETA3),max(dk_ETA3$ETA3), length.out = 512)
  

  ### Prior distributions
  dp_ETA1 <- dnorm(x = xk_ETA1, mean = 0, sd = sqrt(omega[1]) )
  dp_ETA1 <- data.frame(ETA1 = xk_ETA1, freq=dp_ETA1)
  
  dp_ETA2 <- dnorm(x = xk_ETA2, mean = 0, sd = sqrt(omega[2]) )
  dp_ETA2 <- data.frame(ETA2 = xk_ETA2, freq=dp_ETA2)
  
  dp_ETA3 <- dnorm(x = xk_ETA3, mean = 0, sd = sqrt(omega[3]) )
  dp_ETA3 <- data.frame(ETA3 = xk_ETA3, freq=dp_ETA3)

  
  
  
  
  
  #     p1 <- qplot((1:niter),ka,data=p,colour="blue")
  p1 <- ggplot(data=df) + geom_line(aes(x=iteration,y=eta1), colour=colour)+ ylim(min(dk_ETA1$ETA1),max(dk_ETA1$ETA1)) +
    theme_bw()
  
  p1dens <- ggplot(data = dk_ETA1) + geom_line(aes(x=ETA1, y=freq), colour="green") + 
    geom_line(data=dp_ETA1, aes(x=ETA1,y=freq), colour = "red") + 
    coord_flip() + xlim(min(dk_ETA1$ETA1),max(dk_ETA1$ETA1))+
    annotate(geom="text", y=max(dk_ETA1$freq)*1.1, x=dk_ETA1$ETA1[which.max(dk_ETA1$freq)], label="posterior", colour="green") +
    annotate(geom="text", y=max(dp_ETA1$freq)*1.1, x=0, label="prior", colour="red") +
    ylim(0,max(dk_ETA1$freq*1.20))+
    theme_bw() + theme(axis.title.y = element_blank())
  
  p2 <- ggplot(data=df) + geom_line(aes(x=iteration,y=eta2), colour=colour)+ ylim(min(dk_ETA2$ETA2),max(dk_ETA2$ETA2))+
    theme_bw()
  
  p2dens <- ggplot(data = dk_ETA2) + geom_line(aes(x=ETA2, y=freq), colour="green") + 
    geom_line(data=dp_ETA2, aes(x=ETA2,y=freq), colour = "red") + 
    coord_flip() + 
    annotate(geom="text", y=max(dk_ETA2$freq)*1.1, x=dk_ETA2$ETA2[which.max(dk_ETA2$freq)], label="posterior", colour="green") +
    annotate(geom="text", y=max(dp_ETA2$freq)*1.1, x=0, label="prior", colour="red") +
    ylim(0,max(dk_ETA2$freq*1.20))+
    xlim(min(dk_ETA2$ETA2),max(dk_ETA2$ETA2))+
    theme_bw() + theme(axis.title.y = element_blank())
  
  p3 <- ggplot(data=df) + geom_line(aes(x=iteration,y=eta3), colour=colour)+ ylim(min(dk_ETA3$ETA3),max(dk_ETA3$ETA3))+
    theme_bw() 
    
  p3dens <- ggplot(data = dk_ETA3) + geom_line(aes(x=ETA3, y=freq), colour="green") + 
    geom_line(data=dp_ETA3, aes(x=ETA3,y=freq), colour = "red") + 
    coord_flip() + 
    annotate(geom="text", y=max(dk_ETA3$freq)*1.1, x=dk_ETA3$ETA3[which.max(dk_ETA3$freq)], label="posterior", colour="green") +
    annotate(geom="text", y=max(dp_ETA3$freq)*1.1, x=0, label="prior", colour="red") +
    ylim(0,max(dk_ETA3$freq*1.20))+
    xlim(min(dk_ETA3$ETA3),max(dk_ETA3$ETA3))+
    theme_bw() + theme(axis.title.y = element_blank())

  
  
  ret = list(p_iter_ETA1 = p1, 
             p_dens_ETA1 = p1dens,
             p_iter_ETA2 = p2, 
             p_dens_ETA2 = p2dens,
             p_iter_ETA3 = p3, 
             p_dens_ETA3 = p3dens) 
  
  return(ret)
}


