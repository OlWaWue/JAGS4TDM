library('rjags')
library('ggplot2')
library('readxl')
library('shiny')
library('gridExtra')


pk_1cmt_oral_ss<- function(theta, eta, dosing_events, times){
  
  ii <- dosing_events[,3]
  amt <- dosing_events[,2]
  
  
  
  k <- theta[3]*exp(eta[3])
  Vd <- theta[2]*exp(eta[2])
  ka <- theta[1]*exp(eta[1])
  f_oral <- theta[4]*exp(eta[4])
  t_lag <- theta[5]

  IPRED <- vector()
  
  for(i in 1:length(times))
    IPRED[i] <- ifelse(times[i] < t_lag, 
                    f_oral*amt/Vd *(ka/(ka-k))*( (exp(-k*(times[i]+ii-t_lag) )/(1-exp(-k*ii))) - (exp(-ka*(times[i]+ii-t_lag) )/(1-exp(-ka*ii))) ), 
                    f_oral*amt/Vd *(ka/(ka-k))*( (exp(-k*(times[i]-t_lag) )/(1-exp(-k*ii))) - (exp(-ka*(times[i]-t_lag) )/(1-exp(-ka*ii))) ) )
    
  return(IPRED)
}


### PK Model oral 1cmt with lag time
pk_1cmt_oral <- function(theta, eta, dosing_events, times){
  
  dosing_time <- dosing_events[,1]
  amt <- dosing_events[,2]
  
  k <- theta[3]*exp(eta[3])
  Vd <- theta[2]*exp(eta[2])
  ka <- theta[1]*exp(eta[1])
  F_oral <- theta[4]*exp(eta[4])
  t_lag <- theta[5]
  
  IPRED <- vector()
  
  for(t in 1:length(times)){
    temp_conc <- c(nrow = length(amt))
    for(i in 1:length(amt)) {
      
      temp_conc[i] <- ifelse((times[t]-dosing_time[i]) < t_lag, 0, F_oral*amt[i]*ka/(Vd*(ka-k))*(exp(-k*(times[t]-dosing_time[i]-t_lag))-exp(-ka*(times[t]-dosing_time[i]-t_lag))))
      
    }
    IPRED[t] <- sum(temp_conc)
  }
  
  return(IPRED)
}
## NONMEM like dataset
process_data_set <- function(pk_data = data.frame(time=c(0,4,6,12,30,50),
                                                   amt=c(100,".",".",100,".","."),
                                                   conc=c(".", 2, 3, ".", 1.5, 0.72),
                                                   evid=c(1, 0, 0, 1, 0, 0)),
                             n.burn=200, n.iter=1000, 
                             thetas  = c(0.5,
                                          60,
                                          0.09,
                                          0.9,
                                          0.5),
                             omegas = c(0.51,0.54,0.53, 0.53),
                             TIME = seq(0, 72, by=1), sigma=1.15, steady_state=F, n.comp =1) {

    do_plot <- function(jags_result, nburn=n.burn){
      
      df <- as.data.frame(jags_result[[4]])
      
      df <- df[-(1:nburn),]
      
      mcmc_se <- list()
      withProgress(message = "Processing MCMC results", max = nrow(df), {
          for (i in 1:nrow(df)) {
            if(!steady_state & n.comp==1){
                mcmc_se[[i]] <- pk_1cmt_oral(theta=thetas,
                                             eta=c(df$eta1[i],df$eta2[i],df$eta3[i], df$eta4[i]),
                                             dosing_events = dosing_events,
                                             times=TIME)
            } else {
                mcmc_se[[i]] <- pk_1cmt_oral_ss(theta=thetas,
                                             eta=c(df$eta1[i],df$eta2[i],df$eta3[i], df$eta4[i]),
                                             dosing_events = data.frame(time=0, amt=dosing_events$amt[1], ii=dosing_events$ii[1]),
                                             times=TIME)
            }
            
            incProgress(1)
          }
       
      })
      
      
      df_temp <- NULL
      
      for (k in 1:nrow(df)) {
        df_temp <- rbind(df_temp, mcmc_se[[k]])
        
      }
      
      
      s <- apply(df_temp,2,function(x) quantile(x,probs=c(0.05, 0.10, 0.15, 0.20, 0.8, 0.85, 0.9, 0.95, 0.5)))
      
      
      
      
      
    
        pk_data <- data.frame(time=TIME,
                              s1=s[1,],s2=s[8,], 
                              s3=s[2,],s4=s[7,],
                              s5=s[3,],s6=s[6,],
                              s7=s[4,],s8=s[5,],
                              max=s[9,])
  
      
      
      p <- ggplot(pk_data) + 
        geom_ribbon(aes(ymin=s1, ymax=s2, x=time), fill="red", alpha=0.15) + 
        geom_ribbon(aes(ymin=s3, ymax=s4, x=time), fill="red", alpha=0.15) + 
        geom_ribbon(aes(ymin=s5, ymax=s6, x=time), fill="red", alpha=0.15) + 
        geom_ribbon(aes(ymin=s7, ymax=s8, x=time), fill="red", alpha=0.15) + 
        geom_line(aes(y=max, x=time))+
        geom_point(data=tdm_data, aes(x=time, y=conc))
      
      c_at_tlast <- df_temp[,ncol(df_temp)]
      
      ind_y_max <- max(pk_data$s2)
      
      return(list(p, c_at_tlast, ind_y_max))
    }
    
    
    
    ## todo: detect number of etas automatically
    mcmc_diagnosticplots <- function(chain=1, jags_result, nburn = n.burn, omega, colour="red") {
      
      df <- as.data.frame(jags_result[[chain]]) ## Dataframe with etas
      
      nb.etas <- ncol(jags_result[[chain]]) 
      
      niter <- nrow(df)
      
      df$iteration <- (1:niter) 
      
      df <- df[-(1:nburn),]
      

      ## I will never again understand this...
      
      
      ## Problems that COULD happen: never refer aesthetics to columns via df[,n.et]
      ## Do not pass data to ggplot, but to the geom layer => problems with additional
      ## layers that are supposed to use fixed values like geom_vline
      plot_list <- list()
      
      for(n.et in 1:nb.etas) {
        
        dens_post <- density(df[,n.et], adjust=2)
        max_eta <- dens_post$x[which.max(dens_post$y)]
        dens_post <- data.frame(ETA=dens_post$x, freq=dens_post$y, max_eta=max_eta)
        x_this_eta <- seq(min(dens_post$ETA), max(dens_post$ETA), length.out = 512)
        
        dens_prior <- dnorm(x=x_this_eta, mean=0, sd=sqrt(omega[n.et]))
        dens_prior <- data.frame(ETA=x_this_eta, freq=dens_prior)
        
        current_data <- data.frame(iteration=(nburn+1):niter,eta=df[,n.et])
        
        y_name <- paste("ETA", n.et)
        
        ## trace plot for this eta in this chain
        plot_list[[paste("p_iter_ETA", n.et, sep="")]] <- ggplot(data=current_data) + 
          geom_line(aes(x=iteration,y=eta), colour=colour)+ 
          ylim(min(x_this_eta),
               max(x_this_eta)) +
          theme_bw() + ylab(y_name) + xlab("Iteration")
        
        ## 90° flipped density plot for this eta in this chain
        plot_list[[paste("p_dens_ETA", n.et, sep="")]] <- ggplot(data = dens_post) +
          geom_vline(aes(xintercept=max_eta), colour="red", linetype=2) +
          geom_line(aes(x=ETA, y=freq), colour="red") + 
          geom_line(data=dens_prior, aes(x=ETA,y=freq), colour = "blue") + 
          coord_flip() + 
          xlim(min(x_this_eta),
               max(x_this_eta))+
          annotate(geom="text", 
                   y=max(dens_post$freq)+0.2, 
                   x=max_eta, 
                   label="posterior", colour="red") +
          annotate(geom="text", y=max(dens_prior$freq)+0.2, x=0, label="prior", colour="blue") +
          ylim(0,2)+ ylab(paste("Density of", y_name)) +
          theme_bw() + theme(axis.title.y = element_blank()) 
        
     
      }
    

      
      return(plot_list)
    }
    
    
    
    
    dosing_events <- data.frame(time=as.numeric(as.character(pk_data[pk_data$evid==1,]$time)),
                                amt=as.numeric(as.character(pk_data[pk_data$evid==1,]$amt)),
                                ii=as.numeric(as.character(pk_data[pk_data$evid==1,]$ii)))
    
    
    
    tdm_data <- data.frame(conc=as.numeric(as.character(pk_data[pk_data$evid==0,]$conc)),
                           time=as.numeric(as.character(pk_data[pk_data$evid==0,]$time)))
    


    if((n.comp==1) & (!steady_state))  {
    
        jags <- jags.model('1cmt_multiple_dose.bug',
                           data = list('c' = tdm_data$conc,
                                       'amt' = dosing_events$amt, 
                                       'dosing_time' = dosing_events$time,
                                       'ts'= tdm_data$time,
                                       "theta"=thetas,
                                       "omega"=sqrt(omegas),
                                       'sigma'=sqrt(sigma) ),
                           n.chains = 4,
                           n.adapt = n.iter)
        
        d <- coda.samples(jags,
                          c('eta1', 'eta2', 'eta3', 'eta4'),
                          n.iter, thin=1)
    } else if((n.comp==1) & (steady_state)) {
      
          jags <- jags.model('1cmt_ss.bug',
                             data = list('c' = tdm_data$conc,
                                         'amt' = dosing_events$amt[1], 
                                         'ii' = dosing_events$ii[1],
                                         'times'= tdm_data$time,
                                         "theta"=thetas,
                                         "omega"=sqrt(omegas),
                                         'sigma'=sqrt(sigma) ),
                             n.chains = 4,
                             n.adapt = n.iter)
          
          d <- coda.samples(jags,
                            c('eta1', 'eta2', 'eta3', 'eta4'),
                            n.iter, thin=1)
    }
    
    
    
    pk_profile <- (do_plot(d))
    
    
    mcmc_plots_1 <- mcmc_diagnosticplots(1, d, nburn=n.burn, omega=omegas, "red")
    mcmc_plots_2 <- mcmc_diagnosticplots(2, d, nburn=n.burn, omega=omegas, "orange")
    mcmc_plots_3 <- mcmc_diagnosticplots(3, d, nburn=n.burn, omega=omegas, "yellow")
    mcmc_plots_4 <- mcmc_diagnosticplots(4, d, nburn=n.burn, omega=omegas)
    
    
    
    
      result = list(mcmc_plots_1, mcmc_plots_2, mcmc_plots_3, mcmc_plots_4, pk_profile=pk_profile[[1]], c_at_tlast=pk_profile[[2]], ind_y_max=pk_profile[[3]])
      
      return(result)
}

