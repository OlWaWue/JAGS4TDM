library('rjags')
library('ggplot2')
library('readxl')
library('shiny')
library('gridExtra')
library('PerformanceAnalytics')
library('MASS')

## Thetas and omega matrix from the Axitinib Final model in fed state using formulation IV
axi_i_mod_fed <- list(thetas=c(0.523, ## THETA1: ka in fed state formulation IV
                               45.3,  ## THETA2: Vc [L]
                               17.0,  ## THETA3: Cl [L/h]
                               0.465, ## THETA4: F formulation IV in fed state
                               0.457, ## THETA5: Lag time [h]
                               45.9,  ## THETA6: Vp [L]
                               1.74), ## THETA7: Q [L/h]
                       omegas=c(0.506, ## OMEGA1: Variance ka
                                0.0949, ## OMEGA2: Variance Vc
                                0.272, ## OMEGA3: Variance Cl
                                1.07,  ## OMEGA4: Variance Vp
                                0.406,  ## OMEGA5: Variance Q
                                0.141, ## OMEGA6: Covariance Cl~Vc
                                0.619)) ## OMEGA7: Covariance Q~Vp

## Thetas and omega matrix from the Axitinib Final model in !fasted! state using formulation IV
axi_i_mod_fasted <- list(thetas=c(0.523*(1+2.07), ## THETA1: ka in fed state formulation IV
                                  45.3,  ## THETA2: Vc [L]
                                  17.0,  ## THETA3: Cl [L/h]
                                  0.465*(1+0.338), ## THETA4: F formulation IV in fed state
                                  0.457, ## THETA5: Lag time [h]
                                  45.9,  ## THETA6: Vp [L]
                                  1.74), ## THETA7: Q [L/h]
                         omegas=c(0.506, ## OMEGA1: Variance ka
                                  0.0949, ## OMEGA2: Variance Vc
                                  0.272, ## OMEGA3: Variance Cl
                                  1.07,  ## OMEGA4: Variance Vp
                                  0.406,  ## OMEGA5: Variance Q
                                  0.141, ## OMEGA6: Covariance Cl~Vc
                                  0.619)) ## OMEGA7: Covariance Q~Vp

## Thetas and omega matrix from the Axitinib Final model in fed state using formulation XLI
axi_i_mod_fed_form_XLI <- list(thetas=c(0.523, ## THETA1: ka in fed state formulation IV
                                        45.3,  ## THETA2: Vc [L]
                                        17.0,  ## THETA3: Cl [L/h]
                                        0.465*(1-0.150), ## THETA4: F formulation IV in fed state
                                        0.457, ## THETA5: Lag time [h]
                                        45.9,  ## THETA6: Vp [L]
                                        1.74), ## THETA7: Q [L/h]
                               omegas=c(0.506, ## OMEGA1: Variance ka
                                        0.0949, ## OMEGA2: Variance Vc
                                        0.272, ## OMEGA3: Variance Cl
                                        1.07,  ## OMEGA4: Variance Vp
                                        0.406,  ## OMEGA5: Variance Q
                                        0.141, ## OMEGA6: Covariance Cl~Vc
                                        0.619)) ## OMEGA7: Covariance Q~Vp

## Analytical solution of 2cmt model in steady state with lag time
pk_2cmt_oral_ss <- function(theta, eta, dosing_events, times){
  ii <- dosing_events[,3]
  amt <- dosing_events[,2]
  
  dosing_time <- dosing_events[,1]
  
  Cl <- theta[3]*exp(eta[3])
  V1 <- theta[2]*exp(eta[2])
  ka <- theta[1]*exp(eta[1])
  f_oral <- theta[4]*exp(eta[4])
  V2 <- theta[6] * exp(eta[5])
  Q <- theta[7] * exp(eta[6])
  t_lag <- theta[5]
  
  k=Cl/V1
  k12 = Q/V1
  k21 = Q/V2
  
  
  beta = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k)^2 - 4 * k21 * k))
  
  alpha = (k21 * k) / beta
  
  A = ka/V1 * (k21-alpha)/((ka-alpha) * (beta-alpha))
  B = ka/V1 * (k21-beta)/((ka-beta) * (alpha-beta))
  
  IPRED <- vector()
  
  for(i in 1:length(times))
    IPRED[i] <- ifelse(times[i] < t_lag, 
                       f_oral*amt * ( (A*exp(-alpha*(times[i]-dosing_time+ii-t_lag)))/(1-exp(-alpha*ii)) + (B * exp(-beta*(times[i]-dosing_time+ii-t_lag)))/(1-exp(-beta*ii))  - ((A+B)*exp(-ka*(times[i]-dosing_time+ii-t_lag)))/(1-exp(-ka*ii)) ), 
                       f_oral*amt * ( (A*exp(-alpha*(times[i]-dosing_time-t_lag)))/(1-exp(-alpha*ii)) + (B * exp(-beta*(times[i]-dosing_time-t_lag)))/(1-exp(-beta*ii))  - ((A+B)*exp(-ka*(times[i]-dosing_time-t_lag)))/(1-exp(-ka*ii)) )  )
  
  return(IPRED)
  
}

## analytical solution of 2cmt model multiple doses with lag time
pk_2cmt_oral <- function(theta, eta, dosing_events, times){
  
  dosing_time <- dosing_events[,1]
  amt <- dosing_events[,2]
  
  
  Cl <- theta[3]*exp(eta[3])
  V1 <- theta[2]*exp(eta[2])
  ka <- theta[1]*exp(eta[1])
  f_oral <- theta[4]*exp(eta[4])
  V2 <- theta[6] * exp(eta[5])
  Q <- theta[7] * exp(eta[6])
  t_lag <- theta[5]
  
  k=Cl/V1
  k12 = Q/V1
  k21 = Q/V2
  
  
  beta = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k)^2 - 4 * k21 * k))
  
  alpha = (k21 * k) / beta
  
  A = ka/V1 * (k21-alpha)/((ka-alpha) * (beta-alpha))
  B = ka/V1 * (k21-beta)/((ka-beta) * (alpha-beta))
  
  IPRED <- vector()
  
  for(t in 1:length(times)){
    temp_conc <- c(nrow = length(amt))
    for(i in 1:length(amt)) {
      
      temp_conc[i] <- ifelse((times[t]-dosing_time[i]) <= t_lag,
                             0 , 
                             f_oral*amt[i]*(A*exp(-alpha*(times[t]-dosing_time[i]-t_lag)) + B * exp(-beta*(times[t]-dosing_time[i]-t_lag)) - (A+B)*exp(-ka*(times[t]-dosing_time[i]-t_lag))) ) 
      
    }
    IPRED[t] <- sum(temp_conc)
  }
  return(IPRED)
}

## 1cmt model in ss with lag time
pk_1cmt_oral_ss<- function(theta, eta, dosing_events, times){
  
  ii <- dosing_events[,3]
  amt <- dosing_events[,2]
  
  dosing_time <- dosing_events[,1]
  
  Cl <-theta[3]*exp(eta[3])
  
 
  Vd <- theta[2]*exp(eta[2])
  ka <- theta[1]*exp(eta[1])
  f_oral <- theta[4]*exp(eta[4])
  t_lag <- theta[5]

  k <- Cl/Vd
  
  IPRED <- vector()
  
  for(i in 1:length(times))
    IPRED[i] <- ifelse(times[i] < t_lag, 
                    f_oral*amt/Vd *(ka/(ka-k))*( (exp(-k*(times[i]-dosing_time+ii-t_lag) )/(1-exp(-k*ii))) - (exp(-ka*(times[i]-dosing_time+ii-t_lag) )/(1-exp(-ka*ii))) ), 
                    f_oral*amt/Vd *(ka/(ka-k))*( (exp(-k*(times[i]-t_lag) )/(1-exp(-k*ii))) - (exp(-ka*(times[i]-t_lag) )/(1-exp(-ka*ii))) ) )
    
  return(IPRED)
}


### PK Model oral 1cmt with lag time
pk_1cmt_oral <- function(theta, eta, dosing_events, times){
  
  dosing_time <- dosing_events[,1]
  amt <- dosing_events[,2]
  
  Cl <- theta[3]*exp(eta[3])
  Vd <- theta[2]*exp(eta[2])
  ka <- theta[1]*exp(eta[1])
  F_oral <- theta[4]*exp(eta[4])
  t_lag <- theta[5]
  
  k <- Cl/Vd
  
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

## Process a NONMEM like dataset
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

  # inner function to generate individual PK plot
    do_plot <- function(jags_result, nburn=n.burn){
      
      ## Use the fourth chain for the simulation
      df <- as.data.frame(jags_result[[4]])
      
      ## remove burnin iterations
      df <- df[-(1:nburn),]
      
      mcmc_se <- list()
      
      ## Simulate with the etas obtained by sampling from the posterior distribution
      withProgress(message = "Processing MCMC results", max = nrow(df), {
          for (i in 1:nrow(df)) {
            if(!steady_state & n.comp==1){
                mcmc_se[[i]] <- pk_1cmt_oral(theta=thetas,
                                             eta=c(df$eta1[i],df$eta2[i],df$eta3[i], df$eta4[i]),
                                             dosing_events = dosing_events,
                                             times=TIME)
            } else if (steady_state & n.comp==1) {
                mcmc_se[[i]] <- pk_1cmt_oral_ss(theta=thetas,
                                             eta=c(df$eta1[i],df$eta2[i],df$eta3[i], df$eta4[i]),
                                             dosing_events = data.frame(time=0, amt=dosing_events$amt[1], ii=dosing_events$ii[1]),
                                             times=TIME)
            } else if (!steady_state & n.comp==2) {
              mcmc_se[[i]] <- pk_2cmt_oral(theta=thetas,
                                              eta=c(df$eta1[i],df$eta2[i],df$eta3[i], df$eta4[i], df$eta5[i], df$eta6[i]),
                                              dosing_events = dosing_events,
                                              times=TIME)
            } else if (steady_state & n.comp==2) {
              mcmc_se[[i]] <- pk_2cmt_oral_ss(theta=thetas,
                                              eta=c(df$eta1[i],df$eta2[i],df$eta3[i], df$eta4[i], df$eta5[i], df$eta6[i]),
                                              dosing_events = data.frame(time=0, amt=dosing_events$amt[1], ii=dosing_events$ii[1]),
                                              times=TIME)
            } else if (!steady_state & n.comp==3) { ## Axitinib Fed
              mcmc_se[[i]] <- pk_2cmt_oral(theta=thetas, ### ETA 4 is set to zero => no random effect f_oral
                                           eta=c(df$eta1[i],df$eta2[i],df$eta3[i],0, df$eta4[i], df$eta5[i]),
                                           dosing_events = dosing_events,
                                           times=TIME)
            } else if (steady_state & n.comp==3) { ## Axitinib Fed
              mcmc_se[[i]] <- pk_2cmt_oral_ss(theta=thetas, ### ETA 4 is set to zero => no random effect f_oral
                                              eta=c(df$eta1[i],df$eta2[i],df$eta3[i], 0, df$eta4[i], df$eta5[i]),
                                              dosing_events = data.frame(time=0, amt=dosing_events$amt[1], ii=dosing_events$ii[1]),
                                              times=TIME)
            } else if (steady_state & n.comp==4) { ## Axitinib Fasted
              mcmc_se[[i]] <- pk_2cmt_oral_ss(theta=thetas, ### ETA 4 is set to zero => no random effect f_oral
                                              eta=c(df$eta1[i],df$eta2[i],df$eta3[i], 0, df$eta4[i], df$eta5[i]),
                                              dosing_events = data.frame(time=0, amt=dosing_events$amt[1], ii=dosing_events$ii[1]),
                                              times=TIME)
            } else if (steady_state & n.comp==5) { ## Axitinib Fasted formulation XLI
              mcmc_se[[i]] <- pk_2cmt_oral_ss(theta=thetas, ### ETA 4 is set to zero => no random effect f_oral
                                              eta=c(df$eta1[i],df$eta2[i],df$eta3[i], 0, df$eta4[i], df$eta5[i]),
                                              dosing_events = data.frame(time=0, amt=dosing_events$amt[1], ii=dosing_events$ii[1]),
                                              times=TIME)
            } else if (!steady_state & n.comp==4) { ## Axitinib Fasted
              mcmc_se[[i]] <- pk_2cmt_oral(theta=thetas, ### ETA 4 is set to zero => no random effect f_oral
                                           eta=c(df$eta1[i],df$eta2[i],df$eta3[i],0, df$eta4[i], df$eta5[i]),
                                           dosing_events = dosing_events,
                                           times=TIME)
            } else if (!steady_state & n.comp==5) { ## Axitinib Fasted formulation XLI
              mcmc_se[[i]] <- pk_2cmt_oral(theta=thetas, ### ETA 4 is set to zero => no random effect f_oral
                                           eta=c(df$eta1[i],df$eta2[i],df$eta3[i],0, df$eta4[i], df$eta5[i]),
                                           dosing_events = dosing_events,
                                           times=TIME)
            }
            
            incProgress(1)
          }
       
      })
      
      
      df_temp <- NULL
      
      withProgress(message = "Calculating Prediction interval", max = nrow(df), {
          ## bind simulations in a data.frame
          for (k in 1:nrow(df)) {
            df_temp <- rbind(df_temp, mcmc_se[[k]])
            incProgress(1)
          }
      })
      
      ## Generate Quantils
      s <- apply(df_temp,2,function(x) quantile(x,probs=c(0.05, 0.10, 0.15, 0.20, 0.8, 0.85, 0.9, 0.95, 0.5)))

      ## Combine individual PK data and quantils in a data.frame
        pk_data <- data.frame(time=TIME,
                              s1=s[1,],s2=s[8,], 
                              s3=s[2,],s4=s[7,],
                              s5=s[3,],s6=s[6,],
                              s7=s[4,],s8=s[5,],
                              max=s[9,]) # median 
  
      
      ## Build raw individual PK plot
      p <- ggplot(pk_data) + 
        geom_ribbon(aes(ymin=s1, ymax=s2, x=time), fill="red", alpha=0.15) + 
        geom_ribbon(aes(ymin=s3, ymax=s4, x=time), fill="red", alpha=0.15) + 
        geom_ribbon(aes(ymin=s5, ymax=s6, x=time), fill="red", alpha=0.15) + 
        geom_ribbon(aes(ymin=s7, ymax=s8, x=time), fill="red", alpha=0.15) + 
        geom_line(aes(y=max, x=time))+
        geom_point(data=tdm_data, aes(x=time, y=conc))
      
      ## Get the last simulated concentration for the boxplot
      c_at_tlast <- df_temp[,ncol(df_temp)]
      
      ## carry on max value in PK plot for adjusting y-axis
      ind_y_max <- max(pk_data$s2)
      
      ind_y_min <- min(pk_data$s1[pk_data$s1>0])
      
      return(list(p, c_at_tlast, ind_y_max, ind_y_min))
    }

    ### Generate MCMC trace plots and distributions for every eta
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
      
      df <- df[,-ncol(df)]
      
      plot_list[["chain_data"]] <- df
      
      
      
      return(plot_list)
    }
    
    
    
    ### separate dosing events ...
    dosing_events <- data.frame(time=as.numeric(as.character(pk_data[pk_data$evid==1,]$time)),
                                amt=as.numeric(as.character(pk_data[pk_data$evid==1,]$amt)),
                                ii=as.numeric(as.character(pk_data[pk_data$evid==1,]$ii)))
    
    
    ### ... and TDM events
    tdm_data <- data.frame(conc=as.numeric(as.character(pk_data[pk_data$evid==0,]$conc)),
                           time=as.numeric(as.character(pk_data[pk_data$evid==0,]$time)))
    

    ### Use JAGS model to sample from posterior distribution
    ### According to the selected model
    
    #### ---- For some reason, JAGS runs more stable when SD is submitted to the model and is ----#
    #### Calculated back to variance in the model file --- ?? 
    if((n.comp==1) & (!steady_state))  {
    
        jags <- jags.model('1cmt_multiple_dose.bug',
                           data = list('c' = tdm_data$conc,
                                       'amt' = dosing_events$amt, 
                                       'dosing_time' = dosing_events$time,
                                       'times'= tdm_data$time,
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
                                         'dosing_time' = dosing_events$time[1],
                                         'times'= tdm_data$time,
                                         "theta"=thetas,
                                         "omega"=sqrt(omegas),
                                         'sigma'=sqrt(sigma) ),
                             n.chains = 4,
                             n.adapt = n.iter)
          
          d <- coda.samples(jags,
                            c('eta1', 'eta2', 'eta3', 'eta4'),
                            n.iter, thin=1)
    } else if((n.comp==2) & (!steady_state)) {
      
      jags <- jags.model('2cmt_multiple_dose.bug',
                         data = list('c' = tdm_data$conc,
                                     'amt' = dosing_events$amt, 
                                     'dosing_time' = dosing_events$time,
                                     'times'= tdm_data$time,
                                     "theta"=thetas,
                                     "omega"=sqrt(omegas),
                                     'sigma'=sqrt(sigma) ),
                         n.chains = 4,
                         n.adapt = n.iter)
      
      d <- coda.samples(jags,
                        c('eta1', 'eta2', 'eta3', 'eta4', 'eta5', 'eta6'),
                        n.iter, thin=1)
    } else if((n.comp==2) & (steady_state)) {
      jags <- jags.model('2cmt_ss.bug',
                         data = list('c' = tdm_data$conc,
                                     'amt' = dosing_events$amt[1], 
                                     'ii' = dosing_events$ii[1],
                                     'dosing_time' = dosing_events$time[1],
                                     'times'= tdm_data$time,
                                     "theta"=thetas,
                                     "omega"=sqrt(omegas),
                                     'sigma'=sqrt(sigma) ),
                         n.chains = 4,
                         n.adapt = n.iter)
      
      d <- coda.samples(jags,
                        c('eta1', 'eta2', 'eta3', 'eta4', 'eta5', 'eta6'),
                        n.iter, thin=1)
    } else if((n.comp>=3 & n.comp <6 ) & (!steady_state)) {
      
      jags <- jags.model('AXI_I_multiple_dose.bug',
                         data = list('c' = tdm_data$conc,
                                     'amt' = dosing_events$amt, 
                                     'dosing_time' = dosing_events$time,
                                     'times'= tdm_data$time,
                                     "theta"=thetas,
                                     "omega"=sqrt(omegas),
                                     'sigma'=sqrt(sigma) ),
                         n.chains = 4,
                         n.adapt = n.iter)
      
      d <- coda.samples(jags,
                        c('eta1', 'eta2', 'eta3', 'eta4', 'eta5'),
                        n.iter, thin=1)
    } else if((n.comp>=3 & n.comp <6 ) & (steady_state)) {
      jags <- jags.model('AXI_I_ss.bug',
                         data = list('c' = tdm_data$conc,
                                     'amt' = dosing_events$amt[1], 
                                     'ii' = dosing_events$ii[1],
                                     'dosing_time' = dosing_events$time[1],
                                     'times'= tdm_data$time,
                                     "theta"=thetas,
                                     "omega"=sqrt(omegas),
                                     'sigma'=sqrt(sigma) ),
                         n.chains = 4,
                         n.adapt = n.iter)
      
      d <- coda.samples(jags,
                        c('eta1', 'eta2', 'eta3', 'eta4', 'eta5'),
                        n.iter, thin=1)
    }
    
    
    
    pk_profile <- (do_plot(d))
    
    
    mcmc_plots_1 <- mcmc_diagnosticplots(1, d, nburn=n.burn, omega=omegas, "red")
    mcmc_plots_2 <- mcmc_diagnosticplots(2, d, nburn=n.burn, omega=omegas, "orange")
    mcmc_plots_3 <- mcmc_diagnosticplots(3, d, nburn=n.burn, omega=omegas, "yellow")
    mcmc_plots_4 <- mcmc_diagnosticplots(4, d, nburn=n.burn, omega=omegas)
    
    ## Use the fourth chain for the simulation
    df <- as.data.frame(d[[4]])
    
    ## remove burnin iterations
    df <- df[-(1:n.burn),]
    
    
      result = list(mcmc_plots_1, 
                    mcmc_plots_2, 
                    mcmc_plots_3, 
                    mcmc_plots_4, 
                    pk_profile=pk_profile[[1]], 
                    c_at_tlast=pk_profile[[2]], 
                    ind_y_max=pk_profile[[3]], 
                    ind_y_min=pk_profile[[4]], 
                    mcmc_etas = df)
      
      return(result)
}


perform_mc_simulation <- function(n.mc, omegas, pk_mod, is_ss, thetas, app_data, t_from, t_to) {
  
  ## Todo: Outsource the mc simulation 
  if(pk_mod <=2){
    mc_eta2 <- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[2])))
    mc_eta3<- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[3])))
    mc_eta1<- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[1])))
    mc_eta4<- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[4])))
    
    all_etas <- data.frame(ETA1=mc_eta1, ETA2=mc_eta2, ETA3=mc_eta3, ETA4=mc_eta4)
  }
  if(pk_mod==2){
    ### generic 2 compartment models
    ### need two additional etas
    mc_eta5<- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[5])))
    mc_eta6<- (rnorm(n = n.mc, mean=0, sd=sqrt(omegas[6])))
    
    all_etas <- data.frame(ETA1=mc_eta1, ETA2=mc_eta2, ETA3=mc_eta3, ETA4=mc_eta4, ETA5=mc_eta5, ETA6=mc_eta6)
    
  } else if (pk_mod==3){
    
    ### Generate n.mc samples for every random effect
    ### Correlation between V and Cl and Q and V2 is simulated by 
    ### sampling from multivariate normal distributions
    cov_mat_1 <- matrix(c(axi_i_mod_fed$omegas[3], axi_i_mod_fed$omegas[6],
                          axi_i_mod_fed$omegas[6], axi_i_mod_fed$omegas[2]), nrow=2, ncol=2)
    
    cov_mat_2 <- matrix(c(axi_i_mod_fed$omegas[5], axi_i_mod_fed$omegas[7],
                          axi_i_mod_fed$omegas[7], axi_i_mod_fed$omegas[4]), nrow=2, ncol=2)
    
    means <- c(0,0)
    
    etas1 <- mvrnorm(n=n.mc, means, cov_mat_1)
    
    etas2 <- mvrnorm(n=n.mc, means, cov_mat_2)
    
    mc_eta2 <- etas1[,2]
    mc_eta3 <- etas1[,1]
    
    mc_eta1<- (rnorm(n = n.mc, mean=0, sd=sqrt(axi_i_mod_fed$omegas[3])))
    mc_eta4<- etas2[,2]
    mc_eta5<- etas2[,1]
    
    all_etas <- data.frame(ETA1=mc_eta1, ETA2=mc_eta2, ETA3=mc_eta3, ETA4=mc_eta4, ETA5=mc_eta5)
    
  } else if (pk_mod==4){
    
    ### Generate n.mc samples for every random effect
    ### Correlation between V and Cl and Q and V2 is simulated by 
    ### sampling from multivariate normal distributions
    cov_mat_1 <- matrix(c(axi_i_mod_fasted$omegas[3], axi_i_mod_fasted$omegas[6],
                          axi_i_mod_fasted$omegas[6], axi_i_mod_fasted$omegas[2]), nrow=2, ncol=2)
    
    cov_mat_2 <- matrix(c(axi_i_mod_fasted$omegas[5], axi_i_mod_fasted$omegas[7],
                          axi_i_mod_fasted$omegas[7], axi_i_mod_fasted$omegas[4]), nrow=2, ncol=2)
    
    means <- c(0,0)
    
    etas1 <- mvrnorm(n=n.mc, means, cov_mat_1)
    
    etas2 <- mvrnorm(n=n.mc, means, cov_mat_2)
    
    mc_eta2 <- etas1[,2]
    mc_eta3 <- etas1[,1]
    
    mc_eta1<- (rnorm(n = n.mc, mean=0, sd=sqrt(axi_i_mod_fed$omegas[3])))
    mc_eta4<- etas2[,2]
    mc_eta5<- etas2[,1]
    
    all_etas <- data.frame(ETA1=mc_eta1, ETA2=mc_eta2, ETA3=mc_eta3, ETA4=mc_eta4, ETA5=mc_eta5)
    
  } else if (pk_mod==5){
    
    ### Generate n.mc samples for every random effect
    ### Correlation between V and Cl and Q and V2 is simulated by 
    ### sampling from multivariate normal distributions
    cov_mat_1 <- matrix(c(axi_i_mod_fed_form_XLI$omegas[3], axi_i_mod_fed_form_XLI$omegas[6],
                          axi_i_mod_fed_form_XLI$omegas[6], axi_i_mod_fed_form_XLI$omegas[2]), nrow=2, ncol=2)
    
    cov_mat_2 <- matrix(c(axi_i_mod_fed_form_XLI$omegas[5], axi_i_mod_fed_form_XLI$omegas[7],
                          axi_i_mod_fed_form_XLI$omegas[7], axi_i_mod_fed_form_XLI$omegas[4]), nrow=2, ncol=2)
    
    means <- c(0,0)
    
    etas1 <- mvrnorm(n=n.mc, means, cov_mat_1)
    
    etas2 <- mvrnorm(n=n.mc, means, cov_mat_2)
    
    mc_eta2 <- etas1[,2]
    mc_eta3 <- etas1[,1]
    
    mc_eta1<- (rnorm(n = n.mc, mean=0, sd=sqrt(axi_i_mod_fed$omegas[3])))
    mc_eta4<- etas2[,2]
    mc_eta5<- etas2[,1]
    
    all_etas <- data.frame(ETA1=mc_eta1, ETA2=mc_eta2, ETA3=mc_eta3, ETA4=mc_eta4, ETA5=mc_eta5)
  }
  
  dat_mc <- NULL
  
  
  ## Do with progress
  withProgress(message = "Performing Monte Carlo simulation", max = n.mc, {
    for(i in 1:n.mc){
      ### Simulate population PK profiles for every monte carlo sample according to the PK model currently
      ### Selected
      if(pk_mod==1 & !is_ss) {  
        CP_mc <-  pk_1cmt_oral(theta = thetas, 
                               eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i]), 
                               dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                          amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                               times=seq(t_from, t_to, by=0.2))
      } else if (pk_mod==2 & !is_ss) {
        CP_mc <-  pk_2cmt_oral(theta = thetas, 
                               eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i],mc_eta5[i],mc_eta6[i]), 
                               dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                          amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                               times=seq(t_from, t_to, by=0.2))
        
      } else if (pk_mod==1 & is_ss) {
        CP_mc <-  pk_1cmt_oral_ss(theta = thetas, 
                                  eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i]), 
                                  dosing_events = data.frame(time=0,
                                                             amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                             ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                  times=seq(t_from, t_to, by=0.2))
      } else if (pk_mod==2 & is_ss) {
        CP_mc <-  pk_1cmt_oral_ss(theta = thetas, 
                                  eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i],mc_eta5[i],mc_eta6[i]), 
                                  dosing_events = data.frame(time=0,
                                                             amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                             ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                  times=seq(t_from, t_to, by=0.2))
      } else if (pk_mod==3 & !is_ss) {
        ## structural model of axitinib base model is generic oral 2 compartments
        CP_mc <-  pk_2cmt_oral(theta = axi_i_mod_fed$thetas, 
                               eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0,mc_eta4[i],mc_eta5[i]), 
                               dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                          amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                               times=seq(t_from, t_to, by=0.2))
        
      } else if (pk_mod==3 & is_ss) {
        ## structural model of axitinib base model is generic oral 2 compartments
        CP_mc <-  pk_2cmt_oral_ss(theta = axi_i_mod_fed$thetas, 
                                  eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0, mc_eta4[i],mc_eta5[i]), 
                                  dosing_events = data.frame(time=0,
                                                             amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                             ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                  times=seq(t_from, t_to, by=0.2))
        
      } else if (pk_mod==4 & !is_ss) {
        ## structural model of axitinib base model is generic oral 2 compartments
        CP_mc <-  pk_2cmt_oral(theta = axi_i_mod_fasted$thetas, 
                               eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0,mc_eta4[i],mc_eta5[i]), 
                               dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                          amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                               times=seq(t_from, t_to, by=0.2))
        
      } else if (pk_mod==4 & is_ss) {
        ## structural model of axitinib base model is generic oral 2 compartments
        CP_mc <-  pk_2cmt_oral_ss(theta = axi_i_mod_fasted$thetas, 
                                  eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0, mc_eta4[i],mc_eta5[i]), 
                                  dosing_events = data.frame(time=0,
                                                             amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                             ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                  times=seq(t_from, t_to, by=0.2))
        
      } else if (pk_mod==5 & !is_ss) {
        ## structural model of axitinib base model is generic oral 2 compartments
        CP_mc <-  pk_2cmt_oral(theta = axi_i_mod_fed_form_XLI$thetas, 
                               eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0,mc_eta4[i],mc_eta5[i]), 
                               dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                          amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                               times=seq(t_from, t_to, by=0.2))
        
      } else if (pk_mod==5 & is_ss) {
        ## structural model of axitinib base model is generic oral 2 compartments
        CP_mc <-  pk_2cmt_oral_ss(theta = axi_i_mod_fed_form_XLI$thetas, 
                                  eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0, mc_eta4[i],mc_eta5[i]), 
                                  dosing_events = data.frame(time=0,
                                                             amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                             ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                  times=seq(t_from, t_to, by=0.2))
        
      } 
      
      
      
      
      dat_mc <- cbind(dat_mc, CP_mc)
      incProgress(1)
    }
  })
  
  ## transpose data -> rows into columns
  dat_mc <- t(dat_mc)
  
  ## Get quantile from monte carlo result
  s <- apply(dat_mc,2,function(x) quantile(x,probs=c(0.05,0.5,0.95)))
  
  ## Data for population PK Plot
  plot_dat <- data.frame(TIME=seq(t_from, t_to, by=0.2),CP_min=s[1,],CP=s[2,],CP_max=s[3,], DELTA=(s[3,]-s[1,]))
  return(list(plot_dat, dat_mc, all_etas))
}
