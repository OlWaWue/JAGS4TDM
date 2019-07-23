library('rjags')
library('ggplot2')

theo <- Theoph

theo$Dose <- theo$Dose

pk.model1 <- function(psi, t, Dose) {
  D <- Dose
  ka <- psi[1]
  V <- psi[2]
  ke <- psi[3]
  f <- D*ka/V/(ka-ke)*(exp(-ke*t)-exp(-ka*t))
  return(f)
}

pkm1 <- nls(conc ~ pk.model1(psi, Time, Dose), start=list(psi=c(ka=1, V=10, ke=0.1)), data=subset(theo, Subject==1))

coef(pkm1)



la <- lapply(1:12, function(iSubj){
  theo.i <- subset(theo, Subject == iSubj)
  
  print(theo.i)
  
  jags <- jags.model('pk.bug',
                     data = list('c' = theo.i$conc,
                                 'D' = theo.i$Dose[1], 
                                 'ts'= theo.i$Time,
                                 'n' = nrow(theo.i),
                                 "theta"=c(0.054,37,1.78),
                                 "omega"=c(2.51,4.54,5.53)),
                     n.chains = 4,
                     n.adapt = 5000)
  d <- coda.samples(jags,
                    c('eta1', 'eta2', 'eta3'),
                    5000, thin=1)
})


plot(do_plot(2, la))


do_plot <- function(patiend_id =1, la){
    

    
    
    
    dose_i <- subset(theo, Subject == patiend_id)$Dose[1]

    
    df <- as.data.frame(la[[patiend_id]][[4]])
    
    
    
    mcmc_se <- list()
    
    
     
      for (i in 1:nrow(df)) {
        mcmc_se[[i]] <- simulate_profile(D=dose_i, eta1 = df$eta1[i], eta2 = df$eta2[i], eta3 = df$eta3[i]) 
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
    
    
    s <- apply(df_temp,2,function(x) quantile(x,probs=c(0.1/2, 1-0.1/2)))
    
    pk_data <- data.frame(time=seq(0,24,0.1),s1=s[1,],s2=s[2,], avg=simulate_profile(D=dose_i), max=simulate_profile(D=dose_i, eta1 = max_eta1, 
                                                                                                                  eta2 = max_eta2, 
                                                                                                                  eta3 = max_eta3))
    
    
    p <- ggplot(pk_data) + geom_ribbon(aes(ymin=s1, ymax=s2, x=time), fill="red", alpha=0.5) + geom_line(aes(y=avg, x=time), linetype=2) +
          geom_line(aes(y=max, x=time))+
          geom_point(data = subset(theo, Subject == patiend_id), aes(x=Time, y=conc))
    
    return(p)
}


simulate_profile <- function(D, time=seq(0,24,0.1), eta1=0, eta2=0, eta3=0){

    v_F = 37*exp(eta2)
    ka = 1.78*exp(eta3)
    k = 0.054*exp(eta1)
    
    sim <-  (D)/v_F*(ka/(ka-k))*(exp(-1*(k)*time)-exp(-1*ka*time))
    return(sim)
}
