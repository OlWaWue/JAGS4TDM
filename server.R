


shinyServer(function(input, output, session) {
  
  ## Load an Excel file with data
  observeEvent(input$loadPT,{
    inFile <- input$loadPT
    
    if (is.null(inFile))
      return(NULL)
    
    tryCatch({
      app_data$data_set <- read_xlsx(inFile$datapath)
    },
    error = function(e){
      showModal(modalDialog(
        title = label_error,
        HTML(paste("File not recognized!<br>Details:<br><br>",e)),
        easyClose = TRUE,
        footer = NULL
      ))
    })

  })
  
  app_data <- reactiveValues(
    
    ## AppData used in simulation
    
    result = NULL,
    data_set = data.frame(time=c(0,4,6,12,24,30,36,48),
                          amt=c(5,".",".",5,5,".",5,"."),
                          conc=c(".", 0.04, 0.036, ".",".", 0.030,".", 0.0132),
                          evid=c(1, 0, 0, 1,1, 0, 1, 0),
                          ii=c("12",".",".","12","12",".","12",".")),
    
    demo_loaded = FALSE, # Flag shows wether demo simulation has been loaded
    pk_plots = NULL
  )
  
  ## This function takes care of all the simulations
  ## MCMC is outsourced to global.R
  updatePKPlot <- function(){
    

    if(input$choose_PK_mod==2){
      
      app_data$result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                         thetas = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q),
                                         omegas = c(input$omega1, input$omega2, input$omega3, input$omega4, input$omega5, input$omega6),
                                         TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                         steady_state = input$choose_SS, n.comp=input$choose_PK_mod) 
                             
    } else if(input$choose_PK_mod==1) {
        app_data$result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                           thetas = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q),
                                           omegas = c(input$omega1, input$omega2, input$omega3, input$omega4, input$omega5),
                                           TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                           steady_state = input$choose_SS, n.comp=input$choose_PK_mod)
    } else if(input$choose_PK_mod==3) {
      ### Dont get Thetas and Omega matrix from UI, use axi_i_mod_pars
      app_data$result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                         thetas = axi_i_mod_pars$thetas,
                                         omegas = axi_i_mod_pars$omegas,
                                         TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                         steady_state = input$choose_SS, n.comp=input$choose_PK_mod) 
    }
    
    
    ## Todo: Outsource the mc simulation 
    mc_eta2 <- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega2)))
    mc_eta3<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega3)))
    mc_eta1<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega1)))
    mc_eta4<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega4)))
    
    if(input$choose_PK_mod==2){
      ### generic 2 compartment models
      ### need two additional etas
      mc_eta5<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega5)))
      mc_eta6<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega6)))
    } else if (input$choose_PK_mod==3){
      
      ### Generate n.mc samples for every random effect
      ### Correlation between V and Cl and Q and V2 is simulated by 
      ### sampling from multivariate normal distributions
      cov_mat_1 <- matrix(c(axi_i_mod_pars$omegas[3], axi_i_mod_pars$omegas[6],
               axi_i_mod_pars$omegas[6], axi_i_mod_pars$omegas[2]), nrow=2, ncol=2)
      
      cov_mat_2 <- matrix(c(axi_i_mod_pars$omegas[5], axi_i_mod_pars$omegas[7],
               axi_i_mod_pars$omegas[7], axi_i_mod_pars$omegas[4]), nrow=2, ncol=2)
      
      means <- c(0,0)
      
      etas1 <- mvrnorm(n=input$n.mc, means, cov_mat_1)

      etas2 <- mvrnorm(n=input$n.mc, means, cov_mat_2)
      
      mc_eta2 <- etas1[,2]
      mc_eta3 <- etas1[,1]
      
      mc_eta1<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(axi_i_mod_pars$omegas[3])))
      mc_eta4<- etas1[,2]
      mc_eta5<- etas1[,1]
    }
    
    dat_mc <- NULL
    
    
    ## Do with progress
    withProgress(message = "Performing Monte Carlo simulation", max = input$n.mc, {
      for(i in 1:input$n.mc){
        ### Simulate population PK profiles for every monte carlo sample according to the PK model currently
        ### Selected
        if(input$choose_PK_mod==1 & !input$choose_SS) {  
            CP_mc <-  pk_1cmt_oral(theta = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag), 
                                   eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i]), 
                                   dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                              amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                                   times=seq(input$TIME[1], input$TIME[2], by=0.2))
        } else if (input$choose_PK_mod==2 & !input$choose_SS) {
          CP_mc <-  pk_2cmt_oral(theta = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q), 
                                 eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i],mc_eta5[i],mc_eta6[i]), 
                                 dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                            amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                                 times=seq(input$TIME[1], input$TIME[2], by=0.2))
          
        } else if (input$choose_PK_mod==1 & input$choose_SS) {
          CP_mc <-  pk_1cmt_oral_ss(theta = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag), 
                                 eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i]), 
                                 dosing_events = data.frame(time=0,
                                                            amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                            ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                 times=seq(input$TIME[1], input$TIME[2], by=0.2))
        } else if (input$choose_PK_mod==2 & input$choose_SS) {
          CP_mc <-  pk_1cmt_oral_ss(theta = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q), 
                                    eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],mc_eta4[i],mc_eta5[i],mc_eta6[i]), 
                                    dosing_events = data.frame(time=0,
                                                               amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                               ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                    times=seq(input$TIME[1], input$TIME[2], by=0.2))
        } else if (input$choose_PK_mod==3 & !input$choose_SS) {
          ## structural model of axitinib base model is generic oral 2 compartments
          CP_mc <-  pk_2cmt_oral(theta = axi_i_mod_pars$thetas, 
                                 eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0,mc_eta4[i],mc_eta5[i]), 
                                 dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                            amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                                 times=seq(input$TIME[1], input$TIME[2], by=0.2))
          
        } else if (input$choose_PK_mod==3 & input$choose_SS) {
          ## structural model of axitinib base model is generic oral 2 compartments
          CP_mc <-  pk_2cmt_oral_ss(theta = axi_i_mod_pars$thetas, 
                                 eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i],0, mc_eta4[i],mc_eta5[i]), 
                                 dosing_events = data.frame(time=0,
                                                            amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt[1])),
                                                            ii=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$ii[1])), evid=1), 
                                 times=seq(input$TIME[1], input$TIME[2], by=0.2))
          
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
    plot_dat <- data.frame(TIME=seq(input$TIME[1], input$TIME[2], by=0.2),CP_min=s[1,],CP=s[2,],CP_max=s[3,], DELTA=(s[3,]-s[1,]))

    ## Get lowest value (non-zero <- log-scale) and max value in PK plot to adjust y-axis
    pop_y_max <- max(plot_dat$CP_max)
    pop_y_min <- min(plot_dat$CP_min[plot_dat$CP_min >0])
    ind_y_max <- app_data$result[[7]]
    ind_y_min <- app_data$result[[8]]
    
     
    ## get tdm data from the table
    tdm_data <- data.frame(conc=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$conc)),
                           time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$time)))
    
    ## Check if TDM concentration is above upper prediction interval to readjust the plot limits on y-axis
    pop_y_max <- ifelse(max(tdm_data$conc) > pop_y_max, max(tdm_data$conc), pop_y_max)
    ind_y_max <- ifelse(max(tdm_data$conc) > ind_y_max, max(tdm_data$conc), ind_y_max)
    
    ## prepare individual boxplot
    ind_boxplot <- ggplot(data=data.frame(conc=app_data$result[[6]], time="")) + geom_boxplot(aes(x=time, y=conc)) + theme_bw()  +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
          axis.ticks.y = element_blank())+ 
      ggtitle(paste("C last at ", input$TIME[2], " h"), "Individual") + xlab("") + ylim(c(0,ind_y_max*1.2))
    
    ## prepare individual PK plot
    ind_plot <- app_data$result[[5]] + theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") +
      ggtitle("MCMC Result including data (posterior)", "80/85/90/95% PI") + 
      geom_line(data=plot_dat, aes(x=TIME, y=CP), colour="blue", linetype=2) + ylim(c(0,ind_y_max*1.2))

    ## prepare population boxplot
    pop_boxplot <- ggplot(data=data.frame(conc=dat_mc[,ncol(dat_mc)], time="")) + geom_boxplot(aes(x=time, y=conc)) + theme_bw()  +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank())+ ggtitle(paste("C last at ", input$TIME[2], " h") , "Population") + 
      xlab("") + ylim(c(0,pop_y_max*1.2))
    
    ## Prepare population PK plot
    pop_plot <- ggplot(data=plot_dat)  + geom_line(aes(x=TIME, y=CP), colour="blue") +
      geom_ribbon(aes(x=TIME, ymax=CP_max, ymin=CP_min), alpha=0.15, fill="blue") +
      theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") + ggtitle("MC Result without data (prior)", "95% PI") +
      geom_point(data=tdm_data, aes(x=time, y=conc)) + ylim(c(0,pop_y_max*1.2))
    
    
    ## Convert the y-axis to log-scale
    if(input$log_y){
      ind_plot <- ind_plot + scale_y_continuous(trans = "log10", limits=c(ind_y_min,ind_y_max*1.2))
      pop_plot <- pop_plot + scale_y_continuous(trans = "log10", limits=c(pop_y_min,ind_y_max*1.2))
      ind_boxplot <- ind_boxplot + scale_y_continuous(trans = "log10", limits=c(ind_y_min,pop_y_max*1.2))
      pop_boxplot <- pop_boxplot + scale_y_continuous(trans = "log10", limits=c(pop_y_min,pop_y_max*1.2))
    }
    
    plots <- list(ind_plot, pop_plot, ind_boxplot, pop_boxplot)
  }
  
  output$pkPlot <- renderPlot({
    ## Workaround, that allows demo plot without hitting the submit button
    if(!app_data$demo_loaded) {
      app_data$pk_plots <- updatePKPlot()
      app_data$demo_loaded <- TRUE
    }
    ## Combine PK Profile with Boxplot
    grid.arrange(app_data$pk_plots[[2]], app_data$pk_plots[[4]], 
                 app_data$pk_plots[[1]], app_data$pk_plots[[3]], nrow=2, ncol=2,widths=c(4,1))
  })

  output$data_set <- renderTable({
    display_data <-app_data$data_set
    
    ## In Steady state, only the first dosing event is used!
    if(input$choose_SS){
      dosing_data <- display_data[display_data$evid==1,]
      tdm_data <- display_data[display_data$evid==0,]
      
      display_data <- rbind(dosing_data[1,], tdm_data)
   
    }
    
    display_data
    
  })
  
  output$cov_plot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    ## Show correlation matrix - the lazy way
    chart.Correlation(app_data$result[[chain]]$chain_data, histogram=TRUE)

    
  })
  
  output$traceplot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    ## Has to be generalized 
    
    
    ## 
    if(input$choose_PK_mod==1) {
        gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,    
                                app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2, 
                                app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, 
                                app_data$result[[chain]]$p_iter_ETA4, app_data$result[[chain]]$p_dens_ETA4, nrow=4, ncol=2, widths=c(3,1))  
    } else if(input$choose_PK_mod==2) {
        gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,    
                                app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2, 
                                app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, 
                                app_data$result[[chain]]$p_iter_ETA4, app_data$result[[chain]]$p_dens_ETA4,
                                app_data$result[[chain]]$p_iter_ETA5, app_data$result[[chain]]$p_dens_ETA5,
                                app_data$result[[chain]]$p_iter_ETA6, app_data$result[[chain]]$p_dens_ETA6, nrow=6, ncol=2, widths=c(3,1)) 
    } else if(input$choose_PK_mod==3) {
      gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,    
                              app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2, 
                              app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, 
                              app_data$result[[chain]]$p_iter_ETA4, app_data$result[[chain]]$p_dens_ETA4,
                              app_data$result[[chain]]$p_iter_ETA5, app_data$result[[chain]]$p_dens_ETA5,nrow=5, ncol=2, widths=c(3,1)) 
    }
    
    
  })
  
  
  output$modelfile <- renderText({

    ## Show the correct JAGS model file
    
      if(input$choose_PK_mod==1 & !input$choose_SS) { 
          includeText("1cmt_multiple_dose.bug")
        } else if (input$choose_PK_mod==2 & !input$choose_SS) {
          includeText("2cmt_multiple_dose.bug")
        } else if (input$choose_PK_mod==1 & input$choose_SS) {
          includeText("1cmt_ss.bug")
        } else if (input$choose_PK_mod==2 & input$choose_SS) {
          includeText("2cmt_ss.bug")
        } else if (input$choose_PK_mod==3 & input$choose_SS) {
          includeText("AXI_I_ss.bug")
        } else if (input$choose_PK_mod==3 & !input$choose_SS) {
          includeText("AXI_I_multiple_dose.bug")
        }
    
  })
  
  output$info <- renderText({
    
    ## Show a short info about the currently selected PK Model
    if (input$choose_PK_mod==1) {
      paste("Generic one compartment model with first order absorption and linear elimination")
    } else if(input$choose_PK_mod==2){
      paste("Generic two compartment model with first order absorption and linear elimination")
    } else if(input$choose_PK_mod==3) {
      paste("Demo using the Axitinib Base model from Garrett et al. British J Clin Pharmacol, DOI: 10.1111/bcp.12206")
    }
  })
  
  observeEvent(input$submit,
               ## Submit changes
               app_data$pk_plots <- updatePKPlot()
               )
  
  
  
})


