


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
    
    mcmc_result = NULL,
    mc_result= NULL,
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
    
    ### In this very first part: Covariates need to change the THETA in global.R for coded model parameters
    ##
    if(input$choose_PK_mod >=3 & input$choose_PK_mod <=5){
      
      #Change typical value for Vc with WT according to pop PK model
      
      axi_i_mod_fed$thetas[2] <- axi_i_mod_fed$thetas[2] *(input$WT/75.0)^0.758
      axi_i_mod_fasted$thetas[2] <- axi_i_mod_fasted$thetas[2] *(input$WT/75.0)^0.758
      axi_i_mod_fed_form_XLI$thetas[2] <- axi_i_mod_fed_form_XLI$thetas[2] *(input$WT/75.0)^0.758
    }
    
    ## Decide which PK model was choosen and act accordingly
    if(input$choose_PK_mod==2){
      
      app_data$mcmc_result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                         thetas = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q),
                                         omegas = c(input$omega1, input$omega2, input$omega3, input$omega4, input$omega5, input$omega6),
                                         TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                         steady_state = input$choose_SS, n.comp=input$choose_PK_mod) 
                             
    } else if(input$choose_PK_mod==1) {
        app_data$mcmc_result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                           thetas = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q),
                                           omegas = c(input$omega1, input$omega2, input$omega3, input$omega4, input$omega5),
                                           TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                           steady_state = input$choose_SS, n.comp=input$choose_PK_mod)
    } else if(input$choose_PK_mod==3) {
      ### Dont get Thetas and Omega matrix from UI, use axi_i_mod_fed from global.R
      app_data$mcmc_result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                         thetas = axi_i_mod_fed$thetas,
                                         omegas = axi_i_mod_fed$omegas,
                                         TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                         steady_state = input$choose_SS, n.comp=input$choose_PK_mod) 
    } else if(input$choose_PK_mod==4) {
      ### Dont get Thetas and Omega matrix from UI, use axi_i_mod_fed from global.R
      app_data$mcmc_result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                         thetas = axi_i_mod_fasted$thetas,
                                         omegas = axi_i_mod_fasted$omegas,
                                         TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                         steady_state = input$choose_SS, n.comp=input$choose_PK_mod) 
    }
    else if(input$choose_PK_mod==5) {
      ### Dont get Thetas and Omega matrix from UI, use axi_i_mod_fed from global.R
      app_data$mcmc_result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                         thetas = axi_i_mod_fed_form_XLI$thetas,
                                         omegas = axi_i_mod_fed_form_XLI$omegas,
                                         TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                         steady_state = input$choose_SS, n.comp=input$choose_PK_mod) 
    }
    
    
    ## Perform mc simulation to get a prediction interval for the population PK curve
    app_data$mc_result <- perform_mc_simulation(input$n.mc, ## number of simulations
                                                c(input$omega1, input$omega2, input$omega3, input$omega4, input$omega5, input$omega6), ## omegas
                                                input$choose_PK_mod, ## PK model
                                                input$choose_SS, ## Steady state? True or False
                                                c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q), ## thetas
                                                app_data, ## App Data for Dosing / TDM Data
                                                input$TIME[1], input$TIME[2]) ## Time to simulate
    
    plot_dat <- app_data$mc_result [[1]] ## get Plot data ...
    dat_mc <- app_data$mc_result [[2]]  ### .. and raw results from the mc simulation to complete the plots
    
    ## Get lowest value (non-zero <- log-scale) and max value in PK plot to adjust y-axis
    pop_y_max <- max(plot_dat$CP_max)
    pop_y_min <- min(plot_dat$CP_min[plot_dat$CP_min >0])
    ind_y_max <- app_data$mcmc_result[[7]]
    ind_y_min <- app_data$mcmc_result[[8]]
    
     
    ## get tdm data from the table
    tdm_data <- data.frame(conc=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$conc)),
                           time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$time)))
    
    ## Check if TDM concentration is above upper prediction interval to readjust the plot limits on y-axis
    pop_y_max <- ifelse(max(tdm_data$conc) > pop_y_max, max(tdm_data$conc), pop_y_max)
    ind_y_max <- ifelse(max(tdm_data$conc) > ind_y_max, max(tdm_data$conc), ind_y_max)
    
    ## prepare individual boxplot
    ind_boxplot <- ggplot(data=data.frame(conc=app_data$mcmc_result[[6]], time="")) + geom_boxplot(aes(x=time, y=conc)) + theme_bw()  +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
          axis.ticks.y = element_blank())+ 
      ggtitle(paste("C last at ", input$TIME[2], " h"), "Individual") + xlab("") + ylim(c(0,ind_y_max*1.2))
    
    ## prepare individual PK plot
    ind_plot <- app_data$mcmc_result[[5]] + theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") +
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
    chart.Correlation(app_data$mcmc_result[[chain]]$chain_data, histogram=TRUE)

    
  })
  
  output$pop_cov_plot <- renderPlot({
    ## Show correlation matrix - the lazy way
    chart.Correlation(app_data$mc_result[[3]], histogram=TRUE)
  
  })
  
  output$par_dist <- renderPlot({
    dist_plots <- list()

    current_etas <- app_data$mc_result[[3]]
    
    mcmc_etas <- app_data$mcmc_result[[9]]
    
    if(input$choose_PK_mod==1){
        pop_pars <- data.frame(
          pop_ka=input$ka * exp(current_etas[,1]),
          pop_Vc=input$V * exp(current_etas[,2]),
          pop_Cl=input$Cl * exp(current_etas[,3]),
          pop_F_oral=input$F_oral * exp(current_etas[,4]))
        ind_pars <- data.frame(
          ind_ka=input$ka * exp(mcmc_etas[,1]),
          ind_Vc=input$V  * exp(mcmc_etas[,2]),
          ind_Cl=input$Cl * exp(mcmc_etas[,3]),
          ind_F_oral=input$F_oral  * exp(mcmc_etas[,4]))
        
        
        dist_plots[[1]] <- ggplot() +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_ka, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
          geom_density(data=ind_pars, aes(x=ind_ka, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1) +
          xlab("Absorption rate constant (ka) [1/h]") + ylab("Frequency")
        dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_Vc, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
          geom_density(data=ind_pars,aes(x=ind_Vc, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
          xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")
        dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_Cl, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
          geom_density(data=ind_pars,aes(x=ind_Cl, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
          xlab("Clearance (Cl) [L/h]") + ylab("Frequency")
        dist_plots[[4]] <- ggplot(  ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_F_oral, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
          geom_density(data=ind_pars,aes(x=ind_F_oral, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
          xlab("Systemically available fraction") + ylab("Frequency")
        
        grid.arrange(dist_plots[[1]],dist_plots[[2]],dist_plots[[3]],dist_plots[[4]], nrow=2, ncol=2)
      
    }
    else if (input$choose_PK_mod == 2){
      pop_pars <- data.frame(
        pop_ka=input$ka * exp(current_etas[,1]),
        pop_Vc=input$V * exp(current_etas[,2]),
        pop_Cl=input$Cl * exp(current_etas[,3]),
        pop_F_oral=input$F_oral * exp(current_etas[,4]),
        pop_Vp=input$V2 * exp(current_etas[,5]),
        pop_Q=input$Q * exp(current_etas[,6]))
      ind_pars <- data.frame(
        ind_ka=input$ka * exp(mcmc_etas[,1]),
        ind_Vc=input$V  * exp(mcmc_etas[,2]),
        ind_Cl=input$Cl * exp(mcmc_etas[,3]),
        ind_F_oral=input$F_oral  * exp(mcmc_etas[,4]),
        ind_Vp=input$V2 * exp(mcmc_etas[,5]),
        ind_Q=input$Q* exp(mcmc_etas[,5]))
      
      
      dist_plots[[1]] <- ggplot() +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_ka, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
        geom_density(data=ind_pars, aes(x=ind_ka, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1) +
        xlab("Absorption rate constant (ka) [1/h]") + ylab("Frequency")
      dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Vc, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
        geom_density(data=ind_pars,aes(x=ind_Vc, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")
      dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Cl, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Cl, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Clearance (Cl) [L/h]") + ylab("Frequency")
      dist_plots[[4]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Vp, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Vp, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Volume of peripheral compartment (Vp) [L]") + ylab("Frequency")
      dist_plots[[5]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Q, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Q, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Intercompartmental Clearance (Q) [L/h]") + ylab("Frequency")
      dist_plots[[6]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_F_oral, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_F_oral, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Systemically available fraction") + ylab("Frequency")
      
      grid.arrange(dist_plots[[1]],dist_plots[[2]],dist_plots[[3]],dist_plots[[4]],dist_plots[[5]],dist_plots[[6]], nrow=3, ncol=2)
    }
    else if(input$choose_PK_mod==3) {
    
        pop_pars <- data.frame(
                          pop_ka=axi_i_mod_fed$thetas[1] * exp(current_etas[,1]),
                          pop_Vc=axi_i_mod_fed$thetas[2] * exp(current_etas[,2]),
                          pop_Cl=axi_i_mod_fed$thetas[3] * exp(current_etas[,3]),
                          pop_Vp=axi_i_mod_fed$thetas[6] * exp(current_etas[,4]),
                          pop_Q=axi_i_mod_fed$thetas[7] * exp(current_etas[,5]))
        ind_pars <- data.frame(
                          ind_ka=axi_i_mod_fed$thetas[1] * exp(mcmc_etas[,1]),
                          ind_Vc=axi_i_mod_fed$thetas[2] * exp(mcmc_etas[,2]),
                          ind_Cl=axi_i_mod_fed$thetas[3] * exp(mcmc_etas[,3]),
                          ind_Vp=axi_i_mod_fed$thetas[6] * exp(mcmc_etas[,4]),
                          ind_Q=axi_i_mod_fed$thetas[7] * exp(mcmc_etas[,5]))
    
        
        dist_plots[[1]] <- ggplot() +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_ka, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
          geom_density(data=ind_pars, aes(x=ind_ka, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1) +
          xlab("Absorption rate constant (ka) [1/h]") + ylab("Frequency")
        dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_Vc, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
          geom_density(data=ind_pars,aes(x=ind_Vc, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
          xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")
        dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_Cl, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
          geom_density(data=ind_pars,aes(x=ind_Cl, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
          xlab("Clearance (Cl) [L/h]") + ylab("Frequency")
        dist_plots[[4]] <- ggplot(  ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_Vp, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
          geom_density(data=ind_pars,aes(x=ind_Vp, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
          xlab("Volume of peripheral compartment (Vp) [L]") + ylab("Frequency")
        dist_plots[[5]] <- ggplot(  ) +theme_bw()+ 
          geom_density(data=pop_pars, aes(x=pop_Q, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
          geom_density(data=ind_pars,aes(x=ind_Q, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
          xlab("Intercompartmental Clearance (Q) [L/h]") + ylab("Frequency")
        
        grid.arrange(dist_plots[[1]],dist_plots[[2]],dist_plots[[3]],dist_plots[[4]],dist_plots[[5]], nrow=3, ncol=2)
    } else if (input$choose_PK_mod==4){
      pop_pars <- data.frame(
        pop_ka=axi_i_mod_fasted$thetas[1] * exp(current_etas[,1]),
        pop_Vc=axi_i_mod_fasted$thetas[2] * exp(current_etas[,2]),
        pop_Cl=axi_i_mod_fasted$thetas[3] * exp(current_etas[,3]),
        pop_Vp=axi_i_mod_fasted$thetas[6] * exp(current_etas[,4]),
        pop_Q=axi_i_mod_fasted$thetas[7] * exp(current_etas[,5]))
      ind_pars <- data.frame(
        ind_ka=axi_i_mod_fasted$thetas[1] * exp(mcmc_etas[,1]),
        ind_Vc=axi_i_mod_fasted$thetas[2] * exp(mcmc_etas[,2]),
        ind_Cl=axi_i_mod_fasted$thetas[3] * exp(mcmc_etas[,3]),
        ind_Vp=axi_i_mod_fasted$thetas[6] * exp(mcmc_etas[,4]),
        ind_Q=axi_i_mod_fasted$thetas[7] * exp(mcmc_etas[,5]))
      
      
      dist_plots[[1]] <- ggplot() +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_ka, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
        geom_density(data=ind_pars, aes(x=ind_ka, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1) +
        xlab("Absorption rate constant (ka) [1/h]") + ylab("Frequency")
      dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Vc, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
        geom_density(data=ind_pars,aes(x=ind_Vc, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")
      dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Cl, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Cl, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Clearance (Cl) [L/h]") + ylab("Frequency")
      dist_plots[[4]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Vp, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Vp, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Volume of peripheral compartment (Vp) [L]") + ylab("Frequency")
      dist_plots[[5]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Q, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Q, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Intercompartmental Clearance (Q) [L/h]") + ylab("Frequency")
      
      grid.arrange(dist_plots[[1]],dist_plots[[2]],dist_plots[[3]],dist_plots[[4]],dist_plots[[5]], nrow=3, ncol=2)
    } else if (input$choose_PK_mod==5){
      pop_pars <- data.frame(
        pop_ka=axi_i_mod_fed_form_XLI$thetas[1] * exp(current_etas[,1]),
        pop_Vc=axi_i_mod_fed_form_XLI$thetas[2] * exp(current_etas[,2]),
        pop_Cl=axi_i_mod_fed_form_XLI$thetas[3] * exp(current_etas[,3]),
        pop_Vp=axi_i_mod_fed_form_XLI$thetas[6] * exp(current_etas[,4]),
        pop_Q=axi_i_mod_fed_form_XLI$thetas[7] * exp(current_etas[,5]))
      ind_pars <- data.frame(
        ind_ka=axi_i_mod_fed_form_XLI$thetas[1] * exp(mcmc_etas[,1]),
        ind_Vc=axi_i_mod_fed_form_XLI$thetas[2] * exp(mcmc_etas[,2]),
        ind_Cl=axi_i_mod_fed_form_XLI$thetas[3] * exp(mcmc_etas[,3]),
        ind_Vp=axi_i_mod_fed_form_XLI$thetas[6] * exp(mcmc_etas[,4]),
        ind_Q=axi_i_mod_fed_form_XLI$thetas[7] * exp(mcmc_etas[,5]))
      
      
      dist_plots[[1]] <- ggplot() +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_ka, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
        geom_density(data=ind_pars, aes(x=ind_ka, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1) +
        xlab("Absorption rate constant (ka) [1/h]") + ylab("Frequency")
      dist_plots[[2]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Vc, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1) +
        geom_density(data=ind_pars,aes(x=ind_Vc, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Volume of central compartment (Vc) [L]") + ylab("Frequency")
      dist_plots[[3]] <- ggplot( ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Cl, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Cl, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Clearance (Cl) [L/h]") + ylab("Frequency")
      dist_plots[[4]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Vp, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Vp, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Volume of peripheral compartment (Vp) [L]") + ylab("Frequency")
      dist_plots[[5]] <- ggplot(  ) +theme_bw()+ 
        geom_density(data=pop_pars, aes(x=pop_Q, y=..density..), colour="blue", size=.5, fill="blue",alpha=0.25, linetype=1)+
        geom_density(data=ind_pars,aes(x=ind_Q, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)+
        xlab("Intercompartmental Clearance (Q) [L/h]") + ylab("Frequency")
      
      grid.arrange(dist_plots[[1]],dist_plots[[2]],dist_plots[[3]],dist_plots[[4]],dist_plots[[5]], nrow=3, ncol=2)
    }
  })
  
  output$traceplot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    ## Has to be generalized 
    
    
    ## 
    if(input$choose_PK_mod==1) {
        gridExtra::grid.arrange(app_data$mcmc_result[[chain]]$p_iter_ETA1, app_data$mcmc_result[[chain]]$p_dens_ETA1,    
                                app_data$mcmc_result[[chain]]$p_iter_ETA2, app_data$mcmc_result[[chain]]$p_dens_ETA2, 
                                app_data$mcmc_result[[chain]]$p_iter_ETA3, app_data$mcmc_result[[chain]]$p_dens_ETA3, 
                                app_data$mcmc_result[[chain]]$p_iter_ETA4, app_data$mcmc_result[[chain]]$p_dens_ETA4, nrow=4, ncol=2, widths=c(3,1))  
    } else if(input$choose_PK_mod==2) {
        gridExtra::grid.arrange(app_data$mcmc_result[[chain]]$p_iter_ETA1, app_data$mcmc_result[[chain]]$p_dens_ETA1,    
                                app_data$mcmc_result[[chain]]$p_iter_ETA2, app_data$mcmc_result[[chain]]$p_dens_ETA2, 
                                app_data$mcmc_result[[chain]]$p_iter_ETA3, app_data$mcmc_result[[chain]]$p_dens_ETA3, 
                                app_data$mcmc_result[[chain]]$p_iter_ETA4, app_data$mcmc_result[[chain]]$p_dens_ETA4,
                                app_data$mcmc_result[[chain]]$p_iter_ETA5, app_data$mcmc_result[[chain]]$p_dens_ETA5,
                                app_data$mcmc_result[[chain]]$p_iter_ETA6, app_data$mcmc_result[[chain]]$p_dens_ETA6, nrow=6, ncol=2, widths=c(3,1)) 
    } else if(input$choose_PK_mod>=3 & input$choose_PK_mod <=5) {
      gridExtra::grid.arrange(app_data$mcmc_result[[chain]]$p_iter_ETA1, app_data$mcmc_result[[chain]]$p_dens_ETA1,    
                              app_data$mcmc_result[[chain]]$p_iter_ETA2, app_data$mcmc_result[[chain]]$p_dens_ETA2, 
                              app_data$mcmc_result[[chain]]$p_iter_ETA3, app_data$mcmc_result[[chain]]$p_dens_ETA3, 
                              app_data$mcmc_result[[chain]]$p_iter_ETA4, app_data$mcmc_result[[chain]]$p_dens_ETA4,
                              app_data$mcmc_result[[chain]]$p_iter_ETA5, app_data$mcmc_result[[chain]]$p_dens_ETA5,nrow=5, ncol=2, widths=c(3,1)) 
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
        } else if ( (input$choose_PK_mod>=3 & input$choose_PK_mod<=5) & input$choose_SS) {
          includeText("AXI_I_ss.bug")
        } else if ((input$choose_PK_mod>=3 & input$choose_PK_mod<=5) & !input$choose_SS) {
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
      paste("Demo using the Axitinib Final model in <b>fed</b> state from Garrett et al. British J Clin Pharmacol, DOI: 10.1111/bcp.12206")
    } else if(input$choose_PK_mod==4) {
      paste("Demo using the Axitinib Final model in <b>fasted</b> state from Garrett et al. British J Clin Pharmacol, DOI: 10.1111/bcp.12206")
    } else if(input$choose_PK_mod==5) {
      paste("Demo using the Axitinib Final model in <b>fed</b> state and the <b>final formulation</b>  from Garrett et al. British J Clin Pharmacol, DOI: 10.1111/bcp.12206")
    }
  })
  
  observeEvent(input$submit,
               ## Submit changes
               app_data$pk_plots <- updatePKPlot()
               )
  
  
  
})


