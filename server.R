


shinyServer(function(input, output, session) {
   
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
    result = NULL,
    data_set = data.frame(time=c(0,4,6,12,24,30,36,48),
                          amt=c(5,".",".",5,5,".",5,"."),
                          conc=c(".", 0.04, 0.036, ".",".", 0.030,".", 0.0132),
                          evid=c(1, 0, 0, 1,1, 0, 1, 0),
                          ii=c("12",".",".","12","12",".","12",".")),
    
    demo_loaded = FALSE,
    pk_plots = NULL
  )
  
  updatePKPlot <- function(){
    
    ## Remove content of if-statement and the complete if then else check once Steady state and
    ## 2 cmt model are implemented
    if(input$choose_PK_mod==2){
      
      app_data$result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                         thetas = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q),
                                         omegas = c(input$omega1, input$omega2, input$omega3, input$omega4, input$omega5, input$omega6, input$cov1, input$cov2),
                                         TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                         steady_state = input$choose_SS, n.comp=input$choose_PK_mod) 
                             
    } else {
        app_data$result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                           thetas = c(input$ka, input$V, input$Cl, input$F_oral, input$tlag, input$V2, input$Q),
                                           omegas = c(input$omega1, input$omega2, input$omega3, input$omega4, input$omega5),
                                           TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma, 
                                           steady_state = input$choose_SS, n.comp=input$choose_PK_mod)
    }
    
    
    ## Todo: Outsource the mc simulation 
    mc_eta2 <- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega2)))
    mc_eta3<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega3)))
    mc_eta1<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega1)))
    mc_eta4<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega4)))
    
    if(input$choose_PK_mod==2){
      mc_eta5<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega5)))
      mc_eta6<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega6)))
    }
    
    dat_mc <- NULL
    
    
    ## Do with progress
    withProgress(message = "Performing Monte Carlo simulation", max = input$n.mc, {
      for(i in 1:input$n.mc){
        
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
        }
        
        
        dat_mc <- cbind(dat_mc, CP_mc)
        incProgress(1)
      }
    })
    
    dat_mc <- t(dat_mc)
    
    s <- apply(dat_mc,2,function(x) quantile(x,probs=c(0.1,0.5,0.9)))
    
    plot_dat <- data.frame(TIME=seq(input$TIME[1], input$TIME[2], by=0.2),CP_min=s[1,],CP=s[2,],CP_max=s[3,], DELTA=(s[3,]-s[1,]))

    pop_y_max <- max(plot_dat$CP_max)
    pop_y_min <- min(plot_dat$CP_min[plot_dat$CP_min >0])
    ind_y_max <- app_data$result[[7]]
    ind_y_min <- app_data$result[[8]]
    
     
    
    tdm_data <- data.frame(conc=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$conc)),
                           time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$time)))
    
    pop_y_max <- ifelse(max(tdm_data$conc) > pop_y_max, max(tdm_data$conc), pop_y_max)
    ind_y_max <- ifelse(max(tdm_data$conc) > ind_y_max, max(tdm_data$conc), ind_y_max)
    
    ind_boxplot <- ggplot(data=data.frame(conc=app_data$result[[6]], time="")) + geom_boxplot(aes(x=time, y=conc)) + theme_bw()  +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
          axis.ticks.y = element_blank())+ 
      ggtitle(paste("C last at ", input$TIME[2], " h"), "Individual") + xlab("") + ylim(c(0,ind_y_max*1.2))
    
    
    ind_plot <- app_data$result[[5]] + theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") +
      ggtitle("MCMC Result including data (posterior)", "80/85/90/95% PI") + 
      geom_line(data=plot_dat, aes(x=TIME, y=CP), colour="blue", linetype=2) + ylim(c(0,ind_y_max*1.2))

    pop_boxplot <- ggplot(data=data.frame(conc=dat_mc[,ncol(dat_mc)], time="")) + geom_boxplot(aes(x=time, y=conc)) + theme_bw()  +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank())+ ggtitle(paste("C last at ", input$TIME[2], " h") , "Population") + 
      xlab("") + ylim(c(0,pop_y_max*1.2))
    
    pop_plot <- ggplot(data=plot_dat)  + geom_line(aes(x=TIME, y=CP), colour="blue") +
      geom_ribbon(aes(x=TIME, ymax=CP_max, ymin=CP_min), alpha=0.15, fill="blue") +
      theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") + ggtitle("MC Result without data (prior)", "90% PI") +
      geom_point(data=tdm_data, aes(x=time, y=conc)) + ylim(c(0,pop_y_max*1.2))
    
    if(input$log_y){
      ind_plot <- ind_plot + scale_y_continuous(trans = "log10", limits=c(ind_y_min,ind_y_max*1.2))
      pop_plot <- pop_plot + scale_y_continuous(trans = "log10", limits=c(pop_y_min,ind_y_max*1.2))
      ind_boxplot <- ind_boxplot + scale_y_continuous(trans = "log10", limits=c(ind_y_min,pop_y_max*1.2))
      pop_boxplot <- pop_boxplot + scale_y_continuous(trans = "log10", limits=c(pop_y_min,pop_y_max*1.2))
    }
    
    plots <- list(ind_plot, pop_plot, ind_boxplot, pop_boxplot)
  }
  
  output$pkPlot <- renderPlot({
    if(!app_data$demo_loaded) {
      app_data$pk_plots <- updatePKPlot()
      app_data$demo_loaded <- TRUE
    }
    grid.arrange(app_data$pk_plots[[2]], app_data$pk_plots[[4]], 
                 app_data$pk_plots[[1]], app_data$pk_plots[[3]], nrow=2, ncol=2,widths=c(4,1))
  })

  output$data_set <- renderTable({
    display_data <-app_data$data_set
    
    if(input$choose_SS){
      dosing_data <- display_data[display_data$evid==1,]
      tdm_data <- display_data[display_data$evid==0,]
      
      display_data <- rbind(dosing_data[1,], tdm_data)
   
    }
    
    display_data
    
  })
  
  output$cov_plot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    chart.Correlation(app_data$result[[chain]]$chain_data, histogram=TRUE)

    
  })
  
  output$traceplot <- renderPlot({
    chain = as.numeric(input$select_chain)
    
    ## Has to be generalized 
    
    if(input$choose_PK_mod==1) {
        gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,    
                                app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2, 
                                app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, 
                                app_data$result[[chain]]$p_iter_ETA4, app_data$result[[chain]]$p_dens_ETA4, nrow=4, ncol=2)  
    } else if(input$choose_PK_mod==2) {
        gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,    
                                app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2, 
                                app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, 
                                app_data$result[[chain]]$p_iter_ETA4, app_data$result[[chain]]$p_dens_ETA4,
                                app_data$result[[chain]]$p_iter_ETA5, app_data$result[[chain]]$p_dens_ETA5,
                                app_data$result[[chain]]$p_iter_ETA6, app_data$result[[chain]]$p_dens_ETA6, nrow=6, ncol=2) 
    }
    
    
  })
  
  
  output$modelfile <- renderText({

      if(input$choose_PK_mod==1 & !input$choose_SS) {  ## Prepared for the simulation with 2-cmt
          includeText("1cmt_multiple_dose.bug")
        } else if (input$choose_PK_mod==2 & !input$choose_SS) {
          includeText("2cmt_multiple_dose.bug")
        } else if (input$choose_PK_mod==1 & input$choose_SS) {
          includeText("1cmt_ss.bug")
        } else if (input$choose_PK_mod==2 & input$choose_SS) {
          includeText("2cmt_ss.bug")
        }
  })
  
  observeEvent(input$submit,
               app_data$pk_plots <- updatePKPlot()
               )
  
  
  
})


