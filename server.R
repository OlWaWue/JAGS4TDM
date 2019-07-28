


shinyServer(function(input, output) {
   
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
    data_set = data.frame(time=c(0,4,6,12,30,50),
                          amt=c(100,".",".",100,".","."),
                          conc=c(".", 2, 3, ".", 1.5, 0.72),
                          evid=c(1, 0, 0, 1, 0, 0))
  )
  
  output$pkPlot <- renderPlot({
    
    app_data$result = process_data_set(app_data$data_set, n.iter = input$mcmc_n.iter, n.burn = input$mcmc_n.burn,
                                       thetas = c(input$ka, input$V, input$ke, input$tlag),
                                       omegas = c(input$omega1, input$omega2, input$omega3),
                                       TIME =seq(input$TIME[1], input$TIME[2], by=0.2), sigma=input$sigma)
    
    mc_eta2 <- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega2)))
    mc_eta3<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega3)))
    mc_eta1<- (rnorm(n = input$n.mc, mean=0, sd=sqrt(input$omega1)))
    
    
    dat_mc <- NULL
    
    ## Fortschrittsanzeige benutzen

      
    withProgress(message = "Performing Monte Carlo Simulation", max = input$n.mc, {
          for(i in 1:input$n.mc){
    
              CP_mc <-  pk_1cmt_oral(theta = c(input$ka, input$V, input$ke, input$tlag), 
                                                 eta = c(mc_eta1[i], mc_eta2[i],mc_eta3[i]), 
                                                 dosing_events = data.frame(time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$time)),
                                                                             amt=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==1,]$amt))), 
                                                 times=seq(input$TIME[1], input$TIME[2], by=0.2))
              dat_mc <- cbind(dat_mc, CP_mc)
              incProgress(1)
          }
    })
    
    
    s <- apply(t(dat_mc),2,function(x) quantile(x,probs=c(0.1,0.5,0.9)))
    
    plot_dat <- data.frame(TIME=seq(input$TIME[1], input$TIME[2], by=0.2),CP_min=s[1,],CP=s[2,],CP_max=s[3,], DELTA=(s[3,]-s[1,]))
    
    tdm_data <- data.frame(conc=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$conc)),
                           time=as.numeric(as.character(app_data$data_set[app_data$data_set$evid==0,]$time)))
    
    
    ind_plot <- app_data$result[[5]] + theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") +
      ggtitle("MCMC Result including data (posterior)", "80/85/90/95% PI")
    
    pop_plot <- ggplot(data=plot_dat)  + geom_line(aes(x=TIME, y=CP), colour="blue") +
    geom_ribbon(aes(x=TIME, ymax=CP_max, ymin=CP_min), alpha=0.15, fill="blue") +
      theme_bw() + xlab("Time [h]") + ylab("Concentration [mg/L]") + ggtitle("MC Result without data (prior)", "90% PI") +
      geom_point(data=tdm_data, aes(x=time, y=conc))
    
    grid.arrange(pop_plot, ind_plot, nrow=2, ncol=1)
  })

  output$data_set <- renderTable({
    app_data$data_set
  })
  
  output$chain1 <- renderPlot({
    chain = 1
    gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,   
                            app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2,
                            app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, nrow=3, ncol=2)
    
  })
  
  output$chain2 <- renderPlot({
    chain = 2
    gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,   
                            app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2,
                            app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, nrow=3, ncol=2)
    
  })
  
  output$chain3 <- renderPlot({
    chain = 3
    gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,   
                            app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2,
                            app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, nrow=3, ncol=2)
    
  })
  
  output$chain4 <- renderPlot({
    chain = 4
    gridExtra::grid.arrange(app_data$result[[chain]]$p_iter_ETA1, app_data$result[[chain]]$p_dens_ETA1,   
                            app_data$result[[chain]]$p_iter_ETA2, app_data$result[[chain]]$p_dens_ETA2,
                            app_data$result[[chain]]$p_iter_ETA3, app_data$result[[chain]]$p_dens_ETA3, nrow=3, ncol=2)
    
  })
  
  output$modelfile <- renderText({
    includeText("1cmt_multiple_dose.bug")
  })
  
})


