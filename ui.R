




shinyUI(navbarPage("JAGS4TDM - by Oliver Scherf-Clavel (c) 2019 - JMU Wuerzburg",
  
  tabPanel("Main",
      sidebarLayout(
        sidebarPanel(
           checkboxInput ("choose_SS", "Steady state?", value = T),

           h6("Use an Excel File with similar structure (see Table)"),
           fileInput("loadPT", "Load data", multiple = FALSE, 
                     accept = c(".xlsx", ".xls"),                                     
                     width = NULL,buttonLabel = "Browse...", 
                     placeholder = "No file selected"),

           selectInput("choose_PK_mod", "PK Model", list("1 compartment"=1, "2 compartments"=2)),
           numericInput("ka", "Ka [1/h]", value = 0.5),
           numericInput("V", "V [L]", value = 60),
           numericInput("ke", "Ke [1/h]", value = 0.09),
           numericInput("F_oral", "Systemically available fraction (F)", value = 0.9),
           conditionalPanel(condition="input.choose_PK_mod==2",
                            numericInput("k12", "K12 [1/h]", value = 0.09),
                            numericInput("k21", "K21 [1/h]", value = 0.09)
           ),
           numericInput("tlag", "Lagtime [h]", value = 0.5),
           numericInput("omega1", "Variance of ETA1 (ka)", value = 0.51),
           numericInput("omega2", "Variance of ETA2 (V)", value = 0.54),
           numericInput("omega3", "Variance of ETA3 (ke)", value = 0.53),
           numericInput("omega4", "Variance of ETA4 (F)", value = 0.3),
           conditionalPanel(condition="input.choose_PK_mod==2",
                            numericInput("omega5", "Variance of ETA5 (k12)", value = 0.09),
                            numericInput("omega6", "Variance of ETA6 (k21)", value = 0.09)
           ),
           numericInput("sigma", "Additive Error (mg/L)²", value=0.15),
           sliderInput("TIME", "Time to simulate", value = c(0,72), min=0,max=240),
           sliderInput("n.mc",
                       "Number of MC Simulations:",
                       min = 500,
                       max = 5000,
                       value = 1000),
           sliderInput("mcmc_n.iter",
                       "Iterations per chain:",
                       min = 500,
                       max = 5000,
                       value = 1000),
           sliderInput("mcmc_n.burn",
                       "Burn-in iterations per chain:",
                       min = 0,
                       max = 500,
                       value = 200),
           actionButton("submit", label = "Apply Changes", icon = NULL, width = NULL)
        ),
        
        
        
        
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("pkPlot", height = 800),

            wellPanel("Dataset used in the computation",
              tableOutput("data_set")
            )
         )
      )
  ),
  tabPanel("Diagnostic plot MCMC",
           wellPanel("Results of MCMC chain 1",
              plotOutput("chain1")
              ),
           wellPanel("Results of MCMC chain 2",
              plotOutput("chain2")
              ),
           wellPanel("Results of MCMC chain 3",
              plotOutput("chain3")
              ),
           wellPanel("Results of MCMC chain 4",
              plotOutput("chain4")
              )
           ),
  tabPanel("Model File",
           verbatimTextOutput("modelfile"))
        )
)
