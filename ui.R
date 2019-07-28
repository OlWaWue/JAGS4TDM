




shinyUI(navbarPage("JAGS4TDM - by Oliver Scherf-Clavel (c) 2019 - JMU Wuerzburg",
  
  tabPanel("Main",
      sidebarLayout(
        sidebarPanel(
           h6("Use an Excel File with similar structure (see Table)"),
           fileInput("loadPT", "Load data", multiple = FALSE, accept = c(".xlsx", ".xls"),                                     
                     width = NULL,buttonLabel = "Browse...", placeholder = "No file selected"),
           
           
           numericInput("ka", "Ka [1/h]", value = 0.5),
           numericInput("V", "V [L]", value = 60),
           numericInput("ke", "Ke [1/h]", value = 0.09),
           numericInput("tlag", "Lagtime [h]", value = 0.5),
           numericInput("omega1", "Variance of ETA1 (ka)", value = 0.51),
           numericInput("omega2", "Variance of ETA2 (V)", value = 0.54),
           numericInput("omega3", "Variance of ETA3 (ke)", value = 0.53),
           numericInput("sigma", "Residual Variance", value=1.15),
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
           submitButton(text = "Apply Changes", icon = NULL, width = NULL)
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
  tabPanel("Results of MCMC chain 1",
              plotOutput("chain1")),
  tabPanel("Results of MCMC chain 2",
              plotOutput("chain2")),
  tabPanel("Results of MCMC chain 3",
              plotOutput("chain3")),
  tabPanel("Results of MCMC chain 4",
              plotOutput("chain4")),
  tabPanel("Model File",
           verbatimTextOutput("modelfile"))
        )
)
