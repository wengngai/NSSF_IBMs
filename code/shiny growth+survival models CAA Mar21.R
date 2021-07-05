library(shiny)


server <- function(input, output) {
  
  ## MORTALITY MODEL
  d1.2 <- reactive({seq(-5, input$Dthresh.2, len=200)})
  d2.2 <- reactive({seq(input$Dthresh.2, input$Dmax.2, len=200)})
  S1.2 <- reactive({input$K.2 / (1 + exp(-input$r1.2*(d1.2()-input$p1.2)))})
  S2.2 <- reactive({input$K.2 / (1 + exp(input$r2.2*(d2.2()-input$p2.2)))})
  output$Plot2 <- renderPlot({
    plot(S1.2() ~ d1.2(),
         ylab="Survival probability",
         xlab="DBH (log-transformed)",
         xlim=c(-5,input$Dmax.2), ylim=c(0,1),
         cex.lab=1.5, las=1, type="n")
    lines(S1.2() ~ d1.2(),
          lwd=3, col="forestgreen")
    lines(S2.2() ~ d2.2(),
          lwd=3, col="darkgreen")
  })
  
  ## GROWTH MODEL
  dbh1 <- reactive({seq(0, input$Dmax.1, len=500)})
  Y1 <- reactive({(input$alpha * (dbh1() ^ input$beta)) * exp(-input$gamma*dbh1())})
  output$Plot1 <- renderPlot({
    plot(Y1() ~ dbh1(), type="l",
         lwd=3, col="steelblue",
         ylab="Annual growth rate (log-transformed)",
         xlab="DBH (cm)",
         cex.lab=1.5, las=1,
         ylim=c(0,1))
  })
  
}

# Define UI for slider demo app ----
ui <- fluidPage(
    # MORTALITY: Two-part model of Needham et al (2018)
    fluidRow(
      # App title ----
      titlePanel('Survival function (Needham et al 2018)'),
      # Sidebar layout with input definitions ----
      sidebarLayout(
        sidebarPanel(
          # Input: Specification of parameters ----
          sliderInput("r1.2", HTML(paste0("Rate 1 (r", tags$sub("1"), ")")),
                      min = 0, max = 8,
                      step = 0.01, value = 2,
                      animate = animationOptions(interval = 20)),
          sliderInput("p1.2", HTML(paste0("Inflection point 1 (p", tags$sub("1"), ")")),
                      min = -6, max = 0,
                      step = 0.1, value = -4.4,
                      animate = animationOptions(interval = 20)),
          sliderInput("K.2", "K",
                      min = 0.9, max = 1,
                      step = 0.002, value = 0.96,
                      animate = animationOptions(interval = 5)),
          sliderInput("r2.2", HTML(paste0("Rate 2 (-r", tags$sub("2"), "; negative value)")),
                      min = 0, max = 10,
                      step = 0.01, value = 3.8,
                      animate = animationOptions(interval = 20)),
          sliderInput("p2.2", HTML(paste0("Inflection point 2 (p", tags$sub("2"), ")")),
                      min = 2, max = 8,
                      step = 0.1, value = 4,
                      animate = animationOptions(interval = 20)),
          sliderInput("Dmax.2", "Maximum DBH",
                      min = log(40), max = 6,
                      step = 0.01, value = 5),
          sliderInput("Dthresh.2", "DBH threshhold",
                      min = 0, max = 3,
                      step = 0.1, value = 1)
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          # Display equation ----
          withMathJax(),
          titlePanel('$$S = \\frac{K}{1+ e^{-r(D-p)}}$$'),
          
          # Output: Table summarizing the values entered ----
          plotOutput("Plot2")
        )
      )),
    
    # GROWTH: gamma function (Kohyama et al 2020)
    fluidRow(
      # App title ----
      titlePanel('Gamma function (Koyhama et al 2020)'),
      # Sidebar layout with input definitions ----
      sidebarLayout(
        sidebarPanel(
          # Input: Specification of parameters ----
          sliderInput("alpha", HTML("&alpha;"),
                      min = 0, max = 0.5,
                      step = 0.005, value = 0.05,
                      animate = animationOptions(interval = 50)),
          sliderInput("beta", HTML("&beta;"),
                      min = 0, max = 2,
                      step = 0.01, value = 1,
                      animate = animationOptions(interval = 50)),
          sliderInput("gamma", HTML("&gamma;"),
                      min = 0, max = 0.2,
                      step = 0.001, value = 0.05,
                      animate = animationOptions(interval = 50)),
          sliderInput("Dmax.1", "Maximum DBH",
                      min = 10, max = 150,
                      step = 1, value = 60)
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          # Display equation ----
          withMathJax(),
          titlePanel('$$\\alpha D^\\beta e^{-\\gamma D}$$'),
          
          # Output: Table summarizing the values entered ----
          plotOutput("Plot1")
        )
      ))
  
)

shinyApp(ui = ui, server = server)


