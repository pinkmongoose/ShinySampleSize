#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Darren Green
# 28/11/2024

library(shiny)

ui <- fluidPage(
        theme="styles.css",
        titlePanel("Sample size calculation for freedom from disease with imperfect testing"),
        
        tabsetPanel(
                tabPanel("Parameters", 
                        h4("Population parameters"),
                        fluidRow(
                                column(3,"Population size"),
                                column(3,numericInput("N",NULL,50,min=1,max=25000,step=1)),
                                column(3,"Proportion prevalence / number of reactors"),
                                column(3,numericInput("R",NULL,0.1,min=0,max=25000,step=1))
                        ),
                        h4("Test parameters"),
                        fluidRow(
                                column(3,"Test sensitivity"),
                                column(3,numericInput("beta",NULL,1,min=0,max=1,step=0.01)),
                                column(3,"Test specificity"),
                                column(3,numericInput("alpha",NULL,1,min=0,max=1,step=0.01))
                        ),
                        h4("Testing targets"),
                        fluidRow(
                                column(3,"Herd sensitivity"),
                                column(3,numericInput("Breq",NULL,0.95,min=0,max=0.98,step=0.01)),
                                column(3,"Herd specificity"),
                                column(3,numericInput("Areq",NULL,0.95,min=0,max=0.98,step=0.01))
                        ),
                        h4("Model output"),
                        fluidRow(
                                column(3,"Herd sensitivity"),
                                column(3,textOutput("B")),
                                column(3,"Herd specificity"),
                                column(3,textOutput("A"))
                        ),
                        fluidRow(
                                column(3,"Sample size"),
                                column(3,textOutput("n")),
                                column(3,"Cutpoint number of reactors"),
                                column(3,textOutput("cutpoint"))
                        )
                        
                ),
                tabPanel("Hyper-parameters",
                        h4("Test parameters"),
                        fluidRow(
                                column(3,selectInput("Testtype", "Diagnostic test model type",
                                                      choices = c('Binomial'='Binomial',"Beta Binomial"="BetaBinomial"))),
                                column(2,"Sensitivity sample size (if not 'binomial')"),
                                column(2,numericInput("betaN",NULL,1000000,min=1,step=1)),
                                column(2,"Specificity sample size (if not 'binomial')"),
                                column(2,numericInput("alphaN",NULL,1000000,min=1,step=1))
                        ),
                        fluidRow(
                                column(3,selectInput("Sampletype","Sample type",
                                                     choices = c("Without replacement"="Hypergeometric","With replacement"="Binomial")))
                        ),
                        h4("Model tuning"),
                        fluidRow(
                                column(3,"Coverage for convolution"),
                                column(2,numericInput("pConv",NULL,0.999,min=0.99,max=1,step=0.0001))
                        )
                ),
                tabPanel("Presets",
                        h4("Presets"),
                        fluidRow(
                                 column(3,selectInput("preset","Preset",choices = c(
                                         'Stirling MSc default scenario'='MSc',
                                         'Stirling MSc scenario 2'='MSc2',
                                         'Johnson et al. 2003 FMD survey'='Johnson',
                                         'Cameron et al. 1998 FMD survey'='Cameron'
                                         )))
                        ),
                        fluidRow(
                                column(2,actionButton("presetSubmit","SET!"))
                        )
                        
                )
        ),
        hr(),
        actionButton("submit","GO!"),
        textOutput("Ioo"),
        plotOutput("plot",width="75%"),
        img(src='ioa_logo.png',style="width: 256px; align: left; margin-right: 2em"),
        "Darren Green (2021-4)",
        img(src='parasite_2.png',style="width: 64px; align: right; margin-left: 2em")
        
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
        D <- reactiveValues()
        D$active <- F
        D$err <- F
        D$warn <- ""
        source("sampling.r",local=T)
        
        observeEvent(input$stop, {
                showModal(modalDialog(title="Error","Stop"))
        })
        
        observeEvent(input$submit, {
                D$err <- F
                if (!D$err) withProgress(
                        message="Calculating sample size",
                        detail="Getting there...",
                        value=0, min=0, max=input$N,
                        {RunModel()}
                )
                if (!D$err) D$active <- T
        })
        
        observeEvent(input$presetSubmit, {
                showModal(modalDialog(title="Preset","All parameters have been updated"))
                updateNumericInput(session,"Breq",value=0.95)
                updateNumericInput(session,"Areq",value=0.95)
                updateNumericInput(session,"pConv",value=0.999)
                if (input$preset=="MSc") {
                        updateNumericInput(session,"N",value=50)
                        updateNumericInput(session,"R",value=0.1)
                        updateNumericInput(session,"beta",value=1)
                        updateNumericInput(session,"alpha",value=1)
                        updateNumericInput(session,"betaN",value=1000000)
                        updateNumericInput(session,"alphaN",value=1000000)
                        updateSelectInput(session,"Testtype",selected="Binomial")
                } else if (input$preset=="MSc2") {
                        updateNumericInput(session,"N",value=50)
                        updateNumericInput(session,"R",value=0.1)
                        updateNumericInput(session,"beta",value=0.98)
                        updateNumericInput(session,"alpha",value=0.98)
                        updateNumericInput(session,"betaN",value=1000000)
                        updateNumericInput(session,"alphaN",value=1000000)
                        updateSelectInput(session,"Testtype",selected="Binomial")
                } else if (input$preset=="Johnson") {
                        updateNumericInput(session,"N",value=265)
                        updateNumericInput(session,"R",value=0.3)
                        updateNumericInput(session,"beta",value=68.74/(68.74+4.57))
                        updateNumericInput(session,"alpha",value=107.2/(107.2+3.17))
                        updateNumericInput(session,"betaN",value=68.74+4.57)
                        updateNumericInput(session,"alphaN",value=107.2+3.17)
                        updateSelectInput(session,"Testtype",selected="BetaBinomial")
                } else if (input$preset=="Cameron") {
                        updateNumericInput(session,"N",value=265)
                        updateNumericInput(session,"R",value=0.3)
                        updateNumericInput(session,"beta",value=0.95)
                        updateNumericInput(session,"alpha",value=0.98)
                        updateNumericInput(session,"betaN",value=1000000)
                        updateNumericInput(session,"alphaN",value=1000000)
                        updateSelectInput(session,"Testtype",selected="Binomial")
                }
                
        })
        
        output$Ioo <- renderText({
                if (D$active) {
                        if (nchar(D$warn)) {
                                paste(D$warn)
                        } else {
                                paste("If a sample size of",
                                      D$results[1],
                                      "is taken and",
                                      D$results[2],
                                      "or fewer reactors are found, then the probability that the population is free from disease at a prevalence of",
                                      D$results[6],"/",D$results[5],"(",
                                      format(D$results[6]/D$results[5],digits=4),")",
                                      "is",
                                      format(D$results[3]+0,digits=4),
                                      ". The area under the ROC curve is",
                                      format(D$results[7]+0,digits=3),
                                      ".")
                        }
                }
        })
        
        output$n <- renderText({
                if (D$active) D$results[1]
        })
        
        output$cutpoint <- renderText({
                if (D$active) D$results[2]
        })
        
        output$B <- renderText({
                if (D$active) format(D$results[3],digits=4)
        })
        
        output$A <- renderText({
                if (D$active) format(D$results[4],digits=4)
        })
        
        output$plot <- renderPlot({
                if (D$active) DrawROC()
        })
}

# Run the application 
shinyApp(ui = ui, server = server)
