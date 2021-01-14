#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Darren Green
# 14/01/2021

library(shiny)

ui <- fluidPage(
        theme="styles.css",
        titlePanel("Sample size calculation for freedom from disease with imperfect testing"),
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
            column(3,numericInput("TSens",NULL,1,min=0,max=1,step=0.01)),
            column(3,"Test specificity"),
            column(3,numericInput("TSpec",NULL,1,min=0,max=1,step=0.01))
        ),
        h4("Testing targets"),
        fluidRow(
            column(3,"Herd sensitivity"),
            column(3,numericInput("HSens",NULL,0.95,min=0,max=0.98,step=0.01)),
            column(3,"Herd specificity"),
            column(3,numericInput("HSpec",NULL,0.95,min=0,max=0.98,step=0.01))
        ),
        actionButton("submit","GO!"),
        #    actionButton("stop","STOP!"),
        h4("Model output"),
        fluidRow(
            column(3,"Herd sensitivity"),
            column(3,textOutput("oHSens")),
            column(3,"Herd specificity"),
            column(3,textOutput("oHSpec"))
        ),
        fluidRow(
            column(3,"Sample size"),
            column(3,textOutput("on")),
            column(3,"Cutpoint number of reactors"),
            column(3,textOutput("oc"))
        ),
        textOutput("Ioo"),
        plotOutput("plot",width="75%"),
        img(src='ioa_logo.png',style="width: 256px; align: left; margin-right: 2em"),
        "Darren Green (2021), with inspiration from www.ausvet.com.au",
        img(src='parasite_2.png',style="width: 64px; align: right; margin-left: 2em")
        
)

# Define server logic required to draw a histogram
server <- function(input, output) {
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
                          ".")
                }
            }
        })
        
        output$on <- renderText({
            if (D$active) D$results[1]
        })
        
        output$oc <- renderText({
            if (D$active) D$results[2]
        })
        
        output$oHSens <- renderText({
            if (D$active) format(D$results[3],digits=4)
        })
        
        output$oHSpec <- renderText({
            if (D$active) format(D$results[4],digits=4)
        })
        
        output$plot <- renderPlot({
            if (D$active) DrawROC()
        })
}

# Run the application 
shinyApp(ui = ui, server = server)
