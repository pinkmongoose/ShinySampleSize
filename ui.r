library(shiny)
shinyUI(
  
  fluidPage(
    titlePanel("Sample size calculation for freedom from disease with imperfect testing"),
    h4("Population parameters"),
    fluidRow(
      column(3,"Population size"),
      column(3,numericInput("N",NULL,50,min=1,max=10000,step=1)),
      column(3,"Prevalence (number of reactors)"),
      column(3,numericInput("R",NULL,5,min=1,max=10000,step=1))
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
      column(3,numericInput("HSens",NULL,0.95,min=0,max=1,step=0.01)),
      column(3,"Herd specificity"),
      column(3,numericInput("HSpec",NULL,0.95,min=0,max=1,step=0.01))
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
    textOutput("Ioo")
  )
  
)
