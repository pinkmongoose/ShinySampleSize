#ui.r
#Darren Green
#21/05/2019


library(shiny)
shinyUI(
  
  fluidPage(
    titlePanel("Sample size calculation for freedom from disease with imperfect testing"),
    h4("Population parameters"),
    fluidRow(
      column(3,"Population size"),
      column(3,numericInput("N",NULL,50,min=1,max=25000,step=1)),
      column(3,"Percentage prevalence / number of reactors"),
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
    "Darren Green (2019), with inspiration from www.ausvet.com.au",
    img(src='parasite_2.png',style="width: 64px; align: right; margin-left: 2em")
  )
  
)
