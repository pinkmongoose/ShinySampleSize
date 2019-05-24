#ui.r
#Darren Green
#21/05/2019

library(shiny)
shinyUI(
  
  fluidPage(
    titlePanel(
      "Prevalence in pooled testing"),
    sidebarLayout(
      sidebarPanel(
        helpText("Select values for the model run and then click 'GO!'."),
        h4("Sample sizes"),
        sliderInput("m","size of each pool",min=1,max=250,step=1,value=10),
        sliderInput("n","number of pools",min=1,max=250,step=1,value=10),
        h4("Test values"),
        sliderInput("sens","test sensitivity",min=0.5,max=1,step=0.01,value=1),
        sliderInput("spec","test specificity",min=0.5,max=1,step=0.01,value=1),
        h4("Other parameters"),
        sliderInput("ci","confidence interval",min=0.1,max=0.998,step=0.01,value=0.95),
        actionButton("submit","GO!")
      ),
      mainPanel(
        h4("Sample prevalence versus number of positive pools"),
        plotOutput("plot",height="500px"),
        tableOutput("table")
      )
    ),
    img(src='ioa_logo.png',style="width: 256px; align: left; margin-right: 2em"),
    "Darren Green (2019), with inspiration from www.ausvet.com.au",
    img(src='parasite_2.png',style="width: 64px; align: right; margin-left: 2em")
  )
)
