#server.r
#Darren Green
#21/05/2019


library(shiny)
#library(sinib)

shinyServer(
  
  function(input,output) {
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

  }
  
)
