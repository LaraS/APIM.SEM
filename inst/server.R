
shinyServer(function(input, output, session) {
  
  ####### New analysis / reload button ########
  output$reload <- renderUI({
    if (input$newanalysis > 0) {
      tags$script("window.location.reload();")
    }
  })
  
  ######## Reactive Data Input ########
  dataInput <- reactive({
    inFile <- input$file1
    exdata <- input$exdata
    
    if(is.null(inFile)){      
      if(exdata==""){
        return(NULL)        
      }else if(exdata=="campbell"){
        return(campbell)  
      }      
    }
    
    if(!is.null(inFile)){
      
     # Determine document format;
     ptn <- "\\.[[:alnum:]]{1,5}$"
     suf <- tolower(regmatches(inFile$name, regexpr(ptn, inFile$name)))
      
      if(suf == ".csv"){
        return(read.csv(inFile$datapath))  
      }else if(suf == ".txt"){
        return(read.table(inFile$datapath))    
      }else if(suf == ".sav"){
        return(read.spss(inFile$datapath, to.data.frame=TRUE,
                         use.value.labels=input$vallabels))
      }else if(suf == ".xpt"){
        return(read.xport(inFile$datapath))
      }  
    }  
  })
  
  ###### Reactive Run Model #########
  model <- reactive({
        
    ## arguments for apim_sem()
    xname <- input$xname
    yname <- input$yname 
    namex1 <- input$namex1
    namex2 <- input$namex2
    namey1 <- input$namey1
    namey2 <- input$namey2
    lab1 <- input$lab1
    lab2 <- input$lab2
    labs1 <- input$labs1
    labs2 <- input$labs2
    dname <- input$dname
    data <- dataInput()
              
    tryCatch(
      apim_sem(xname=xname,
               yname=yname, 
               namex1=namex1, 
               namex2=namex2, 
               namey1=namey1, 
               namey2=namey2,
               lab1=lab1, 
               lab2=lab2, 
               labs1 = labs1, 
               labs2=labs2,
               dname=dname, 
               data=data, 
               relx = .8, 
               rely = .8, 
               alpha = .05, 
               dvarr=TRUE, 
               fix1k = TRUE, 
               boottrials=50L, 
               fix2k = TRUE)
    )
  })

  
  
  
  ###### Output Data Table #########  
  output$mytable1 = renderDataTable({ 
    d <- dataInput()
    dprint <- format(d, digits=3)
    dprint
  })

  #### Main Output ###########    
  output$summary = renderPrint({
    
    m1 <- model()
    print(m1)  
  })

  #### Table Output ###########    
  output$tableoutput = renderPrint({
    
    m1 <- model()
    printTables(m1)  
  })
  
  

})