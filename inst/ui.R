
shinyUI(fluidPage(
  titlePanel(title="Acter Partner Interdependence Model (based on SEM)"),
  
  sidebarLayout(
  sidebarPanel(
    tabsetPanel(
      ######### Data ############
      tabPanel('Data',
        h5("Start a New Analysis"),
        actionButton("newanalysis","New Analysis"),
        uiOutput("reload"),
        br(),
        hr(),        
        h5("Example Data"),
        selectizeInput(inputId="exdata", label="", selected="",
                       choices= c("","campbell"),
                       options = list(placeholder = 'choose example data')),
        br(),        
        hr(),
        h5("Data File"),
        tryCatch(
          fileInput("file1", "", 
                    accept=c(".csv", ".txt", ".sav", ".xpt", 
                              ".CSV", ".TXT", ".SAV", ".XPT"))
        ),
        helpText('Select either a .csv, .txt, .sav or a .xpt file to be uploaded. The corresponding R function (read.csv, read.table, read.spss, or read.xport) will be chosen automatically with the default settings for arguments. To read in data from inside R into the shiny interface, it is easiest to save your dataset using write.csv with the default settings. If reading your SPSS file does not work, please try saving it in SAS Transport format with file ending .xpt. Causes for errors may be special characters in file names and/or path names.')
      ),
      ########## Variable Names ############
      tabPanel('Variable Names',
#                conditionalPanel(
#                   condition = "!input.latenty",
#                   h5("Dependent Variable Y"),
#                   selectizeInput("variabley", "", "",
#                         options = list(placeholder = 'select dependent variable'))
#                ),
              textInput("dname", "Name of Distinguishable Variable", value = "gender"),
              textInput("xname", "Name of Independent Variable", value = "neurotocism"),
              textInput("yname", "Name of Dependent Variable", value = "distress"),
              textInput("namex1", "Name of Independent variabel first role", value = "mneuro"),
              textInput("namex2", "Name of Independent variabel second role", value = "fneuro"),
              textInput("namey1", "Name of dependent variabel first role", value = "mdistr"),
              textInput("namey2", "Name of Name of dependent variabel second role", value = "fdistr"),
              textInput("labs1", "Label first role", value = "husband"),
              textInput("labs2", "Label second role", value = "wive"),
              textInput("lab1", "Label first role (plural)", value = "husbands"),
              textInput("lab2", "Label second role", value = "wives")
      )
  )),
  
  mainPanel(
    tabsetPanel(
      ######### Data Table ##########
      tabPanel('Data', dataTableOutput("mytable1")),
      
      ######### Main Output ##########
      tabPanel("Main Output", verbatimTextOutput("summary")),
      
      ######### Table Output ##########
      tabPanel("Table Output", verbatimTextOutput("tableoutput"))
  ))

)))

