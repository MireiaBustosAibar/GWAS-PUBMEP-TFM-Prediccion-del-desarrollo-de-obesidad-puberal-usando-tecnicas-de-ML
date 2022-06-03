## Features names 
data = read.table("www/datos_score.csv",row.names=1,header=TRUE)
names1 = rownames(data)
names = names1[2:length(names1)]

credentials <- data.frame(
  user = ("mireiaTFM2022"), # mandatory
  password = ("PRS_obesity_2022"), # mandatory
  start = c("2022-05-16"), # optinal (all others)
  expire = c(NA, "2022-09-30"),
  admin = c(FALSE, TRUE),
  comment = "Simple and secure authentification mechanism 
  for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
)


# UI ----
ui <- dashboardPage( 
  ## HEADER ----
  dashboardHeader(
    title="PRSity",
    titleWidth = 130),
  ## MENU ----
  dashboardSidebar(
    ## SECTIONS ----
    menuItem("Introduction", tabName = "intro", icon=icon("info-circle"), 
      menuSubItem("What is PRSity?", tabName = "what"),
      menuSubItem("How to use it?", tabName = "use")),
    menuItem("Importing external data", tabName = "data", icon = icon("database")), # database icon
    menuItem("PRS Evaluation", tabName = "importance", icon = icon("table")),  # table icom
    menuItem("ML Classification results", tabName = "classification", icon = icon("poll")), # poll icon
    menuItem("Plots", tabName = "rf", icon = icon("server"), badgeLabel = "new", badgeColor = "green")
    ), 
  ## BODY ----
  dashboardBody(
    ## EVERY ITEMS ----
    tabItems(
      # ITEM: INTRODUCTION
      # SUBITEM1: -----
      tabItem(tabName = "what",
        h1("What is PRSity?"),
        br(),
        p("PRSity (short for construction of machine learning models based on PRSs and GRSs and environmental variables for prediction of severe obesity in children) 
        is an interactive web application to predict the puberal development of obesity and its comorbidities,
          using polygenic risk scores (PRSs) computed through the analysis of genomic data from children. This application is based on a Machine Learning model known as C4.5. 
          Notably, we have used", tags$code("Caret"), "R package to compute this model.", tags$code("Caret"), "stands for Classification And Regression Training,
          more information of its use can be found in", tags$a(href="https://topepo.github.io/caret/", 
                                                                                            "Caret"), 
          "is a popular R package which could use several algorithms though their",  
          tags$a(href= "https://www.rdocumentation.org/packages/caret/versions/4.47/topics/train", "train()"), "function."),
        p("The design and development of PRSity is a part of the Master thesis", tags$em("Predicción de la susceptibilidad hereditaria del desarrollo de obesidad puberal usando técnicas de Machine Learning"), "by Mireia Bustos-Aibar. During this work, it has been selected 
          different powerful and understandable Machine Learning models such as:"), 
        tags$ol(
          tags$li("Decision trees: C4.5 (also known as J48 in R)."),
          tags$li("Rule based models: C5.0Rules."),
          tags$li("Plots models: rf and parRF.")
        ),
        p("These algorithms, which are able to learn about the data and create predictive models, are implemented in 
          ", tags$a(href="https://topepo.github.io/caret/", "caret"), "R package. Moreover, these selected models are 
          understandable because of their explainability and/or interpretability."),
        p("The Machine learning models have been created using a genomic dataset from children in pubertal stages. These data belong to the PUBMEP project.
          The main objetive of this study is the use of Machine Learning models that could learn about the genomic data to predict 
          the risk of develop obesity and/or metabolic syndrome during puberty."),
        p("These data have been collected during the development of a research project that is named",
          tags$a(href="https://uceens.ugr.es/proyecto/pubertad-y-riesgo-metabolico-en-ninos-obesos-alteraciones-epigeneticas-e-implicaciones-fisiopatologicas-y-diagnosticas-estudio-pubmep/","PUBMEP"),".
          These data have common problems such as high heteronenity, high dimensionality, imbalance class and 
          missing values presence. To solve the previous problems it has been implemented a repeated k fold cross validation (5x5) 
          to minimize the error estimation, it has been used a ", tags$a(href="https://www.rdocumentation.org/packages/missForest/versions/1.4/topics/missForest", "missForest"),  
          "R package to impute missing values and 
           it has been used an oversampling method to improve the model performance."),
        #img(src = 'omics_approaches.png', height = 300, width = 900),
        #HTML('<div>
        #<img src="genobox.png" style="float: left">
        #<img src="pubmep.png" style="float: right">
        #<img src="omics_approaches.png" width="1400" style"float: center">
        #</div></pre>'),
#        HTML('<center><img src="omics_approaches.png"width="1200" height="500"></center>')
      ),
      # SUBITEM2-----
      tabItem(tabName = "use",
         h1("How to use PRSity?"),  
         h2("Concept of this web application"),
         p("Further information on the concept and objectives of this web application can be found in the introductory section", tags$b("What is PRSity?")),
        
         h2("ML Classification results"), 
         p("This section have been created to report about the obtained classification metrics
           using the selected Machine Learning algorithms and oversampling methods.
           The user can select in the", tags$b("Option panel"), "seeing models results (Classification metrics by models)." 
           ),
         h2("Importing external data"),
         p("In order to ensure the repeatability of the calculation of the models developed in this work, the users can import their own data to make predictions using the C4.5 model. This is the function of
           the", tags$em("Importing external data"), "section. Therefore, user data can be uploaded in the panel", tags$em("Uploading data"),". In this same place, the user 
           can also choose different different options related to the formatting of input data; such as header, separator or quote.",tags$b("The input file must be a csv file"),".  
           When the data have been loaded, the first features and samples can be visualized. The exact number of features to be visualized can be selected by the user using the option bar."
           ), 
         h2("PRS Evaluation"), 
         p("The computation and evaluation of PRSs and GRSs are performed with PRSice-2 as is described in the documentation of its main web-site:"
           ),
         HTML('<center><img src="PRS_calculation.png"</center>')
         
         ),
      
      # ITEM0: CLASSIFICATION RESULTS ----
      tabItem(tabName = "classification",
              h1("ML Classification results using PRSs"),
              # ROW FROM ITEM0
              fluidRow(
              # ROW FROM ITEM0=BOX1 
                box(
                  title = "Classification metrics by models",
                  tableOutput("results")                        # OUTPUT: table
                ),
                # ROW FROM ITEM0=BOX2
                box(
                  title = "Options",
                  # Input1: Select the options
                  radioButtons("under", "Select the undersampling method:",                #radioButtons
                               choices = c("1) Without overrsampling" = 1,
                                           "2) With oversampling" = 2),
                               selected = 2),
              
                collapsible = TRUE,))
      ),
      # ITEM1: IMPORT YOUR DATA ----
      tabItem(tabName = "data",
              h1("Importing external data"),
              # ROW FROM ITEM1
              fluidRow(
                # ROW FROM ITEM1=BOX1 (head(data))
                box(title="Your first examples",
                    tableOutput("head_data")),                  # OUTPUT1: table
                # ROW FROM ITEM1=BOX2 (importdata))
                box(
                  title = "Uploading data",
                  
                  # Input1: Box to upload the file
                  fileInput("dataset", label="Upload user data", #fileInput
                      buttonLabel = "Upload local file", 
                      placeholder = "File to import", 
                      accept = c(".csv")),
                  # Validation of data importation
                  verbatimTextOutput("validation"), 
                  # Input2: Checkbox if file has header 
                  checkboxInput("header", "Header", TRUE),      #checkboxinput
                  
                  # Input3: Select separator
                  radioButtons("sep", "Separator",              #radioButtons
                               choices = c(Comma = ",",
                                           Semicolon = ";",
                                           Tab = "\t"),
                               selected = ","),
                  # Input4: Select quotes 
                  radioButtons("quote", "Quote",                #radioButtons
                               choices = c(None = "",
                                           "Double Quote" = '"',
                                           "Single Quote" = "'"),
                               selected = '"'),
                  # Horizontal line
                  tags$hr(),
                  
                  # Input5: Number of examples which user want to see
                  sliderInput("n_examples",                     # sliderInput
                              "Number of examples", 1, 30, 10), 
                collapsible = TRUE)
              )
      ),
      # ITEM3: Plots predictions ----
      tabItem(tabName = "rf",
              h1("Final Plots"),
              # ROW FROM ITEM2
              fluidRow(
                # ROW FROM ITEM2 = BOX1 
                box(title="Plots",
                    tableOutput("predict")), # OUTPUT1: table
                # ROW FROM ITEM2 = BOX2
                box(title="Plotting options",
                    # Input1
                    radioButtons("predict_options", "How do you want to get the prediction?",                #radioButtons
                                 choices = c("1) Bar-plot" = 1,
                                             "2) Box-plot" = 2), 
                                 selected = 1),
                    # Input2
              )
              )
      )
    )
))


# Packages ----
library(shiny)
library(shinydashboard)
library(caret)
library(randomForest)
library(ggplot2)
library(openxlsx)
library(shinymanager)
library(rsconnect)

# Wrap your UI with secure_app
ui <- secure_app(ui)

#SERVER ----
server <- function(input, output, session){
  
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )
  
  output$auth_output <- renderPrint({
    reactiveValuesToList(res_auth)
  })
  # ITEM1: OUTPUT1 ----
  output$head_data = renderTable({
    if(is.null(input$dataset)){
      return(NULL)
    }
    else{ # ITEM1: input1,2,3,4,5 (dataset, header, sep, quote, n_examples)
      dat = read.csv(input$dataset$datapath, 
                     header=input$header, 
                     sep=input$sep, 
                     quote=input$quote)}
      head(dat[, 1:8], input$n_examples)
  })
  output$validation = renderPrint(
    if(is.null(input$dataset)){
      print("")
    }else{
       dat = read.csv(input$dataset$datapath, 
                     header=input$header, 
                     sep=input$sep, 
                     quote=input$quote)
       correct = colnames(dat)==names1
       correct2 = unique(correct)
       paste("Correct data importation:", correct2, collapse="")
       }
    )
  
  # ITEM3: OUTPUT1,2 ----
  output$predict = renderTable({
    dat = read.csv(input$dataset$datapath, 
                   header=input$header, 
                   sep=input$sep, 
                   quote=input$quote)
    predict1 = predict(model_RF, dat, type="response")
    predict1 = as.character(predict1)
    predict2 = predict(model_RF, dat, type="prob")
    predict3 = cbind(Response = predict1, predict2)
    if(input$predict_options==1){
      head(cbind(ID = dat$X, Predict = predict1), input$n_examples2)
    }else{
      head(cbind(ID = dat$X, predict2), n=input$n_examples2)
    }
  })
  

}

# INTEGRATION INTERFACE (UI) and INTERACTIVE PART (SERVER) ----
shinyApp(ui, server)
