## ui.R ##
library(shinydashboard)
library(shinyBS)  

#################### HEADER ####################
header = dashboardHeader(
  title = "ProteoExplorer",
  # Button pointing to thelab Homepage
  tags$li(class = "dropdown",id = "Saezlab",
          tags$a(href="http://saezlab.org/", target="_blank", 
                 tags$img(icon("file-text"))
          ),
          bsTooltip("Saezlab", "Go to the Saez lab Homepage",
                    "bottom", options = list(container = "body"))),
  
  # Button pointing to the orginal HeLa publiation
  tags$li(class = "dropdown", id="publication",
          tags$a(href="XXXXXX",
                 target="_blank", tags$img(icon("file-text"))
          ),
          bsTooltip("publication", "Go to the publication",
                    "bottom", options = list(container = "body"))),
  
  # Button pointing to the GitHub respository
  tags$li(class = "dropdown", id="github",
          tags$a(href="XXXXXX",
                 target="_blank", tags$img(icon("github"))
          ),
          bsTooltip("github", "View the code of this webside on GitHub",
                    "bottom", options = list(container = "body")))
)

#################### SIDEBAR ####################
sidebar = dashboardSidebar(
  sidebarMenu(id = "menu",
              menuItem("Home", tabName = "Home", icon=icon("home")),
              menuItem("Data Navigation", tabName = "Data_Navigation", icon=icon("upload")),
              menuItem("Results", tabName = "Results", icon=icon("bar-chart")),
              menuItem("Contact", tabName = "Contact", icon=icon("group"))
  )
)

#################### BODY ####################
body = dashboardBody(heigth="auto",width=1, 
                     tabItems(
                       tabItem(tabName = "Home",
                               h1("Welcome to the ProteoExplorer"),
                               h3(tags$i("A web application to visualize Protein abundance prediction based on mRNA and genomic data")),
                               br(),
                               h4("Description:"),
                               p("Even though cancer is driven by genomic alterations, the aberrant cellular phenotypes are frequently the result of dysregulated protein functions. To answer fundamental questions about how different levels of biological signal relate to one another, we launched a collaborative competition: The NCI-CPTAC DREAM Proteogenomics Challenge. The best performance for predicting protein abundances based on mRNA and genetic data was achieved by an ensemble of models including transcript level as proxy, interaction between genes, and conservation across tumor types. From the models in the ensemble, we identified common regulators, which are more essential and predictive of patient outcome. To assess the utility of the protein abundance predictions, we applied them to patient stratification of breast tumors, and identified biomarkers, involved in the immune system, not found using mRNA alone. Those results underline the potential application of these models to large scale proteogenomics characterization of cancer."),
                               br(),
                               h4("The original DREAM Proteogenomics publication can be cited as:"),
                               p("Yang. et al.", 
                                 tags$a("XXXXXX",
                                        href="XXXXXXX",
                                        target="_blank"))
                       ),
                       tabItem(tabName = "Data_Navigation",
                               fluidRow(
                                 box(width = 10 ,
                                     title = tagList(icon("photo"),"User inputs"), status="primary", solidHeader=TRUE,
                                     
                                     tabsetPanel(
                                       tabPanel(
                                         title = tagList(icon("edit"), "Settings"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         helpText(strong("1. Select gene identifier")),
                                         selectInput("select_Identifier", label=NULL, 
                                                     choices = list("HGNC" = 1, "SwissProt ID" = 2), selected = 1)
                                         
                                       ),
                                      
                                       tabPanel(
                                         title = tagList(icon("edit"), "Comparison observed vs predicted"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         helpText(strong("1. Enter 1 gene name")), 
                                         textInput("F1_1_gene", label=NULL),
                                         
                                         helpText(strong("2. Select Testing data")), 
                                         selectInput("F1_layer_X", label=NULL,
                                                     choices = c("True Proteomics (ovarian)"),  ##  "True Proteomics (breast)",
                                                     selected = "True Proteomics (ovarian)", selectize = FALSE),
                                         
                                         helpText(strong("3. Select Prediction output")),
                                         selectInput("F1_layer_Y", label=NULL,
                                                     choices = c("Team Guan, Predicted Proteomics (ovarian)","Ensemble, Predicted Proteomics (ovarian)"),
                                                     selected = "Team Guan, Predicted Proteomics (ovarian)", selectize = FALSE), #  "Team Guan, Predicted Proteomics (breast)",  "Ensemble, Predicted Proteomics (breast)",
                   
                                         actionButton("go_scatter_plot", label="Scatter Plot",icon=icon("bar-chart"))
                                         
                                       ), 
                                       
                                       tabPanel(
                                         title = tagList(icon("edit"), "Prediction performance for a given gene"),
                                         status = "primary", solidHeader = TRUE,
                                         
                                         # Select 1 gene
                                         helpText(strong("1. Enter 1 gene name")),
                                         textInput("F2_1_gene", label=NULL),
                                         
                                         # select omics
                                         helpText(strong("2. Select a team")),
                                         selectInput("F2_1_layer", label=NULL,
                                                     choices = c("Team Guan (breast)","Team Guan (ovarian)"),
                                                     selected = "Team Guan (breast)", selectize = FALSE),
                                         
                                         actionButton("go_density", label="Density plot",icon=icon("bar-chart"))
                                  
                                       )
                                       
                                     ) 
                                     
                                 ),  ## end of box
                                 
                                 tabBox(
                                   width = 10,
                                   title = tagList(icon("table"), "Data available for each input module"),
                                   id = "matrices",
                                   tabPanel("Comparison observed vs predicted", DT::dataTableOutput("F1_layer_X_common")),
                                   tabPanel("Prediction performance for a given gene", DT::dataTableOutput("F2_1_layer"))
                                   
                               )
                               
                       )
                       ),
                       
                       tabItem(tabName = "Results",
                               
                               box(width = 6 ,
                                   title = tagList(icon("photo"),"Scatter plot (observed vs predicted)"), status="primary",
                                   solidHeader=TRUE, plotOutput("scatter_plot", height="500px")
                               ),
                               
                               box(width = 6 ,
                                   title = tagList(icon("photo"),"Prediction performance for a given gene"), status="primary",
                                   solidHeader=TRUE, plotOutput("density_plot", height="500px")
                               )
                               
                       ),
                       
                       tabItem(tabName = "Contact",
                               h3("Contact Us"),
                               p("Please do not hesitate to", 
                                 a("contact us", href="mailto:mi.yang0586@gmail.com"),
                                 " for feedback, questions or suggestions."),
                               br(),
                               h3("On Twitter"),
                               p("You can follow institute news and updates on",
                                 a("@sysbiomed", 
                                   href="https://twitter.com/sysbiomed",
                                   target="_blank"),"."),
                               br(),
                               br(),
                               br(),
                               br(),
                               a(img(src = "logo_saezlab.png", height = 72, width = 72),
                                 href="http://saezlab.org/", target="_blank"
                               ), "Institute for Computational Biomedicine 
                                   BIOQUANT-Zentrum BQ 0053 AG Saez-Rodriguez 
                                   Im Neuenheimer Feld 267, 69120 Heidelberg

                                   Tel: +49 6221 5451334,"
                       )
                     )
                       )

#################### JOINING ####################
tagList(dashboardPage(header, sidebar, body),
        tags$footer("ProteoExplorer, version 0.1 (2019) developed by MI YANG", align = "center", style = "
              position:absolute;
              bottom:0;
              width:100%;
              height:50px;   /* Height of the footer */
              color: black;
              padding: 10px;
              #background-color: black;
              z-index: 1000;"))