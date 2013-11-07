#Upgraded Shiny-Server
#sudo su - -c "R -e \"install.packages('devtools', repos=c(RStudio='http://rstudio.org/_packages', CRAN='http://cran.rstudio.com'))\""
#sudo su - -c "R -e \"devtools::install_github('shiny', 'rstudio')\""
#restart shiny-server

library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Melanoma Rapid Learning Utility")
  
  # Sidebar Panel
   ,sidebarPanel(
     
    # Group Plot
    selectInput("groupBy", "Compare Outcomes By:",
                list("Drug Class" = "DRUG_CLASS",
                     "Drug Name" = "MEDICATION",
                     "Mutation Type" = "MUTATION",
                     "BRAF Status" = "BRAF",
                     "NRAS Status" = "NRAS",
                     "Sex" = "SEX"
                ))    
    ,br(),br()
    # Minimum Patients
    ,checkboxInput("includeMinGroupSize", "Minimum Patients Per Outcome Group", FALSE),
      conditionalPanel(
        condition = "input.includeMinGroupSize == true",
        numericInput("minGroupSize", "Require n >=", 12)
      )
    # INSTITUTION
    ,checkboxInput("includeSource", "Filter By Institution of Origin", FALSE),
    conditionalPanel(
      condition = "input.includeSource == true",
      selectInput("selectSource", "Include Patients From:",
                  list("Stanford" = "STANFORD",
                       "Vanderbilt" = "VANDERBILT"))
    )

    # SEX
    ,conditionalPanel(
      condition = "input.groupBy != 'SEX'",
      checkboxInput("includeSex", "Filter By Sex", FALSE),
      conditionalPanel(
        condition = "input.includeSex == true",
        uiOutput("selectSex")
      )
    )
    ,
    # AGE
    conditionalPanel(
    condition = "input.groupBy != 'AGE'",
    checkboxInput("includeAge", "Filter By Age", FALSE),
    conditionalPanel(
      condition = "input.includeAge == true",
      uiOutput('selectAge'),
      uiOutput("age_range_slider")
    )
    ),
    
    # BRAF
    conditionalPanel(
    condition = "input.groupBy != 'BRAF'",
    checkboxInput("includeBRAF", "Filter By BRAF Status", FALSE),
     conditionalPanel(
      condition = "input.includeBRAF == true",
      uiOutput('selectBRAF')
    )
    ),
    
    # NRAF
    conditionalPanel(
    condition = "input.groupBy != 'NRAS'",
    checkboxInput("includeNRAS", "Filter By NRAS Status", FALSE),
    conditionalPanel(
      condition = "input.includeNRAS == true",
      uiOutput('selectNRAS')
    )
    ),
    
    # Filter Drug Class
    checkboxInput("includeClasses", "Filter By Drug Class", FALSE),
    conditionalPanel(
      checkboxInput("twoClass", "Compare First Two Treatments", FALSE),
      condition = "input.includeClasses == true",
      uiOutput("drug_classes")
      ,conditionalPanel(
      condition = "input.twoClass == true"
      ,radioButtons("selectAllNoneClasses", "Quick Check/Uncheck All Classes: ",
                   list("Check All" = 'all', 
                        "Uncheck All" = "none"
                   ))
      )
    ),
 
    # Filter Drug Name
    checkboxInput("includeDrug", "Filter By Drug Name", FALSE),
    conditionalPanel(
      condition = "input.includeDrug == true",
      uiOutput("drug_names"),
      radioButtons("selectAllNoneDrugs", "Quick Check/Uncheck All Drugs: ",
                   list("Check All" = 'all', 
                        "Uncheck All" = "none"
                   )),
      br()
    )
    )
  
   ,mainPanel(
    
     tabsetPanel(
       tabPanel("Cohort Summary"
                 ,plotOutput("cohortAgeSex", height="auto")
                 ,plotOutput("cohortGen", height="auto")
                 ,plotOutput("cohortDrugs", height="auto")),
       tabPanel("Outcomes"
                ,tabsetPanel(
                  tabPanel("Kaplan Meier",
                           plotOutput("survCurv", height="auto")
                           ,h5("Model Summary")
                           ,verbatimTextOutput("survSummary")
                           #,h5("Differences Statistics")
                           #,verbatimTextOutput("survDiffSummary")
                           )
                  ,tabPanel("Time to Next Treatment",
                            conditionalPanel(
                              condition = "input.groupBy == 'MEDICATION'",
                              HTML("<center>"),
                              plotOutput("boxResponseMed", height="auto",width="100%"),
                              HTML("</center>")
                              )
                           ,conditionalPanel(
                             condition = "input.groupBy != 'MEDICATION'",
                             HTML("<center>"),
                             plotOutput("boxResponse", height="auto",width="60%"),
                             HTML("</center>")
                           )
                          ,verbatimTextOutput("boxSummary")
                           )
                  ,tabPanel("Tumor Response", plotOutput("barResponse", height="auto"))
                )
                ),
       tabPanel("Raw Data (Cohort Specific)",
                tabsetPanel(
                  tabPanel("Clinical", helpText('All MRNS are de-identified.\n
                                    "Days to Second Treatment" only
                                    relevant to patients who received more than one drug regimine, and represents
                                    the duration (in days) of the first treatment protocol used for the patient
                                    before a new drug was used.\n
                                    Days to Death measured from date of administration of first drug,
                                                and is calculated to the date of data collection for
                                                those patients alive at the end of the study.  
                                                These patients are marked as "Alive at End of Treatment"\n
                                                IMPORTANT: Tumor Change is randomly generated'),
                           downloadButton("downloadDataClinical","Download as .csv"), br(),br(),
                           dataTableOutput("cohortTableClinical")
                           ),
                  tabPanel("Drugs (Distinct)",helpText('All MRNS are de-identified, and all dates are shifted.
                                            Lists includes only the first date in which each distinct cancer drug was used
                                            on each patient.'),
                           downloadButton("downloadDataDrug","Download as .csv"), br(),br(),
                           dataTableOutput("cohortTableDrugs") 
                ),
                  tabPanel("Drugs (Full Data)",helpText('All MRNS are de-identified, and all dates are shifted.
                                            Lists includes only the first date in which each distinct cancer drug was used
                                            on each patient.'),
                           dataTableOutput("cohortTableDrugsFull") 
                  
                )
     )
    )
     )
  )
)
)
