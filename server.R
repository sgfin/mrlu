require('survival')

#filepath <- "~/Documents/0 Research/Rubin/0 mrlu/"
filepath <- "/home/samfin/ShinyApps/mrlu/"

data <- read.csv(paste(filepath, "clinical_data.csv",sep=""))
pat.drugs <- read.csv(paste(filepath, "pat_drugs.csv", sep=""))
pat.drugs.full <- read.csv(paste(filepath, "drug_data_full.csv", sep=""))

# Function to identify the first "number.treats" drugs used on patient 'MRN'
# Combines PACLI and Carbo when used together
extract_drugs <- function(MRN, number.treats) {
  drug.list <- pat.drugs[which(pat.drugs$MRN == MRN), c("NAME", "ORDERING_DATE")]
  if("PACLITAXEL" %in% drug.list && "CARBOPLATIN" %in% drug.list){
    drug.list <- gsub("PACLITAXEL", "PACL/CARBO", x=drug.list[,"NAME"])
    drug.list <- gsub("CARBOPLATIN", "PACL/CARBO",x=drug.list[,"NAME"])
  }
  
  drug.list <- drug.list[which(!duplicated(drug.list[,"ORDERING_DATE"])) ,"NAME"]
  drug.list <- paste(drug.list[1:min(number.treats,length(drug.list))], collapse=", ")
  return(drug.list)
}

# Identifies the first drug class used on patient 'MRN'
extract_drug_class <- function(MRN, number.treats) {
  d_class <- pat.drugs[which(pat.drugs$MRN == MRN), c("DRUG_CLASS", "ORDERING_DATE")]
  d_class <- as.character(d_class[which(!duplicated(d_class[,"ORDERING_DATE"])) ,"DRUG_CLASS"])
  return(paste( d_class[1:min(number.treats,length(d_class))] , collapse=", "))
}

MEDICATION <- sapply(data$MRN, extract_drugs, 1)
DRUG_CLASS <- sapply(data$MRN, extract_drug_class, 1) # REPLACE WITH SOME MAP
DRUG_CLASS_TWO <- sapply(data$MRN, extract_drug_class, 2) # REPLACE WITH SOME MAP
TUMOR_CHANGE <- runif(length(data$MRN),-1,1)

# Finalizes Data
data <- cbind(data, MEDICATION, DRUG_CLASS, TUMOR_CHANGE)
data2 <- data
data2$DRUG_CLASS <- DRUG_CLASS_TWO

shinyServer(function(input, output) {
    
# Renders UI for Sex Filter
  output$selectSex <- renderUI({
    selectInput("sex", "Sex:",
                list("MALE" = "MALE", 
                     "FEMALE" = "FEMALE"))
  })

# Renders UI for Race Filter
  output$selectRace <- renderUI({
    selectInput("race", "Race:",
                list("White" = 'WHITE', 
                     "Black or African American" = "BLACK OR AFRICAN AMERICAN",
                     "Asian" = "ASIAN"
                ))
  })

# Renders UI for Age Filter
  output$selectAge <- renderUI({
    numericInput("POI.age", "Patient Age:", median(data$AGE))
  })
  
  # Creates Slider for Selecting Age Ranges
  output$age_range_slider <- renderUI({
    age.min <- min(data$AGE)
    age.max <- max(data$AGE)
    
    sliderInput(inputId = "age_range",
                label = paste("Age range"),
                min = age.min, max = age.max,
                value= c(max(age.min, input$POI.age - 15), min(age.max, input$POI.age + 15))
    )
  })
  
# Renders UI for BRAF Filter
  output$selectBRAF <- renderUI({
    selectInput("braf", "BRAF Test Result:",
                list("POSITIVE" = "1", 
                     "NEGATIVE" = "0"))
  })

# Renders UI for NRAS Filter
  output$selectNRAS <- renderUI({
      selectInput("nras", "NRAS Test Result:",
                  list("POSITIVE" = "1", 
                       "NEGATIVE" = "0"))
  })
  
# Renders UI For Drug Name Filter  
  output$drug_names <- renderUI({   
    if(input$selectAllNoneDrugs == 'none') {
      return(checkboxGroupInput("drug_name", "Include:", sort(unique(data$MEDICATION))))
    }
    else {
      if(input$selectAllNoneDrugs == 'all') {
        return(checkboxGroupInput("drug_name", "Include:", sort(unique(data$MEDICATION)), sort(unique(data$MEDICATION))))
      }
    }
    checkboxGroupInput("drug_name", "Include:", sort(unique(data$MEDICATION)), sort(unique(data$MEDICATION)))
  })
  
# Renders UI For Drug Class Filter  
  output$drug_classes <- renderUI({
    if(input$selectAllNoneClasses == 'none') {
      if(input$twoClass){
        return(checkboxGroupInput("drug_class", "Include:", sort(unique(data2$DRUG_CLASS))))
      }
      else{
        return(checkboxGroupInput("drug_class", "Include:", sort(unique(data$DRUG_CLASS))))
      }
          }
    else {
      if(input$twoClass){
        return(checkboxGroupInput("drug_class", "Include:", sort(unique(data2$DRUG_CLASS)), sort(unique(data2$DRUG_CLASS))))
      }
      else{
        return(checkboxGroupInput("drug_class", "Include:", sort(unique(data$DRUG_CLASS)), sort(unique(data$DRUG_CLASS))))
      }
    }
  })
  
#  Function to Narrow Down Dataset According to UI inputs
  selectPats <- function(data, subgroups = TRUE, two_class = FALSE){
    if(two_class){
      plot.data <- data2
    }
    else{
      plot.data <- data
    }
    #EXCLUDE SURVIVORS
    if(input$includeSource){
      plot.data <- subset(plot.data, SOURCE==input$selectSource, drop=T)
    }
    
    if(input$includeSex & input$groupBy != 'sex') {
      plot.data <- subset(plot.data, SEX==input$sex,drop=T)      
    }
    if(input$includeAge) {
      plot.data <- subset(plot.data, AGE >= min(input$age_range) & AGE <= max(input$age_range),drop=T)      
    }
    if(input$includeBRAF) {
      plot.data <- subset(plot.data, BRAF==as.numeric(input$braf),drop=T)      
    }
    if(input$includeNRAS) {
      plot.data <- subset(plot.data, NRAS==as.numeric(input$nras),drop=T)         
    }
    if(input$includeClasses){
      plot.data <- subset(plot.data, DRUG_CLASS %in% input$drug_class, drop=T)
    }
    if(input$includeDrug){
      plot.data <- subset(plot.data, MEDICATION %in% input$drug_name, drop=T)
    }
    if(input$includeMinGroupSize){
      group.sizes <- table(factor(plot.data[,input$groupBy]))
      big.groups <- names(group.sizes)[which(group.sizes > input$minGroupSize)]
      plot.data <- subset(plot.data, plot.data[,input$groupBy] %in% big.groups)
    }
    plot.data
  }
  
  # Build Cox Model, Including Right-Censoring for those patients alive at end of study
  coxModel <- function(plot.data, groupBy){
    switch(groupBy,
           MEDICATION = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(MEDICATION), data=plot.data),
           DRUG_CLASS = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(DRUG_CLASS), data=plot.data),
           SEX = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(SEX), data=plot.data),
           BRAF = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(BRAF), data=plot.data),
           MUTATION = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(MUTATION), data=plot.data),
           NRAS = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(NRAS), data=plot.data)
    )
  }

  # Calculates desired height of plot
  plotHeight <- function(){
    plot.data <- selectPats(data, two_class = input$twoClass)
    distinct = unique(plot.data[,input$groupBy])
    return(length(distinct)*200+100)
  }
  
  nDistinct <-3
  # Generates Survival Plot  
  output$barResponse <- renderPlot({
    plot.data <- selectPats(data, two_class = input$twoClass)
    plot.data <- subset(plot.data, !is.na(plot.data[,input$groupBy]))
    
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
      distinct = unique(plot.data[,input$groupBy])
      if(length(distinct)){
        nDistinct <- length(distinct)
        par(mar=rep(3,4),mgp = c(2, 3, 0), mfrow=c(length(distinct),1))
        for(i in 1:length(distinct)){
          inds <-which(plot.data[,input$groupBy] == sort(distinct)[i])
          barplot(sort(data$TUMOR_CHANGE[inds]),
                  main=paste("Change in Tumor Burden for", sort(distinct)[i]," Patients"),
                  ylab="% Change",
                  xlab="Each Bar Represents a Distinct Patient's Response"
          )
        }
      }
    },height=plotHeight)
  
  

  # Generates Box Plot  
  output$boxResponse <- renderPlot({
    plot.data <- selectPats(data, two_class = input$twoClass)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    name.list <- sort(names(table(factor(plot.data[,input$groupBy]))))
    name.list <- paste(name.list,' (', table(factor(plot.data[,input$groupBy])), ')',sep='')
      if(input$groupBy == 'SEX') {
        plot.data$SEX <- factor(plot.data$SEX)
        boxplot(DAYS_TO_NEXT_RX ~ SEX , plot.data, cex.axis=.65,
                main = 'Time to Next Treatment by Sex',
                ylab = "Days",
                names = name.list,
                xlab="Sex")
      }
      else if(input$groupBy == 'BRAF') {
        plot.data$BRAF <- factor(plot.data$BRAF)
        boxplot(DAYS_TO_NEXT_RX ~ BRAF , plot.data, cex.axis=.65,
                main = 'Time to Next Treatment by BRAF Result',
                names = name.list,
                ylab = "Days",
                xlab="BRAF Result\n(0 = Neg, 1 = Pos)")
      }
      else if(input$groupBy == 'NRAS') {
        plot.data$NRAS <- factor(plot.data$NRAS)
        boxplot(DAYS_TO_NEXT_RX ~ NRAS , plot.data, cex.axis=.65,
                main = 'Time to Next Treatment by NRAS Result',
                names = name.list,
                ylab = "Days",
                xlab="NRAS Result\n(0 = Neg, 1 = Pos)")
      }
    
      else if(input$groupBy == 'DRUG_CLASS') {
        plot.data$DRUG_CLASS <- factor(plot.data$DRUG_CLASS)
        boxplot(DAYS_TO_NEXT_RX ~ DRUG_CLASS , plot.data, cex.axis=.65,
                main = 'Time to Next Treatment by Drug Class',
                names = name.list,
                ylab = "Days",
                xlab="Drug Class")
      }
      else if(input$groupBy == 'MUTATION') {
        plot.data$MUTATION <- factor(plot.data$MUTATION)
        boxplot(DAYS_TO_NEXT_RX ~ MUTATION , plot.data, cex.axis=.65,
                main = 'Time to Next Treatment by Mutation Type',
                names = name.list,
                ylab = "Days",
                xlab="Drug Class")
      }
      else {
        plot.data$MEDICATION <- factor(plot.data$MEDICATION)
        boxplot(DAYS_TO_NEXT_RX ~ MEDICATION , plot.data, cex.axis=.65,
                main = 'Time to Next Treatment by Drug Name',
                names = name.list,
                ylab = "Days",
                xlab="Drug Name")
      }
    
  },height=700)
  
  # Renders boxlplot for time to next treatment
  output$boxResponseMed <- renderPlot({
    plot.data <- selectPats(data, two_class = input$twoClass)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    plot.data$MEDICATION <- factor(plot.data$MEDICATION)
    name.list <- sort(unique(factor(plot.data[,input$groupBy])))
    name.list <- paste(name.list,' (', table(plot.data[,input$groupBy]), ')',sep='')
    boxplot(DAYS_TO_NEXT_RX ~ MEDICATION , plot.data, cex.axis=.65,
            main = 'Time to Next Treatment by Drug Name',
            names = name.list,
            ylab = "Days",
            xlab="Drug Name")
    },height=900)
  
  
  output$boxSummary <- renderPrint({  
    plot.data <- selectPats(data, two_class = input$twoClass)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    plot.data[,input$groupBy] <- factor(plot.data[,input$groupBy])
    
    tapply(plot.data$DAYS_TO_NEXT_RX, plot.data[,input$groupBy], summary)
  })
  
  # Generates Survival Plot  
  output$survCurv <- renderPlot({
    plot.data <- selectPats(data, two_class = input$twoClass)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    
    name.list <- sort(unique(factor(plot.data[,input$groupBy])))
    name.list <- paste(name.list,' (', table(factor(plot.data[,input$groupBy])), ')',sep='')
       
      ml.surv <- coxModel(plot.data, input$groupBy)
      mfit <- survfit(ml.surv)
      plot(mfit, xlab="Days", ylab="Percent Living", main="Kaplain Meier",
           col=1:length(unique(plot.data[,input$groupBy])), conf.int = FALSE)
        legend(x="topright", unique(plot.data[,input$groupBy]),
               leg=name.list,
               fill=1:length(unique(plot.data[,input$groupBy])),
               inset=c(0,0))

  },height=500)
  
  output$survSummary <- renderPrint({  
    plot.data <- selectPats(data, two_class = input$twoClass)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    
    ml.surv <- coxModel(plot.data, input$groupBy)
    mfit <- survfit(ml.surv)
    list(ml.surv, mfit)
  })
  
#  THIS CODE IS NOT YET TO BE DEPLOYED.  ONLY EXPLORING SURVDIFF FUNCTION.  
#   output$survDiffSummary <- renderPrint({
#     plot.data <- selectPats(data)
#     if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
#       plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
#       plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
#     }
# 
#     diff_output <- switch(input$groupBy,
#                           MEDICATION = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ MEDICATION, data=plot.data),
#                           DRUG_CLASS = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ DRUG_CLASS, data=plot.data),
#                           SEX = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ SEX, data=plot.data),
#                           BRAF = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ BRAF, data=plot.data),
#                           MUTATION = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ MUTATION, data=plot.data),
#                           NRAS = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ NRAS, data=plot.data)
#     )
#      diff_output  
#   })

  
  # Displays Age and Sex Distributions of Patients in Cohort
  output$cohortAgeSex <- renderPlot({
    plot.data <- selectPats(data)
    layout(matrix(c(1,2),1,2,byrow=T),widths=c(1,2))
    barplot(table(plot.data$SEX), main = "Cohort Sex Distribution")
    hist(plot.data$AGE, main="Cohort Age Distribution", xlab="Age (Years)", ylab="# Patients",col='grey')
  }, height = 250)
  
  # Displays Mutation Distributions of Patients in Cohort
  output$cohortGen <- renderPlot({
    plot.data <- selectPats(data)
    par(mfrow=c(1,2))
    barplot(table(plot.data$MUTATION), main="Cohort Positive Mutation Results Distribution")
    barplot(table(plot.data$DRUG_CLASS), main="Cohort Treatment Drug Class Distribution",
            names = c("CHEMO.","IMMUNO.","KINASE INH.","ANTIBODY"))
    
  }, height = 250)
  
  # Displays Drugs Used on Patients in Cohort
  output$cohortDrugs <- renderPlot({
    plot.data <- selectPats(data)
   # par(mar=c(5,4,2,2),mgp = c(1, .5, 0))
    par(mar=c(9,4,2,2))
    barplot(table(plot.data$MEDICATION), main="Cohort Drug Distribution",las=3)                              
  }, height = 400)  
  
  # Creates Table of Clinical Data for Display/Download
  output$cohortTableClinical <- renderDataTable({
    plot.data <- selectPats(data)[,c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                                     "DRUG_CLASS", "DAYS_TO_NEXT_RX",
                                     "DAYS_TO_DEATH","R.CENSORED","TUMOR_CHANGE",
                                     "SOURCE")]
    colnames(plot.data)<- c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                            "DRUG_CLASS", "DAYS TO 2ND TREATMENT",
                            "DAYS TO DEATH","ALIVE AT END OF STUDY","TUMOR CHANGE",
                            "INSTITUTION")
    plot.data
  },options=list(iDisplayLength = 25))  
  
  # Creates Table of Drug Data for Display/Download
  output$cohortTableDrugsFull <- renderDataTable({
    plot.data <- selectPats(data)
    MRNS <- unique(plot.data$MRN)
    data <- pat.drugs.full[which(pat.drugs.full$MRN %in% MRNS),]
    data
  },options=list(iDisplayLength = 25))  
  
  # Creates Table of Drug Data for Display/Download
  output$cohortTableDrugs <- renderDataTable({
    plot.data <- selectPats(data)
    MRNS <- unique(plot.data$MRN)
    data <- pat.drugs[which(pat.drugs$MRN %in% MRNS),]
    data
  },options=list(iDisplayLength = 25))
  
  # Dowload Handler for Clinical Data
  output$downloadDataClinical <- downloadHandler(
    filename = function() { paste("mrlu_patient_cohort_clinical", '.csv', sep='') },
    content = function(file) {
      plot.data <- selectPats(data)[,c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                                       "DRUG_CLASS", "DAYS_TO_NEXT_RX",
                                       "DAYS_TO_DEATH","R.CENSORED","TUMOR_CHANGE",
                                       "SOURCE")]
      colnames(plot.data)<- c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                              "DRUG_CLASS", "DAYS TO 2ND TREATMENT",
                              "DAYS TO DEATH","ALIVE AT END OF STUDY","TUMOR CHANGE",
                              "INSTITUTION")
      plot.data
      write.csv(plot.data, file)
    }
  )
  
  # Dowload Handler for Drug Data
  output$downloadDataDrug <- downloadHandler(
    filename = function() { paste("mrlu_patient_cohort_drugs", '.csv', sep='') },
    content = function(file) {
      plot.data <- selectPats(data)
      MRNS <- unique(plot.data$MRN)
      data <- pat.drugs[which(pat.drugs$MRN %in% MRNS),]
      write.csv(data, file)
    }
  )
  
  })
