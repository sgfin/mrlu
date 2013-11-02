
require('survival')

#filepath <- "~/Documents/0 Research/Rubin/0 mrlu/"
filepath <- "/home/samfin/ShinyApps/mrlu/"

data <- read.csv(paste(filepath, "clinical_data.csv",sep=""))
pat.drugs <- read.csv(paste(filepath, "pat_drugs.csv", sep=""))

extract_drugs <- function(MRN, number.treats) {
  drug.list <- pat.drugs$NAME[which(pat.drugs$MRN == MRN)]
  drug.list <- as.character(drug.list)
  
  if("PACLITAXEL" %in% drug.list && "CARBOPLATIN" %in% drug.list){
    drug.list <- gsub("PACLITAXEL", "PACL/CARBO", x=drug.list)
    drug.list <- gsub("CARBOPLATIN", "PACL/CARBO",x=drug.list)
  }
  drug.list <- paste(drug.list[1:min(number.treats,length(drug.list))], collapse=", ")
}

extract_drug_class <- function(MRN) {
  d_class <- as.character(pat.drugs$DRUG_CLASS[which(pat.drugs$MRN == MRN)])
  d_class <- d_class[!is.na(d_class)]
  return(d_class[1])
}

MEDICATION <- sapply(data$MRN, extract_drugs, 1)
DRUG_CLASS <- sapply(data$MRN, extract_drug_class) # REPLACE WITH SOME MAP
TUMOR_CHANGE <- runif(length(data$MRN),-1,1)

data <- cbind(data, MEDICATION, DRUG_CLASS, TUMOR_CHANGE)


shinyServer(function(input, output) {
    
#   output$inputPat <- renderUI({
#     textInput("POI.ID", "Enter Patient ID:", value="MR0375" )
#   })
  
  output$selectSex <- renderUI({
# #     if(input$specifyPat == 'import'){
# #       return(selectInput("sex", "Sex:",
# #                          list("MALE" = "MALE", 
# #                               "FEMALE" = "FEMALE"),
# #                          selected = as.character(ehr.data[input$POI.ID,]$sex)))
# #     }
# #     else{
      selectInput("sex", "Sex:",
                  list("MALE" = "MALE", 
                       "FEMALE" = "FEMALE"))
# #     }
  })
#   
  output$selectRace <- renderUI({
# #     if(input$specifyPat == 'import') {
# #       return(selectInput("race", "Race:",
# #                          list("White" = 'WHITE', 
# #                               "Black or African American" = "BLACK OR AFRICAN AMERICAN",
# #                               "Asian" = "ASIAN"),
# #                          selected = as.character(ehr.data[input$POI.ID,]$race)))
# #     }
    selectInput("race", "Race:",
                list("White" = 'WHITE', 
                     "Black or African American" = "BLACK OR AFRICAN AMERICAN",
                     "Asian" = "ASIAN"
                ))
  })
#   
  output$selectAge <- renderUI({
# #     if(input$specifyPat == 'import'){
# #       return(numericInput("POI.age", "Patient Age:", as.numeric(ehr.data[input$POI.ID,]$age)))
# #     }
    numericInput("POI.age", "Patient Age:", median(data$AGE))
  })
  
  output$selectBRAF <- renderUI({
# #     if(input$sp == 'import'){
# #       return(selectInput("braf_status", "BRAF STATUS:",
# #                          list("POSITIVE" = "POSITIVE", 
# #                               "NEGATIVE" = "NEGATIVE"),
# #                          selected = as.character(ehr.data[input$POI.ID,]$sex)))
# #     }
# #     else{
    selectInput("braf", "BRAF Test Result:",
                list("POSITIVE" = "1", 
                     "NEGATIVE" = "0"))
# #    }
  })
  
  output$selectNRAS <- renderUI({
# #     if(input$specifyPat == 'import'){
# #       return(selectInput("sex", "Sex:",
# #                          list("MALE" = "MALE", 
# #                               "FEMALE" = "FEMALE"),
# #                          selected = as.character(ehr.data[input$POI.ID,]$sex)))
# #     }
# #     else{
      selectInput("nras", "NRAS Test Result:",
                  list("POSITIVE" = "1", 
                       "NEGATIVE" = "0"))
# #     }
  })
  
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
  
  output$drug_classes <- renderUI({   
    checkboxGroupInput("drug_class", "Include:", sort(unique(data$DRUG_CLASS)), sort(unique(data$DRUG_CLASS)))
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
  
#  Function to Narrow Down Dataset According to UI inputs
  selectPats <- function(data, subgroups = TRUE){
    plot.data <- data
    
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
    plot.data
  }
  
  # Build Cox Model
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

  plotHeight <- function(){
    plot.data <- selectPats(data)
    distinct = unique(plot.data[,input$groupBy])
    return(length(distinct)*200+100)
  }
  
  nDistinct <-3
  # Generates Survival Plot  
  output$barResponse <- renderPlot({
    plot.data <- selectPats(data)
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
    plot.data <- selectPats(data)
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
  
  output$boxResponseMed <- renderPlot({
    plot.data <- selectPats(data)
    plot.data$MEDICATION <- factor(plot.data$MEDICATION)
    name.list <- sort(unique(factor(plot.data[,input$groupBy])))
    name.list <- paste(name.list,' (', table(plot.data[,input$groupBy]), ')',sep='')
    boxplot(DAYS_TO_NEXT_RX ~ MEDICATION , plot.data, cex.axis=.65,
            main = 'Time to Next Treatment by Drug Name',
            names = name.list,
            ylab = "Days",
            xlab="Drug Name")
    },height=900)
  
  # Generates Survival Plot  
  output$survCurv <- renderPlot({
    plot.data <- selectPats(data)
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
  
  output$cohortAges <- renderPlot({
    plot.data <- selectPats(data)
    layout(matrix(c(1,2),1,2,byrow=T),widths=c(1,2))
    barplot(table(plot.data$SEX), main = "Cohort Sex Distribution")
    hist(plot.data$AGE, main="Cohort Age Distribution", xlab="Age (Years)", ylab="# Patients",col='grey')
  }, height = 250)
  
  output$cohortSexGen <- renderPlot({
    plot.data <- selectPats(data)
    par(mfrow=c(1,2))
    #barplot(table(plot.data$SEX), main = "Cohort Sex Distribution")
    barplot(table(plot.data$MUTATION), main="Cohort Positive Mutation Results Distribution")
    barplot(table(plot.data$DRUG_CLASS), main="Cohort Treatment Drug Class Distribution",
            names = c("CHEMO.","IMMUNO.","KINASE INH.","ANTIBODY"))
    
  }, height = 250)
  
  output$cohortDrugs <- renderPlot({
    plot.data <- selectPats(data)
   # par(mar=c(5,4,2,2),mgp = c(1, .5, 0))
    par(mar=c(9,4,2,2))
    barplot(table(plot.data$MEDICATION), main="Cohort Drug Distribution",las=3)                              
  }, height = 400)  
  
  deidentifyStan <- function(MRNS){
    MRNS_deidentified <- MRNS
    MRNS_list <- unique(MRNS)
    for(i in 1:length(MRNS_list)){
      MRNS_deidentified[which(MRNS_deidentified==MRNS_list[i])] <- paste("STAN",i)
    }
    MRNS_deidentified
  }
  
  output$cohortTableClinical <- renderDataTable({
    plot.data <- selectPats(data)[,c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                                     "DRUG_CLASS", "DAYS_TO_NEXT_RX",
                                     "DAYS_TO_DEATH","R.CENSORED","TUMOR_CHANGE",
                                     "SOURCE")]
    stanford.inds <- which(plot.data$SOURCE == "STANFORD")
    plot.data$MRN[stanford.inds] <- deidentifyStan(plot.data$MRN[stanford.inds])
    vand.inds <- which(plot.data$SOURCE == "VANDERBILT")
    plot.data$MRN[vand.inds] <- paste("VAND", plot.data$MRN[vand.inds])
    colnames(plot.data)<- c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                            "DRUG_CLASS", "DAYS TO 2ND TREATMENT",
                            "DAYS TO DEATH","ALIVE AT END OF STUDY","TUMOR CHANGE",
                            "INSTITUTION")
    plot.data
  },options=list(iDisplayLength = 25))  
  
  
  output$cohortTableDrugs <- renderDataTable({
    plot.data <- selectPats(data)
    MRNS <- unique(plot.data$MRN)
    data <- pat.drugs[which(pat.drugs$MRN %in% MRNS),]
    stanford.inds <- which(data$INSTITUTION != "VANDERBILT")
    data$MRN[stanford.inds] <- deidentifyStan(data$MRN[stanford.inds])
    data$MRN[which(data$INSTITUTION == "VANDERBILT")] <- paste("VAND", data$MRN[which(data$INSTITUTION == "VANDERBILT")])    
    data
  },options=list(iDisplayLength = 25))  
  
  output$downloadDataClinical <- downloadHandler(
    filename = function() { paste("mrlu_patient_cohort_clinical", '.csv', sep='') },
    content = function(file) {
      plot.data <- selectPats(data)[,c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                                       "DRUG_CLASS", "DAYS_TO_NEXT_RX",
                                       "DAYS_TO_DEATH","R.CENSORED","TUMOR_CHANGE",
                                       "SOURCE")]
      stanford.inds <- which(plot.data$SOURCE == "STANFORD")
      plot.data$MRN[stanford.inds] <- deidentifyStan(plot.data$MRN[stanford.inds])
      vand.inds <- which(plot.data$SOURCE == "VANDERBILT")
      plot.data$MRN[vand.inds] <- paste("VAND", plot.data$MRN[vand.inds])
      colnames(plot.data)<- c("MRN","SEX","AGE","MUTATION", "MEDICATION",
                              "DRUG_CLASS", "DAYS TO 2ND TREATMENT",
                              "DAYS TO DEATH","ALIVE AT END OF STUDY","TUMOR CHANGE",
                              "INSTITUTION")
      plot.data
      write.csv(plot.data, file)
    }
  )
  
  output$downloadDataDrug <- downloadHandler(
    filename = function() { paste("mrlu_patient_cohort_drugs", '.csv', sep='') },
    content = function(file) {
      plot.data <- selectPats(data)
      MRNS <- unique(plot.data$MRN)
      data <- pat.drugs[which(pat.drugs$MRN %in% MRNS),]
      stanford.inds <- which(data$SOURCE == "STANFORD")
      data$MRN[stanford.inds] <- deidentifyStan(data$MRN[stanford.inds])
      data$MRN[which(data$SOURCE == "VANDERBILT")] <- paste("VAND", data$MRN[which(plot.data$SOURCE == "VANDERBILT")])
      write.csv(data, file)
    }
  )
  
  })