require('survival')
library("ggplot2")


source("multiplot.R")

connect.db <- function(username, db.name, host = "127.0.0.1") {
  require(DBI)
  require(RMySQL)
  m <- dbDriver("MySQL")
  
  con <- dbConnect(m, username=username, dbname = db.name, host=host)
  con
}

.libPaths("/home/samfin/R/x86_64-redhat-linux-gnu-library/3.0")
con = connect.db(username = "", db.name="test")

data <- dbGetQuery(con, "SELECT * from clinical_data")
data$MRN <- factor(data$MRN); data$SEX <- factor(data$SEX); data$RX1 <- factor(data$RX1);
data$SOURCE <- factor(data$SOURCE); data$MUTATION <- factor(data$MUTATION); data$R.CENSORED <- data$R.CENSORED=="T"

pat.drugs <- dbGetQuery(con, "select * from drug_data group by MRN, NAME, DRUG_CLASS, INSTITUTION
order by MRN, STR_TO_DATE(ORDERING_DATE, '%m/%d/%Y'), NAME")
pat.drugs.full <- dbGetQuery(con, "SELECT * from drug_data order by MRN, STR_TO_DATE(ORDERING_DATE, '%m/%d/%Y'), NAME")

#filepath <- "~/Documents/0 Research/Rubin/0 mrlu/"
filepath <- "/home/samfin/ShinyApps/mrlu/"

#data <- read.csv(paste(filepath, "clinical_data.csv",sep=""))
#pat.drugs <- read.csv(paste(filepath, "pat_drugs.csv", sep=""))
#pat.drugs.full <- read.csv(paste(filepath, "drug_data_full.csv", sep=""))

# Function to identify the first "number.treats" drugs used on patient 'MRN'
# Combines PACLI and Carbo when used together
extract_drugs <- function(MRN, number.treats) {
  drug.list <- pat.drugs[which(pat.drugs$MRN == MRN), c("NAME", "ORDERING_DATE")]
  
  drug.list$NAME <- gsub("/", "/\n", x=drug.list[,"NAME"])
  
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
MEDICATION_TWO <- sapply(data$MRN, extract_drugs, 2)
DRUG_CLASS <- sapply(data$MRN, extract_drug_class, 1) # REPLACE WITH SOME MAP
DRUG_CLASS_TWO <- sapply(data$MRN, extract_drug_class, 2) # REPLACE WITH SOME MAP
TUMOR_CHANGE <- runif(length(data$MRN),-1,1)

TTNT.CENSORED <- data$DAYS_TO_NEXT_RX == 0
data$DAYS_TO_NEXT_RX[TTNT.CENSORED] <- data$DAYS_TO_DEATH[TTNT.CENSORED]
data <- cbind(data, TTNT.CENSORED)

# Finalizes Data
data <- cbind(data, MEDICATION, DRUG_CLASS, TUMOR_CHANGE)
data2 <- data
data2$DRUG_CLASS <- DRUG_CLASS_TWO
data2$MEDICATION <- MEDICATION_TWO

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
      if(input$twoDrug){
        return(checkboxGroupInput("drug_name", "Include:", sort(unique(data2$MEDICATION))))
      }
      else{
        return(checkboxGroupInput("drug_name", "Include:", sort(unique(data$MEDICATION))))
      }
    }
    else {
      if(input$twoDrug){
        return(checkboxGroupInput("drug_name", "Include:", sort(unique(data2$MEDICATION)), sort(unique(data2$MEDICATION))))
      }
      else{
        return(checkboxGroupInput("drug_name", "Include:", sort(unique(data$MEDICATION)), sort(unique(data$MEDICATION))))
      }
    }
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
  selectPats <- function(data, subgroups = TRUE, two_class = input$twoClass, two_drug = input$twoDrug){
    if(two_class || two_drug){
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
    if(input$outcomeVar=="TTNT"){
      switch(groupBy,
             MEDICATION = coxph(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ 1 + strata(MEDICATION), data=plot.data),
             DRUG_CLASS = coxph(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ 1 + strata(DRUG_CLASS), data=plot.data),
             SEX = coxph(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ 1 + strata(SEX), data=plot.data),
             BRAF = coxph(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ 1 + strata(BRAF), data=plot.data),
             MUTATION = coxph(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ 1 + strata(MUTATION), data=plot.data),
             NRAS = coxph(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ 1 + strata(NRAS), data=plot.data)
      )
    }
    else{
      switch(groupBy,
             MEDICATION = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(MEDICATION), data=plot.data),
             DRUG_CLASS = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(DRUG_CLASS), data=plot.data),
             SEX = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(SEX), data=plot.data),
             BRAF = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(BRAF), data=plot.data),
             MUTATION = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(MUTATION), data=plot.data),
             NRAS = coxph(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ 1 + strata(NRAS), data=plot.data)
      )
    }
  }

  # Calculates desired height of plot
  plotHeight <- function(){
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
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
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    plot.data <- subset(plot.data, TTNT.CENSORED==F)
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
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    plot.data$MEDICATION <- factor(plot.data$MEDICATION)
    plot.data <- subset(plot.data, TTNT.CENSORED==F)
    name.list <- sort(unique(factor(plot.data[,input$groupBy])))
    name.list <- paste(name.list,' (', table(plot.data[,input$groupBy]), ')',sep='')
    boxplot(DAYS_TO_NEXT_RX ~ MEDICATION , plot.data, cex.axis=.65,
            main = 'Time to Next Treatment by Drug Name',
            names = name.list,
            ylab = "Days",
            xlab="Drug Name")
    },height=900)
  
  output$boxSummary <- renderTable({
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    plot.data[,input$groupBy] <- factor(plot.data[,input$groupBy])
    plot.data <- subset(plot.data, TTNT.CENSORED==F)
    tab <- tapply(plot.data$DAYS_TO_NEXT_RX, plot.data[,input$groupBy], summary)
    nms.cols <- names(tab[[1]])
    nms.rows <- names(tab)
    tab <- matrix(unlist(tab), ncol=length(nms.cols), byrow=T)
    colnames(tab) <- nms.cols
    rownames(tab) <- nms.rows
    tab
  })
  
  output$t2ntANOVA <- renderPrint({
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    plot.data <- subset(plot.data, TTNT.CENSORED==F)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    plot.data <- subset(plot.data, plot.data$DAYS_TO_NEXT_RX != 0)
    if(input$groupBy == 'SEX') {
      summary(aov(DAYS_TO_NEXT_RX ~ SEX, plot.data))
    }
    else if(input$groupBy == 'BRAF') {
      summary(aov(DAYS_TO_NEXT_RX ~ BRAF, plot.data))
    }
    else if(input$groupBy == 'NRAS') {
      summary(aov(DAYS_TO_NEXT_RX ~ NRAS, plot.data))
    }
    else if(input$groupBy == 'DRUG_CLASS') {
      summary(aov(DAYS_TO_NEXT_RX ~ DRUG_CLASS, plot.data))
    }
    else if(input$groupBy == 'MUTATION') {
      summary(aov(DAYS_TO_NEXT_RX ~ MUTATION, plot.data))
    }
    else {
      summary(aov(DAYS_TO_NEXT_RX ~ MEDICATION, plot.data))
    }
  })
  
  output$t2ntPairwiseTT <- renderTable({
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$DAYS_TO_NEXT_RX)),]
    plot.data <- subset(plot.data, TTNT.CENSORED==F)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    plot.data <- subset(plot.data, plot.data$DAYS_TO_NEXT_RX != 0)
    pairwise.t.test(plot.data[,"DAYS_TO_NEXT_RX"], plot.data[,input$groupBy])$p.value
  })
  
  output$tumorBurdenSummary <- renderPrint({
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$TUMOR_CHANGE)),]
    plot.data[,input$groupBy] <- factor(plot.data[,input$groupBy])
    tapply(plot.data$TUMOR_CHANGE, plot.data[,input$groupBy], summary)
  })
  
  output$tumorBurdenANOVA <- renderPrint({
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$TUMOR_CHANGE)),]
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    if(input$groupBy == 'SEX') {
      summary(aov(TUMOR_CHANGE ~ SEX, plot.data))
    }
    else if(input$groupBy == 'BRAF') {
      summary(aov(TUMOR_CHANGE ~ BRAF, plot.data))
    }
    else if(input$groupBy == 'NRAS') {
      summary(aov(TUMOR_CHANGE ~ NRAS, plot.data))
    }
    else if(input$groupBy == 'DRUG_CLASS') {
      summary(aov(TUMOR_CHANGE ~ DRUG_CLASS, plot.data))
    }
    else if(input$groupBy == 'MUTATION') {
      summary(aov(TUMOR_CHANGE ~ MUTATION, plot.data))
    }
    else {
      summary(aov(TUMOR_CHANGE ~ MEDICATION, plot.data))
    }
  })
  
  output$tumorBurdenPairwiseTT <- renderPrint({
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    plot.data <- plot.data[which(!is.na(plot.data$TUMOR_CHANGE)),]
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    pairwise.t.test(plot.data[,"TUMOR_CHANGE"], plot.data[,input$groupBy])
  })

  # Generates Survival Plot  
  output$survCurv <- renderPlot({
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    
    name.list <- sort(unique(factor(plot.data[,input$groupBy])))
    name.list <- paste(name.list,' (', table(factor(plot.data[,input$groupBy])), ')',sep='')
       
      ml.surv <- coxModel(plot.data, input$groupBy)
      mfit <- survfit(ml.surv)
      plot(mfit, xlab="Days", ylab="Survival", main="Kaplain Meier",
           col=1:length(unique(plot.data[,input$groupBy])), conf.int = FALSE)
        legend(x="topright", unique(plot.data[,input$groupBy]),
               leg=name.list,
               fill=1:length(unique(plot.data[,input$groupBy])),
               inset=c(0,0))

  },height=500)
  
  output$survCall <- renderPrint({  
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    coxModel(plot.data, input$groupBy)$call
  })
  
  output$survSummary <- renderTable({  
    plot.data <- selectPats(data, two_class = input$twoClass, two_drug = input$twoDrug)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }
    summary(survfit(coxModel(plot.data, input$groupBy)))$table
  })
  
  output$survDiffSummary <- renderTable({
    plot.data <- selectPats(data)
    if(input$groupBy=='NRAS' || input$groupBy=='BRAF'){
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='0')] <- paste(input$groupBy, " Negative", sep="")
      plot.data[,input$groupBy][which(plot.data[,input$groupBy]=='1')] <- paste(input$groupBy, " Positive",sep="")
    }

    if(input$outcomeVar == "TTNT"){
      diff_output <- switch(input$groupBy,
                            MEDICATION = survdiff(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ MEDICATION, data=plot.data),
                            DRUG_CLASS = survdiff(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ DRUG_CLASS, data=plot.data),
                            SEX = survdiff(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ SEX, data=plot.data),
                            BRAF = survdiff(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ BRAF, data=plot.data),
                            MUTATION = survdiff(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ MUTATION, data=plot.data),
                            NRAS = survdiff(Surv(DAYS_TO_NEXT_RX, event=(!plot.data$TTNT.CENSORED), type="right") ~ NRAS, data=plot.data)
      )
    }
    else{
      diff_output <- switch(input$groupBy,
                            MEDICATION = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ MEDICATION, data=plot.data),
                            DRUG_CLASS = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ DRUG_CLASS, data=plot.data),
                            SEX = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ SEX, data=plot.data),
                            BRAF = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ BRAF, data=plot.data),
                            MUTATION = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ MUTATION, data=plot.data),
                            NRAS = survdiff(Surv(DAYS_TO_DEATH, event=(!plot.data$R.CENSORED), type="right") ~ NRAS, data=plot.data)
      )
    }
    
    out <- t(data.frame(paste(deparse(diff_output$call), collapse=""),
    diff_output$chisq,
    (1 - pchisq(diff_output$chisq, length(diff_output$n) - 1))))
    colnames(out) <- c("")
    rownames(out) <- c("R Call", "Chi Squared Statistic", "P-Value")
    out
  }, include.colnames=FALSE)

  
  # Displays Age and Sex Distributions of Patients in Cohort
  output$cohortAgeSex <- renderPlot({
    plot.data <- selectPats(data)
    p.sex <- qplot(SEX, data = plot.data, fill = SOURCE, main = "Cohort Sex Distribution") +
      scale_x_discrete(drop = F) + theme(axis.text=element_text(size=16),
                                         title=element_text(size=12),
                                         axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=16))
    p.age <- qplot(AGE, data = plot.data, fill = SOURCE, main="Cohort Age Distribution") +
      theme(axis.text=element_text(size=16),
            title=element_text(size=12))
    multiplot(p.sex, p.age, layout = matrix(c(1,2,2),1,3,byrow=T))
    #print(p.sex)
  }, height = 250)
  
  # Displays Mutation Distributions of Patients in Cohort
  output$cohortGen <- renderPlot({
    plot.data <- selectPats(data)
    #par(mfrow=c(1,2))
    #barplot(table(plot.data$MUTATION), main="Cohort Positive Mutation Results Distribution")
    #barplot(table(plot.data$DRUG_CLASS), main="Cohort Treatment Drug Class Distribution",
    #        names = c("CHEMO-\nTHERAPY","IMMUNO-\nTHERAPY","KINASE\nINHIBITOR","THERAPEUT.\nANTIBODY"))
    p.mut <- qplot(MUTATION, data=plot.data, fill = SOURCE,
                   main="Cohort Positive Mutation Results Distribution") +
      scale_x_discrete(drop = F) + theme(axis.text=element_text(size=16),
                                         title=element_text(size=12),
                                         axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=16))
    p.drug <- qplot(DRUG_CLASS, data =plot.data, fill = SOURCE,
                    main="Cohort Treatment Drug Class Distribution") +
                    scale_x_discrete(labels = c("CHEMO","IMMUNO","KINASE","ANTIBOD"), drop=F) +
      theme(axis.text=element_text(size=16),
            title=element_text(size=12),
            axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=16))
    multiplot(p.mut, p.drug, layout = matrix(c(1,2),1,2,byrow=T))
  }, height = 250)
  
  # Displays Drugs Used on Patients in Cohort
  output$cohortDrugs <- renderPlot({
    plot.data <- selectPats(data) 
    p.med <- qplot(MEDICATION, data=plot.data, fill = SOURCE,
          main="Cohort Drug Distribution", drop=F) +
      theme(axis.text.y= element_text(size=16),
            axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=16),
            title=element_text(size=12)
            )
    print(p.med)
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
  },options=list(iDisplayLength = 10))  
  
  # Creates Table of Drug Data for Display/Download
  output$cohortTableDrugsFull <- renderDataTable({
    plot.data <- selectPats(data)
    MRNS <- unique(plot.data$MRN)
    data <- pat.drugs.full[which(pat.drugs.full$MRN %in% MRNS),]
    data
  },options=list(iDisplayLength = 10))  
  
  # Creates Table of Drug Data for Display/Download
  output$cohortTableDrugs <- renderDataTable({
    plot.data <- selectPats(data)
    MRNS <- unique(plot.data$MRN)
    data <- pat.drugs[which(pat.drugs$MRN %in% MRNS),]
    data
  },options=list(iDisplayLength = 10))
  
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
