testRun001 <- function() {
  #SampleID layout: [YYYYMMDD]-[INITIALS]-[SAMPLENAME]-TP[##]-R[##]
  data <- readTxt(dir=system.file("extdata", "AD001", package = "vicellR"), extended=T, purgeDate=T, sampleIDSearch=c("_unit1","_unit2"), sampleIDReplace=c("-R01","-R02"))
  df["replicate"] <- 1
  groups <- regmatches(df$sampleName, regexpr("^.{2,3}-", df$sampleName, perl=T))
  groups <- substr(groups, 1, nchar(groups)-1)
  df$group1 <- groups
  groups <- gsub("^.{2,3}-", "", df$sampleName, perl=T)
  group2 <- regmatches(groups, regexpr("^.{2,5}-", groups, perl=T))
  group2 <- substr(group2, 1, nchar(group2)-1)
  df$group2 <- group2
  df$group3 <- gsub("^.{2,3}-.{2,5}-", "", df$sampleName, perl=T)
  df <- calculateGrowth(df)
  drawConcVia(df, save=T, height=210, width=390)
  drawGrowth(df,save=T,height=210,width)
}

testRun002 <- function() {
  #SampleID layout: [YYYYMMDD]-[INITIALS]-[SAMPLENAME]-TP[##]-R[##]
  data <- readTxt(dir=system.file("extdata", "1710_Reactor", package = "vicellR"), extended=T, purgeDate=T, sampleIDSearch=c("_unit1","_unit2"), sampleIDReplace=c("-R01","-R02"))

  #ID-Adjustment-----------------------------------------------------------------------
  timepoints1<-sapply(data,'[[',"runDate") [grepl("R01",sapply(data,'[[',"sampleName"))]
  timepoints2<-sapply(data,'[[',"runDate") [grepl("R02",sapply(data,'[[',"sampleName"))]
  timepoints1<-timepoints1 [order(timepoints1)]
  timepoints2<-timepoints2 [order(timepoints2)]
  for(j in 1:length(unique(timepoints1))){
    data[[which(sapply(data,'[[',"runDate") %in% timepoints1 [j])]] [["sampleName"]]<-gsub("_[0-9][0-9]-",paste("-TP",sprintf("%02d",j),"-",sep=""),data[[which(sapply(data,'[[',"runDate") %in% timepoints1 [j])]] [["sampleName"]])
  }
  for(j in 1:length(unique(timepoints2))){
    data[[which(sapply(data,'[[',"runDate") %in% timepoints2 [j])]] [["sampleName"]]<-gsub("_[0-9][0-9]-",paste("-TP",sprintf("%02d",j),"-",sep=""),data[[which(sapply(data,'[[',"runDate") %in% timepoints2 [j])]] [["sampleName"]])
  }
  #-------------------------------------------------------------------------------------

  data<-avgViabDiam(data)   #Calculates and adds viable cell diameter
  data<-dataTransformation(data) #Calculates and adds viable cell volume (VCV) and ln-transformed VCD and VCV
  data<-assignReplicates(data,"Reactor-output") #Assignes replicates
  data<-addSd(data) #Adds STDEV for different factors
  df<-getDF(data) #Creates Dataframe out of list
  df<-drawBatchGraphs(df,path="Reactor-output/batch") #Draws vs-t-plots
  exportCSV(df,path="Reactor-output/batch")
}

testRun003 <- function() {
  #SampleID layout: [YYYYMMDD]-[INITIALS]-[SAMPLENAME]-TP[##]-R[##]
  #ReplicateData
  data <- readTxt(dir=system.file("extdata", "replicateData", package = "vicellR"), extended=T, purgeDate=F, sampleIDSearch=c(" D0-1"," D1-1"," D1-2"," D2-1"," D2-2"," D3-1"," D3-2"," D4-1"," D4-2"," D5-1"," D6-1"," D7-1"," D7-2"," D8-1"," D9-1","-1","-2","GK"), sampleIDReplace=c("-TP01","-TP02","-TP03","-TP04","-TP05","-TP06","-TP07","-TP08","-TP09","-TP10","-TP11","-TP12","-TP13","-TP14","-TP15","-R01","-R02",""))

  #ID-Adjustment--------------------------------------------------------------------
  for(i in 1:length(data)){
      data[[i]] [["sampleName"]]<-gsub("^[0-9]+[-]","",data[[i]] [["sampleName"]])
  }
  #---------------------------------------------------------------------------------

  data<-avgViabDiam(data) #Calculates and adds viable cell diameter
  data<-dataTransformation(data) #Calculates and adds viable cell volume (VCV) and ln-transformed VCD and VCV
  data<-assignReplicates(data,"Rep-output") #Assignes replicates
  data<-addSd(data) #Adds STDEV for different factors
  df<-getDF(data) #Creates Dataframe out of list
  df<-drawBatchGraphs(df,path="Rep-output/batch") #Draws vs-t-plots
  corrcombgrowthRates<-calculateGrowthCorr(df,comb=TRUE,corrSeed=TRUE) #Calculate growthrates, replicates combined and seed-corrected
  combgrowthRates<-calculateGrowthCorr(df,comb=TRUE) #Calculate growthrates, replicates combined, not seed-corrected
  growthRates<-calculateGrowthCorr(df,comb=FALSE) #Calculate growthrates, replicates separated
  outputVariable(variable=combgrowthRates,parameter="growthRates",path="Rep-output/batch/growth/comb/")  #writes growth-data into file
  outputVariable(variable=corrcombgrowthRates,parameter="growthRates",path="Rep-output/batch/growth/comb/seed-corrected/")
  outputVariable(variable=growthRates,parameter="growthRates",path="Rep-output/batch/growth/sep/")  #writes growth-data into file
  drawGrowthCharts(combgrowthRates,path="Rep-output/batch/growth/comb") #Draws growth charts
  drawGrowthCharts(corrcombgrowthRates,path="Rep-output/batch/growth/comb/seed-corrected") #Draws growth charts
  drawGrowthCharts(growthRates,path="Rep-output/batch/growth/sep") #Draws growth charts
  exportCSV(df,path="Rep-output/batch")
}

testRun004 <- function() {
  #SampleID layout: [YYYYMMDD]-[INITIALS]-[SAMPLENAME]-TP[##]-R[##]
  #NN
  data <- readTxt(dir=system.file("extdata", "NN", package = "vicellR"), extended=T, purgeDate=F, sampleIDSearch=c("_D09","_D0","_D1","_D2","_D3","_D4","_D5","_D6","_D7","_D8","_D9","_D10"), sampleIDReplace=c("-TP10","-TP01","-TP02","-TP03","-TP04","-TP05","-TP06","-TP07","-TP08","-TP09","-TP10","-TP11"))
  data<-avgViabDiam(data) #Calculates and adds viable cell diameter
  data<-dataTransformation(data) #Calculates and adds viable cell volume (VCV) and ln-transformed VCD and VCV
  data<-assignReplicates(data,"NN-output") #Assignes replicates
  data<-addSd(data) #Adds STDEV for different factors
  df<-getDF(data) #Creates Dataframe out of list
  grouped_df<-df
  #Group EmptyVector and 29ABneg together
  for(i in 1:length(unique(grouped_df[grep("EmptyVector",rownames(grouped_df)),"replicateID"]))){
    grouped_df[grouped_df[,"replicateID"]==unique(grouped_df[grep("EmptyVector",rownames(grouped_df)) ,"replicateID"]) [i],"replicateNumber"]<-i
  }

  for(i in 1:length(unique(grouped_df[grep("29ABneg",rownames(grouped_df)),"replicateID"]))){
    grouped_df[grouped_df[,"replicateID"]==unique(grouped_df[grep("29ABneg",rownames(grouped_df)),"replicateID"]) [i],"replicateNumber"]<-i
  }
  grouped_df[grep("EmptyVector",rownames(grouped_df)),"replicateID"]<-1
  grouped_df[grep("29ABneg",rownames(grouped_df)),"replicateID"]<-2
  grouped_df[grep("EmptyVector",rownames(grouped_df)),"replicateName"]<-"EmptyVector"
  grouped_df[grep("29ABneg",rownames(grouped_df)),"replicateName"]<-"29ABneg"
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------
  df<-drawBatchGraphs(df,path="NN-output/batch",perRow=6) #Draws vs-t-plots
  growthRates<-calculateGrowthCorr(df,comb=FALSE)
  outputVariable(variable=growthRates,parameter="growthRates",path="NN-output/batch/growth/sep/")  #writes growth-data into file
  drawGrowthCharts(growthRates,path="NN-output/batch/growth/sep") #Draws growth charts
  exportCSV(df,path="NN-output/batch")
  #Grouped----------------------------------------------------------------------------------------------------------------------------------------
  grouped_df<-drawBatchGraphs(grouped_df,path="NN-output/grouped") #Draws vs-t-plots
  grouped_growthRates<-calculateGrowthCorr(grouped_df,comb=TRUE,corrSeed = TRUE)
  outputVariable(variable=grouped_growthRates,parameter="growthRates",path="NN-output/grouped/growth")  #writes growth-data into file
  drawGrowthCharts(grouped_growthRates,path="NN-output/grouped/growth") #Draws growth charts
  exportCSV(grouped_df,path="NN-output/batch")
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------


}

testRun005 <- function() {
  #SampleID layout: [YYYYMMDD]-[INITIALS]-[SAMPLENAME]-TP[##]-R[##]
  #ReplicateData
  data <- readTxt(dir=system.file("extdata", "CS", package = "vicellR"), extended=T, purgeDate=TRUE,sampleIDSearch = c("_","-CS-","_CS_","-GK-","MOK","PRDM_","PRDM-","CHOM","-D0","-D1","-D2","-D3","-D4","-D5","-D6","_D5","MOCK1","PRDM2-R1","D3E7","R1","R2"),sampleIDReplace=c("-","-","-","-","MOCK","PRDM1-","PRDM1-","M","-TP01","-TP02","-TP03","-TP04","-TP05","-TP06","-TP07","-TP06","MOCK","PRDM1-R02","D1E7","R01","R02"))
  #ID-Adjustment--------------------------------------------------------------------
  data[["20170820-C1D3-MOCK-R2-CS-D5"]] [["sampleName"]]<-"C1D3-MOCK-R02-TP05" # Wrong day in File name
  #---------------------------------------------------------------------------------

  data<-avgViabDiam(data) #Calculates and adds viable cell diameter
  data<-dataTransformation(data) #Calculates and adds viable cell volume (VCV) and ln-transformed VCD and VCV
  data<-assignReplicates(data,"CS-output") #Assignes replicates
  data<-addSd(data) #Adds STDEV for different factors
  df<-getDF(data) #Creates Dataframe out of list
  df<-drawBatchGraphs(df,path="CS-output/batch") #Draws vs-t-plots
  corrcombgrowthRates<-calculateGrowthCorr(df,comb=TRUE,corrSeed=TRUE) #Calculate growthrates, replicates combined and seed-corrected
  combgrowthRates<-calculateGrowthCorr(df,comb=TRUE,corrSeed=FALSE) #Calculate growthrates, replicates combined and seed-corrected
  growthRates<-calculateGrowthCorr(df,comb=FALSE,lagDelay=0) #Calculate growthrates, replicates separated
  outputVariable(variable=corrcombgrowthRates,parameter="growthRates",path="CS-output/batch/growth/comb/")
  outputVariable(variable=growthRates,parameter="growthRates",path="CS-output/batch/growth/sep/")  #writes growth-data into file
  drawGrowthCharts(corrcombgrowthRates,path="CS-output/batch/growth/comb/") #Draws growth charts
  drawGrowthCharts(growthRates,path="CS-output/batch/growth/sep") #Draws growth charts
  df<-assignTiter(df=df,inputFile=system.file("extdata", "CS","Titer_input_R.txt", package = "vicellR"),outputPath="CS-output/batch/prod" )
  df<-calculateCumulative(df=df,outputPath="CS-output/batch")
  specProd<-calculateProdCorr(df,nmin=5,corrSeed=TRUE,comb=FALSE) #nmin=5, because there are no values at TP1
  combspecProd<-calculateProdCorr(df,nmin=5,corrSeed=TRUE,comb=TRUE)
  drawProdCharts(df,specProd,path="CS-output/batch/prod/sep")
  drawProdCharts(df,combspecProd,path="CS-output/batch/prod/comb")
  outputVariable(variable=specProd,parameter="productivity",path="CS-output/batch/prod/sep/")
  outputVariable(variable=combspecProd,parameter="productivity",path="CS-output/batch/prod/comb/")
  exportCSV(df,path="CS-output/batch")
}

testRun006 <- function() {
  data <- readTxt(dir = system.file("extdata", "NiM", package = "vicellR"),  extended = T, sampleIDSearch = c("dil","_CHOEpoFcst6_b_st","-CHOEpoFc_b_","_CHOEpoFcbatch_","_CHOEpoFcst6_b","_CHEpoFcst6_b_","_CHOEpoFc_b_","_CHOEpoFc_batch_","_CHOEpoFcst6_b_","_CHOEpoFc_st6_b_st","NIM","_1","_2"), sampleIDReplace = c("","-st","-st","-","-st","-st","-","-","-st","-st","","-R01","-R02"))
  #Adjust sampleName------------------------------------------------------------------------------------------------------------------------------------------
  data$NIM20171013_CHOEpoFcst6_b_12_3$sampleName<-"20171013-st14-R01"
  data$NIM20171013_CHOEpoFcst6_b_14_1$sampleName<-"20171013-st14-R02"
  for(i in 1:length(data)){
    data[[i]][["sampleName"]]<-paste(substr(data[[i]][["sampleName"]],10,nchar(data[[i]][["sampleName"]])),"-TP",sprintf("%02d",as.numeric(substr(data[[i]][["sampleName"]],7,8))-8),sep="")
  }
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------
  data <- avgViabDiam(data)
  #Adjust the two txt-files which were prepared manually---------------------------------------------------------------------------------
  data$NIM20171013_CHEpoFcst6_b_5_1$avgViableDiameter<-data$NIM20171013_CHEpoFcst6_b_5_1$avgViableDiameter*data$NIM20171013_CHEpoFcst6_b_5_1$viableCells/data$NIM20171014_CHEpoFcst6_b_5_1$viableCells
  data$NIM20171013_CHEpoFcst6_b_5_2$avgViableDiameter<-data$NIM20171013_CHEpoFcst6_b_5_2$avgViableDiameter*data$NIM20171013_CHEpoFcst6_b_5_2$viableCells/data$NIM20171014_CHEpoFcst6_b_5_2$viableCells
  data$NIM20171013_CHEpoFcst6_b_5_1$viablePerMlData<-data$NIM20171013_CHEpoFcst6_b_5_1$viablePerMlData*data$NIM20171013_CHEpoFcst6_b_5_1$viableCells/data$NIM20171014_CHEpoFcst6_b_5_1$viableCells
  data$NIM20171013_CHEpoFcst6_b_5_2$viablePerMlData<-data$NIM20171013_CHEpoFcst6_b_5_2$viablePerMlData*data$NIM20171013_CHEpoFcst6_b_5_2$viableCells/data$NIM20171014_CHEpoFcst6_b_5_2$viableCells
  #----------------------------------------------------------------------------------------------------------------------------------------------------------
  data <- dataTransformation(data)
  data <- assignReplicates(data, "NiM-output")
  data <- addSd(data)
  df <- getDF(data)
  df <- drawBatchGraphs(df, path = "NiM-output/batch")
  corrcombgrowthRates <- calculateGrowthCorr(df, comb = TRUE,corrSeed = TRUE)
  combgrowthRates <- calculateGrowthCorr(df, comb = TRUE)
  growthRates <- calculateGrowthCorr(df)
  outputVariable(combgrowthRates,parameter = "growthRates", path = "NiM-output/batch/growth/comb")
  outputVariable(corrcombgrowthRates,parameter = "growthRates", path =  "NiM-Output/batch/growth/comb/corr/")
  outputVariable(growthRates,parameter = "growthRates", path =  "NiM-Output/batch/growth/sep")
  drawGrowthCharts(combgrowthRates, path = "NiM-output/batch/growth/comb")
  drawGrowthCharts(corrcombgrowthRates, path = "NiM-Output/batch/growth/comb/corr/")
  drawGrowthCharts(growthRates, path = "NiM-Output/batch/growth/sep")
  df <- assignTiter(df = df, inputFile = system.file("extdata", "NiM","Titer_Input_R.txt", package = "vicellR"), outputPath = "NiM-Output/batch/prod")
  df <- calculateCumulative(df = df, outputPath = "NiM-Output/batch/prod")
  specProd <- calculateProdCorr(df, nmin = 4, corrSeed = FALSE,comb = FALSE)
  combspecProd <- calculateProdCorr(df, nmin = 4, corrSeed = TRUE,comb = TRUE)
  drawProdCharts(df=df,specProd=specProd, path = "NiM-Output/batch/prod/sep",perRow=3)
  drawProdCharts(df=df,specProd=combspecProd, path = "NiM-Output/batch/prod/comb")
  outputVariable(variable = specProd, parameter = "productivity",path = "NiM-Output/batch/prod/sep")
  outputVariable(variable = combspecProd, parameter = "productivity",path = "NiM-Output/batch/prod/comb")

  #As the median of t post seed for both growth rate calculations and spec productivity is between two measurement, it takes for some samples different timepoints--> It make sense to prepare figures where the same timepoints is used -->better comparability
  #Only prepared for one parameter (combined and seed-corrected), only necessary for VCD, VCV is at one measurement
  #Growth
  png("NiM-Output/batch/growth/comb/corr/growthRatesconcViable_TP6.png",width=600,height=400) #Creates png-file
    barcenters<-barplot(t(corrcombgrowthRates[["growthRatesconcViable"]] [,5]),ylim=c(0,1.15*max(corrcombgrowthRates[["growthRatesconcViable"]] [,5]+corrcombgrowthRates[["cigrowthRatesconcViable"]] [,5],na.rm=TRUE)), beside=TRUE,ylab="µ[h^-1]") #Draws barchart
    arrows(barcenters,t(corrcombgrowthRates[["growthRatesconcViable"]] [,5]-corrcombgrowthRates[["cigrowthRatesconcViable"]] [,5]), barcenters,t(corrcombgrowthRates[["growthRatesconcViable"]] [,5]+corrcombgrowthRates[["cigrowthRatesconcViable"]] [,5]), lwd = 1.5, angle = 90,code = 3, length = 0.05) #Adds Confidence interval
  dev.off()
  png("NiM-Output/batch/growth/comb/corr/growthRatesconcViable_TP7.png",width=600,height=400) #Creates png-file
    barcenters<-barplot(t(corrcombgrowthRates[["growthRatesconcViable"]] [,6]),ylim=c(0,1.15*max(corrcombgrowthRates[["growthRatesconcViable"]] [,6]+corrcombgrowthRates[["cigrowthRatesconcViable"]] [,6],na.rm=TRUE)), beside=TRUE,ylab="µ[h^-1]") #Draws barchart
    arrows(barcenters,t(corrcombgrowthRates[["growthRatesconcViable"]] [,6]-corrcombgrowthRates[["cigrowthRatesconcViable"]] [,6]), barcenters,t(corrcombgrowthRates[["growthRatesconcViable"]] [,6]+corrcombgrowthRates[["cigrowthRatesconcViable"]] [,6]), lwd = 1.5, angle = 90,code = 3, length = 0.05) #Adds Confidence interval
  dev.off()
  #SpecificProductivity
  png("NiM-Output/batch/prod/comb/specProdTRMaxsumccd_TP6.png",width=600,height=400) #Creates png-file
    barcenters<-barplot(t(combspecProd[["specProdsumccd"]] [,5]),ylim=c(0,1.15*max(combspecProd[["specProdsumccd"]] [,5]+combspecProd[["cispecProdsumccd"]] [,5],na.rm=TRUE)), beside=TRUE,ylab="qp[pg/(cell*d)]") #Draws barchart
    arrows(barcenters,t(combspecProd[["specProdsumccd"]] [,5]-combspecProd[["cispecProdsumccd"]] [,5]), barcenters,t(combspecProd[["specProdsumccd"]] [,5]+combspecProd[["cispecProdsumccd"]] [,5]), lwd = 1.5, angle = 90,code = 3, length = 0.05) #Adds Confidence interval
  dev.off()
  png("NiM-Output/batch/prod/comb/specProdTRMaxsumccd_TP7.png",width=600,height=400) #Creates png-file
    barcenters<-barplot(t(combspecProd[["specProdsumccd"]] [,6]),ylim=c(0,1.15*max(combspecProd[["specProdsumccd"]] [,6]+combspecProd[["cispecProdsumccd"]] [,6],na.rm=TRUE)), beside=TRUE,ylab="qp[pg/(cell*d)]") #Draws barchart
    arrows(barcenters,t(combspecProd[["specProdsumccd"]] [,6]-combspecProd[["cispecProdsumccd"]] [,6]), barcenters,t(combspecProd[["specProdsumccd"]] [,6]+combspecProd[["cispecProdsumccd"]] [,6]), lwd = 1.5, angle = 90,code = 3, length = 0.05) #Adds Confidence interval
  dev.off()
  #----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  exportCSV(df, path = "NiM-Output/batch/")
}
