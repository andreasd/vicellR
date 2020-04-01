########################################################################################################################################
#Creator: GK
#Dependencies: readTxt;
#Does: Returns the original list containing the average viable diameter(s) and the viableDiameterPerCell
#Use: The Viable cell diameter is not included in the ViCell data, and is calculated here. This is important for biomass estimations
#via cell volume and cell size-corrected growth calculations
########################################################################################################################################
avgViabDiam <- function(data){
    slot_size<-c() #Slot size for each data set
    slot_diameter<-c() #Slot diameter, for each slot
    ViabDiamDataPerCell<-c() # Transforms data  so that the diameter for each cell is saved
    tmp_diameter_per_slot<-c() #temporary slot diameter
    sum_viable_diameter<-0 # sum of slot diameters
    avg_diameter_viable_cells<-c() #sum of slot diameters divided by sum of viable cells
    #Calculates diameter-slotsize
    for(i in 1:length(data)){
      slot_size[i]<-((data[[i]] [["setMaxDiameter"]]) - (data[[i]] [["setMinDiameter"]]))/(data[[i]] [["setNumberOfBins"]]) #(Maximum - minimum diameter)/slots
    }
    #Calculates sum of diameters for each slot, sums it up, and calculates the average per cell
    for(i in 1:length(data)){
      ViabDiamDataPerCell<-0
      z<-1
      slot_diameter<-data[[i]] [["setMinDiameter"]]+(slot_size[i]/2) #Adds half of the slot size --> new slot diameter
      sum_viable_diameter[i]<-0
      for(j in 1:length(data[[i]] [["viableSizeData"]])){
        tmp_diameter_per_slot<-slot_diameter*data[[i]] [["viableSizeData"]] [j] #Number of cells per slot * slot diameter
        #Adds number of cells with the same size to a new parameter called avgViabDiamData
        if(data[[i]] [["viableSizeData"]] [j]>0){ #If at least one cell in this specific slot
          ViabDiamDataPerCell[z:(z+data[[i]] [["viableSizeData"]] [j]-1)]<-c(rep(slot_diameter,data[[i]] [["viableSizeData"]] [j])) #Adds as many diameters as cells in the slot
          z<-z+data[[i]] [["viableSizeData"]] [j] #sets counter for pointing into vector
        }
        sum_viable_diameter[i]<-sum_viable_diameter[i] + tmp_diameter_per_slot #Sums up the whole diameter for the slot + previous ones
        slot_diameter<-slot_diameter+slot_size[i] #Increases diameter for next slot
      }
      data[[i]] [["avgViableDiameter"]]<-sum_viable_diameter[i]/data[[i]] ["viableCells"] [[1]] #Calculates and adds average viable Diameter by dividing cumulative diameter by the number of viable cells
      data[[i]] [["ViableDiamDataPerCell"]]<-ViabDiamDataPerCell #Adds Viable diameter data for each cell as element
    }
    return(data) #returns list + avgViableDiameter + ViableDiamDataPerCell
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; avgViabDiam;
#Does: Returns the original list containing the average viable cell volume, Viable cell volume (VCV), and ln-transformed VCD and VCV
#Use: The Viable cell volume can be used for biomass estimations. The ln-transformed values are used for data presentations
#and growth calculations
########################################################################################################################################
dataTransformation <- function(data){
  #Calculates the average volume of a cell (avgViableCellVolume) and the volume of all cells/ml (concViableCellVolume) based on volume of a sphere
  for(i in 1:length(data)){
        data[[i]] [["avgViableCellVolume"]]<-1/6*pi*(data[[i]] [["avgViableDiameter"]]*10^(-6))^3 #Calculates the average volume of a cell (avgViableCellVolume)  based on volume of a sphere
        data[[i]] [["concViableCellVolume"]]<-data[[i]] [["concViable"]]*10^6 *data[[i]] [["avgViableCellVolume"]] #Viable cell density * avgViableCellVolume per cell = Viable Cell Volume/ml [m³/ml] (here:concViableCellVolume)
        data[[i]] [["lnconcViable"]]<-log(data[[i]] [["concViable"]]*10^6) #LN-Transforms the Viable cell density
        data[[i]] [["lnconcViableCellVolume"]]<-log(data[[i]] [["concViableCellVolume"]]*10^6) #LN-transforms the Viable Cell Volume
  }
  return(data) #returns list + avgViableCellVolume + concViableCellVolume + lnconcViable + lnconcViableCellVolume
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt
#Does: assign replicate number, replicate IDs, replicate Names and TP numbers, according to sample IDs
#Use: Is used to group data for presentation, and growth calculations
########################################################################################################################################
assignReplicates<-function(data,path){
  data<-data[order(sapply(data,'[[',"sampleName"))]
  repAssign<-matrix(NA,nrow=length(data),ncol=99) #Initializes Matrix which contains positions of samples with replicate number = column number
  repAssignNames<-matrix(NA,nrow=length(data),ncol=99) #Initializes Matrix which contains sample names with replicate number = column number, but replicate info withdrawn from sample names
  tpAssign<-matrix(NA,nrow=length(data),ncol=99) #See repAssign, but for TP
  tpAssignNames<-matrix(NA,nrow=length(data),ncol=99) #See repAssignNames, but for TP
  #Comment: what if no replicate 3, but a replicate 4 (eg one of the replicates was terminated earlier)?
  for(i in 1:99){ #Possible to assign 1 to 99 replicates
    if(length(grep(paste("-R",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName")))>0){ #If replicat with number i exist
      repAssign[1:length(grep(paste("-R",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName"))),i]<-grep(paste("-R",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName")) #The positions of all samples with a replicate number i will be written in repAssign-column i
      repAssignNames[1:length(grep(paste("-R",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName"))),i]<-mgsub(paste("-R",sprintf("%02d",i),sep=""),"",grep(paste("-R",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName"),value=TRUE)) #The names of all samples with a replicate number i will be written in repAssign-column i, but without replicate-information
    }else{ #if no replicate with number i exists
      if(length(grep(paste("-R",sprintf("%02d",1),sep=""),sapply(data,'[[',"sampleName")))==0){ #if no replicate info is present
        repAssign[1:length(data),1]<-1:length(data) #No replicate-info present-> all sample positions get assigned to column 1 of repAssign
        repAssignNames[repAssign[!is.na(repAssign[,1]),1],1]<-sapply(data,'[[',"sampleName") #No replicate-info present-> all samples names get assigned to column 1 of repAssign
      }
      break
    }
  }
  for(i in 1:99){#Possible to assign 1 to 99  timepoints (TP). Comments see replicates
    if(length(grep(paste("-TP",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName")))>0){
      tpAssign[1:length(grep(paste("-TP",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName"))),i]<-grep(paste("-TP",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName"))
      tpAssignNames[1:length(grep(paste("-TP",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName"))),i]<-mgsub(paste("-TP",sprintf("%02d",i),sep=""),"",grep(paste("-TP",sprintf("%02d",i),sep=""),sapply(data,'[[',"sampleName"),value=TRUE))
    }else{
      if(length(grep(paste("-TP",sprintf("%02d",1),sep=""),sapply(data,'[[',"sampleName")))==0){
        tpAssign[1:length(data),1]<-1:length(data)
        tpAssignNames[tpAssign[!is.na(tpAssign[,1]),1],1]<-sapply(data,'[[',"sampleName")
      }
      break
    }
  }
  #Trims Info
  repAssign<-repAssign[,apply(repAssign,2,function(x) !all(is.na(x)))] #Removes all empty columns
  if(length(dim(repAssign))==0){
    repAssign<-matrix(repAssign,ncol=1) #Generates a matrix if a vector is generated
  }else{
    repAssign<-repAssign[apply(repAssign,1,function(x) !all(is.na(x))),] #Removes all empty rows
    if(length(dim(repAssign))==0){
      repAssign<-matrix(repAssign,nrow=1) #Generates a matrix if a vector is generated
    }
  }
    repAssignNames<-repAssignNames[,apply(repAssignNames,2,function(x) !all(is.na(x)))]#Removes all empty columns
  if(length(dim(repAssignNames))==0){
    repAssignNames<-matrix(repAssignNames,ncol=1) #Generates a matrix if a vector is generated
  }else{
    repAssignNames<-repAssignNames[apply(repAssignNames,1,function(x) !all(is.na(x))),]#Removes all empty rows
    if(length(dim(repAssignNames))==0){
      repAssignNames<-matrix(repAssignNames,nrow=1) #Generates a matrix if a vector is generated
    }
  }
  tpAssign<-tpAssign[,apply(tpAssign,2,function(x) !all(is.na(x)))] #Removes all empty columns
  if(length(dim(tpAssign))==0){
    tpAssign<-matrix(tpAssign,ncol=1) #Generates a matrix if a vector is generated
  }else{
    tpAssign<-tpAssign[apply(tpAssign,1,function(x) !all(is.na(x))),]#Removes all empty rows
    if(length(dim(tpAssign))==0){
      tpAssign<-matrix(tpAssign,nrow=1) #Generates a matrix if a vector is generated
    }
  }
  tpAssignNames<-tpAssignNames[,apply(tpAssignNames,2,function(x) !all(is.na(x)))] #Removes all empty columns
  if(length(dim(tpAssignNames))==0){
    tpAssignNames<-matrix(tpAssignNames,ncol=1) #Generates a matrix if a vector is generated
  } else{
    tpAssignNames<-tpAssignNames[apply(tpAssignNames,1,function(x) !all(is.na(x))),] #Removes all empty rows
    if(length(dim(tpAssignNames))==0){
      tpAssignNames<-matrix(tpAssignNames,nrow=1) #Generates a matrix if a vector is generated
    }
  }
  #Calculates number of R and TP in total
  repNumber<-ncol(repAssign) #Calculate number of replicates based on columnes of repAssign
  tpNumber<-ncol(tpAssign) #Calculate number of timepoints based on columnes of repAssign
  repsampleNames<-unique(c(repAssignNames[!is.na(repAssignNames)])) #Extracts unique sample names without replicate info
  repsampleAssign<-matrix(NA,nrow=length(repsampleNames),ncol=repNumber) #Initializes a matrix which will contain info about sample replicate assignment
  tpsampleNames<-unique(c(tpAssignNames[!is.na(tpAssignNames)])) #Extracts unique sample names without timepoint info
  tpsampleAssign<-matrix(NA,nrow=length(tpsampleNames),ncol=tpNumber) #Initializes a matrix which will contain info about sample timepoint assignment
  #Extracts positions of samples with the same IDs in the list
  for(i in 1:length(repsampleNames)){
    for(j in 1:repNumber){
      if(length(grep(repsampleNames[i],repAssignNames[,j]))>0){
        repsampleAssign[i,j]<-repAssign[grep(repsampleNames[i],repAssignNames[,j]),j] #Extracts sample positions in the list, and puts all samples with the same ID (without rep-info) into a row
      }else{
        repsampleAssign[i,j]<-NA
      }
    }
  }
  for(i in 1:length(tpsampleNames)){
    for(j in 1:tpNumber){
      if(length(grep(tpsampleNames[i],tpAssignNames[,j]))>0){
        tpsampleAssign[i,j]<-tpAssign[grep(tpsampleNames[i],tpAssignNames[,j]),j] #Extracts sample positions in the list, and puts all samples with the same ID (without tp-info) into a row
      }else{
        tpsampleAssign[i,j]<-NA
      }
    }
  }
  #Add replicate-info to data-list
  for(i in 1:nrow(repsampleAssign)){
    for(j in 1:ncol(repsampleAssign)){
      if(!is.na(repsampleAssign[i,j])){
        data[[repsampleAssign[i,j]]] [["replicateNumber"]]<-j #Adds replicate number to list
      }
    }
  }
  for(i in 1:nrow(tpsampleAssign)){
    for(j in 1:ncol(tpsampleAssign)){
      if(!is.na(tpsampleAssign[i,j])){
        data[[tpsampleAssign[i,j]]] [["TPNumber"]]<-j #Adds tp number to list
      }
    }
  }
  assignNames<-gsub("-TP[0-9][0-9]","",repAssignNames) #Generates sample names without rep and tp info
  sampleNames<-unique(c(assignNames)) #unique sample names without rep and tp info
  for(i in 1:length(sampleNames)){
    for(j in 1:ncol(assignNames)){
      for(k in which(assignNames[,j] %in% sampleNames[i])){
        if(!is.na(repAssign[k,j])){
          data[[repAssign[k,j]]] [["replicateID"]]<-i #Adds replicateID number to all samples with the same sample name without rep and tp info
          data[[repAssign[k,j]]] [["replicateName"]]<-sampleNames[i] #Adds sample name without rep and tp info to all samples with the same sample name without rep and tp info
        }
      }
    }
  }
  #Comparison between Replicate 1...n
  if(repNumber>1){ #If more than one replicate
    #Initialize matrices which contain information to be compared
    compViability<-matrix(NA,ncol=repNumber,nrow=length(repsampleAssign[,1]))
    compConcTotal<-matrix(NA,ncol=repNumber,nrow=length(repsampleAssign[,1]))
    compConcViable<-matrix(NA,ncol=repNumber,nrow=length(repsampleAssign[,1]))
    compAvgDiameter<-matrix(NA,ncol=repNumber,nrow=length(repsampleAssign[,1]))
    compAvgViableDiameter<-matrix(NA,ncol=repNumber,nrow=length(repsampleAssign[,1]))
    compConcViableCellVolume<-matrix(NA,ncol=repNumber,nrow=length(repsampleAssign[,1]))
    for(j in 1:length(repsampleAssign[,1])){ #Add data to matrices which should be compared
      for(i in 1:repNumber){
        if(!is.na(repsampleAssign[j,i])){
          compViability[j,i]<-data[[repsampleAssign[j,i]]] [["viability"]]
          compConcTotal[j,i]<-data[[repsampleAssign[j,i]]] [["concTotal"]]
          compConcViable[j,i]<-data[[repsampleAssign[j,i]]] [["concViable"]]
          compAvgDiameter[j,i]<-data[[repsampleAssign[j,i]]] [["avgDiameter"]]
          compAvgViableDiameter[j,i]<-data[[repsampleAssign[j,i]]] [["avgViableDiameter"]]
          compConcViableCellVolume[j,i]<-data[[repsampleAssign[j,i]]] [["concViableCellVolume"]]
        }
      }
    }
    pdf(file=paste(path,"/replicate_comparison.pdf",sep="")) #Generates a pdf-file with all dot plots between R1 and R2 to RN-1 and RN
    for(i in 2:repNumber){
      plot(compViability[,i-1]~compViability[,i],xlab=paste("R",i-1,sep="-"),ylab=paste("R",i,sep="-"),main="viability [%]")
      legend("topleft",paste("R=",round(cor(compViability[,i-1],compViability[,i],use="pairwise.complete.obs"),digits=4)),bty="n")
      plot(compConcTotal[,i-1]~compConcTotal[,i],xlab=paste("R",i-1,sep="-"),ylab=paste("R",i,sep="-"),main="TCD [e6 cells/ml]")
      legend("topleft",paste("R=",round(cor(compConcTotal[,i-1],compConcTotal[,i],use="pairwise.complete.obs"),digits=4)),bty="n")
      plot(compConcViable[,i-1]~compConcViable[,i],xlab=paste("R",i-1,sep="-"),ylab=paste("R",i,sep="-"),main="VCD [e6 cells/ml]")
      legend("topleft",paste("R=",round(cor(compConcViable[,i-1],compConcViable[,i],use="pairwise.complete.obs"),digits=4)),bty="n")
      plot(compAvgDiameter[,i-1]~compAvgDiameter[,i],xlab=paste("R",i-1,sep="-"),ylab=paste("R",i,sep="-"),main="Average diameter [µm]")
      legend("topleft",paste("R=",round(cor(compAvgDiameter[,i-1],compAvgDiameter[,i],use="pairwise.complete.obs"),digits=4)),bty="n")
      plot(compAvgViableDiameter[,i-1]~compAvgViableDiameter[,i],xlab=paste("R",i-1,sep="-"),ylab=paste("R",i,sep="-"),main="Average viable diameter [µm]")
      legend("topleft",paste("R=",round(cor(compAvgViableDiameter[,i-1],compAvgViableDiameter[,i],use="pairwise.complete.obs"),digits=4)),bty="n")
      plot(compConcViableCellVolume[,i-1]~compConcViableCellVolume[,i],xlab=paste("R",i-1,sep="-"),ylab=paste("R",i,sep="-"),main="VCV [m³/ml]")
      legend("topleft",paste("R=",round(cor(compConcViableCellVolume[,i-1],compConcViableCellVolume[,i],use="pairwise.complete.obs"),digits=4)),bty="n")
    }
    dev.off()
  }
  return(data) #returns data list + replicateNumber + TPnumber + replicateID + replicateName
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; avgViabDiam; dataTransformation
#Does: Adds confidence intervals for VCD, ln(VCD), viableDiameter, average viable cell volume, VCV and ln(VCV)
#Use: Used for error representations at plots
########################################################################################################################################
addSd<-function(data){
  for(i in 1:length(data)){ # For all samples in the list
    concFactorViable<-sum(data[[i]] [["viablePerMlData"]])/(data[[i]] [["viableCells"]]) #Determines the factor between number of viable cells and cell density
    avgViableCellsImage<-data[[i]] [["viableCells"]]/data[[i]] [["images"]] #Determines the average number of viable cells per image (needed for poisson-distribution)
    data[[i]] [["ciconcViable"]]<-paste0(c(((avgViableCellsImage+1.96*sqrt(avgViableCellsImage/data[[i]] [["images"]]))*concFactorViable),((avgViableCellsImage-1.96*sqrt(avgViableCellsImage/data[[i]] [["images"]])))*concFactorViable),collapse="||") #Adds Confidence interval (95%) based on Poisson-distribution and normal-approximation
    data[[i]] [["cilnconcViable"]]<-paste0(c(log((avgViableCellsImage+1.96*sqrt(avgViableCellsImage/data[[i]] [["images"]]))*concFactorViable*10^6),log((avgViableCellsImage-1.96*sqrt(avgViableCellsImage/data[[i]] [["images"]]))*concFactorViable*10^6)),collapse="||") #LN-transforms Confidence interval

    data[[i]] [["civiability"]] <- paste0(c(((data[[i]] [["viability"]]/100 +1.96*sqrt(1/data[[i]] [["totalCells"]]*(data[[i]] [["viability"]]/100)*(1-data[[i]] [["viability"]]/100)))*100),((data[[i]] [["viability"]]/100-1.96*sqrt(1/data[[i]] [["totalCells"]]*data[[i]] [["viability"]]/100*(1-data[[i]] [["viability"]]/100)))*100)),collapse="||") #Add confidence interval for Bernoulli-distribution (percentage data), using normal approximation

    data[[i]] [["ciavgViableDiameter"]]<-paste0(c((data[[i]] [["avgViableDiameter"]]+1.96*(sd(data[[i]] [["ViableDiamDataPerCell"]])/sqrt(data[[i]] [["viableCells"]]))),(data[[i]] [["avgViableDiameter"]]-1.96*(sd(data[[i]] [["ViableDiamDataPerCell"]])/sqrt(data[[i]] [["viableCells"]])))),collapse="||") #Add confidence interval for Normal distribution, if no better solution
    data[[i]] [["ciavgViableCellVolume"]]<-paste0(c((1/6*pi*(as.numeric(strsplit(data[[i]] [["ciavgViableDiameter"]],split="||",fixed=TRUE) [[1]] [1])*10^(-6))^3),1/6*pi*(as.numeric(strsplit(data[[i]] [["ciavgViableDiameter"]],split="||",fixed=TRUE) [[1]] [2])*10^(-6))^3),collapse="||") #Transform by same treatment as avgViableCellVolume

    ciconcViableCellVolume<-data[[i]] [["concViableCellVolume"]]*((as.numeric(strsplit(data[[i]] [["ciavgViableCellVolume"]],split="||",fixed=TRUE) [[1]] [1])-as.numeric(strsplit(data[[i]] [["ciavgViableCellVolume"]],split="||",fixed=TRUE) [[1]] [2]))/data[[i]] [["avgViableCellVolume"]]+(as.numeric(strsplit(data[[i]] [["ciconcViable"]],split="||",fixed=TRUE) [[1]] [1])-as.numeric(strsplit(data[[i]] [["ciconcViable"]],split="||",fixed=TRUE) [[1]] [2]))/data[[i]] [["concViable"]]) #Calculate new confidence intervals by transforming ci to a relative range, and do error propagation, using gy=g1+g2
    data[[i]] [["ciconcViableCellVolume"]]<-paste0(c((data[[i]] [["concViableCellVolume"]]+ciconcViableCellVolume/2),(data[[i]] [["concViableCellVolume"]]-ciconcViableCellVolume/2)),collapse="||") #Transfer calculated confidence interval into list
    data[[i]] [["cilnconcViableCellVolume"]]<-paste0(c((log(as.numeric(strsplit(data[[i]] [["ciconcViableCellVolume"]],split="||",fixed=TRUE) [[1]] [1])*10^6)),(log(as.numeric(strsplit(data[[i]] [["ciconcViableCellVolume"]],split="||",fixed=TRUE) [[1]] [2])*10^6))),collapse="||") #LN-transforms confidence interval
  }
  return(data) #returns list + added standard deviations
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; avgViabDiam; dataTransformation; assignReplicates; addSd; getDF;
#Does: Creates vs-t figures of all samples in the batch for VCD, ln(VCD), VCV, ln(VCV), viability, average viable cell Volume
#Use: For culture representations and comparisons
########################################################################################################################################
drawBatchGraphs<-function(df,seedTime=NA,path,col=NA,perRow=NA){
  if(is.na(seedTime)){ # Calculates time post seed in hours, and also determines maximum x-axis-values for all graphs
    xmaxTime<-(as.numeric(max(df[,"runDate"]))-as.numeric(min(df[,"runDate"])))/3600
    tPostSeed<-(as.numeric(df[,"runDate"])-as.numeric(min(df[,"runDate"])))/3600
  } else { #For offset between seeding and first measurement
    xmaxTime<-(as.numeric(max(df[,"runDate"]))-as.numeric(as.POSIXct(seedTime,origin='1970-01-01 00:00:00 UTC')))/3600
    tPostSeed<-(as.numeric(df[,"runDate"])-as.numeric(as.POSIXct(seedTime,origin='1970-01-01 00:00:00 UTC')))/3600
  }
  df<-cbind(df,tPostSeed) #Adds time post seed to dataframe
  if(is.na(col)){ #Defines colours for plots
    col<-1:length(unique(df[,"replicateNumber"]))
  }
  drawParvsTime(df=df,parameter="concViable",col=col,path=path,xmaxTime=xmaxTime,ylabel="VCD [1E6 cells/ml]",legend="topleft",perRow=perRow) #Executes function that draws VCDvstime
  drawParvsTime(df=df,parameter="lnconcViable",col=col,path=path,xmaxTime=xmaxTime,ylabel="ln(VCD [cells/ml])",legend="topleft",perRow=perRow)  #Executes function that draws lnVCDvstime
  drawParvsTime(df=df,parameter="avgViableCellVolume",col=col,path=path,xmaxTime=xmaxTime,ylabel="Average Viable Cell Volume [m³/cell]",legend="topright",perRow=perRow)  #Executes function that draws Average VCVvstime
  drawParvsTime(df=df,parameter="concViableCellVolume",col=col,path=path,xmaxTime=xmaxTime,ylabel="VCV [m³/ml]",legend="topleft",perRow=perRow)  #Executes function that draws VCVvstime
  drawParvsTime(df=df,parameter="lnconcViableCellVolume",col=col,path=path,xmaxTime=xmaxTime,ylabel="ln (VCV [m³/ml])",legend="topleft",perRow=perRow)  #Executes function that draws lnVCVvstime
  drawParvsTime(df=df,parameter="viability",col=col,path=path,xmaxTime=xmaxTime,ylabel="Viability [%]",legend="bottomleft",perRow=perRow) #Executes function that draws Viabilityvstime
  return(df) #Return df + tpostseed
}
#Used function from drawBatchGraphs
drawParvsTime<-function(df,parameter,col,path,xmaxTime,ylabel,legend,suppCI=FALSE,perRow=NA){ #Draws vs time-plot
  if(suppCI==FALSE){
    ymax<-max(as.numeric(unlist(strsplit(df[,paste("ci",parameter,sep="")],split="||",fixed=TRUE)))) #Determines max y-axis-value
    ymin<-min(as.numeric(unlist(strsplit(df[,paste("ci",parameter,sep="")],split="||",fixed=TRUE)))) #Determines min y-axis-value
  }else{
    ymax<-max(df[,paste(parameter)],na.rm=TRUE)
    ymin<-min(df[,paste(parameter)],na.rm=TRUE)
  }
  if(!is.na(perRow)){
    png(paste(path,"/",parameter,"vstime.png",sep=""),width=400*perRow,height=(400*length(unique(df[,"replicateID"]))/perRow))
    par(mfrow=c(ceiling(length(unique(df[,"replicateID"]))/perRow),perRow)) #Determines number of plots
  }else if(length(unique(df[,"replicateID"]))>1){
    png(paste(path,"/",parameter,"vstime.png",sep=""),width=800,height=(400*length(unique(df[,"replicateID"]))/2))
    par(mfrow=c(ceiling(length(unique(df[,"replicateID"]))/2),2)) #Determines number of plots
  } else{
    png(paste(path,"/",parameter,"vstime.png",sep=""),width=800,height=400) #Creates png-file
  }
  for(i in unique(df[,"replicateID"])){ #Creates plot for each replicate ID
    par(new=FALSE)
    temp<-df[df[,"replicateID"]==i,] #Subsample
    gtitle<-temp[1,c("replicateName")] # Title selection for each sample
    for(j in unique(temp[,"replicateNumber"])){
      temp2<-temp[temp[,"replicateNumber"]==j,] #Subsample
      for(k in 1:length(parameter)){
        plot(temp2[,c("tPostSeed")],temp2[,paste(parameter)[k]],xlab="t[h]",ylab= paste(ylabel),main=gtitle,ylim=c(ymin,ymax),xlim=c(0,xmaxTime),col=col[j],pch=20,cex.axis=1.5,cex.lab=1.5,cex.main=2.5)# Plots data
        lines(temp2[,c("tPostSeed")],temp2[,paste(parameter)[k]],type="l",col=col[j]) #Adds lines between points
        if(suppCI==FALSE){
          arrows(temp2[,c("tPostSeed")],as.numeric(sapply(strsplit(temp2[,paste("ci",parameter,sep="")],split="||",fixed=TRUE),'[[',2)),temp2[,c("tPostSeed")],as.numeric(sapply(strsplit(temp2[,paste("ci",parameter,sep="")],split="||",fixed=TRUE),'[[',1)), lwd = 1.5, angle = 90,code = 3, length = 0.05,col=col[j]) #Adds confidence interval
        }
        par(new=TRUE)#So that replicate will be added into same plot
      }
    }
    legend(paste(legend),legend=paste(sort(unique(temp[,"replicateNumber"]))),col=col,pch=20,cex=1.5) #Adds replicate legend
  }
  dev.off() #Closes file

}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; avgViabDiam; dataTransformation; getDF;
#Does: Calculates growth rates based on maximum correlation coefficients and median of timepoints of maximum correlations between time and log-transformed values (VCD, VCV )
#Use: For growth representations and culture comparisons
########################################################################################################################################
calculateGrowthCorr<-function(df,nmin=4,corrSeed=FALSE,comb=FALSE,lagDelay=0){
  if(comb==FALSE){ #If replicates should not be combined
    sepreplicateID<-nrow(df[,1]) #New replicateID is initialized
    for(i in 1:nrow(unique(df[,c("replicateID","replicateNumber")]))){
      sepreplicateID[df[,c("replicateID")] %in% unique(df[,c("replicateID","replicateNumber")]) [i,1] & df[,c("replicateNumber")] %in% unique(df[,c("replicateID","replicateNumber")]) [i,2]]<-i #New replicateID is created, based on differences in replicateID and replicateNumber
    }
    for(i in 1:nrow(df)){
      df[i,c("replicateName")]<-paste0(df[i,c("replicateName","replicateNumber")],collapse="-") #New replicateName is created (replicateID+replicateNumber)
    }
    df[,c("replicateID")]<-sepreplicateID #New replicateID replaces old one
  } else{
    #Check whether samples overlap in their confidence interval
    if(corrSeed==TRUE & max(df[,"replicateNumber"])>1){ #if replicates should be combined, and seed-corrected
      corrFactorlnconcViable<-c() #Initializes correctionfactors for ln(VCD)
      corrFactorlnconcViableCellVolume<-c()  #Initializes correctionfactors for ln(VCV)
      for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID separately
        for(j in 2:max(df[,"replicateNumber"])){ #For each replicateNumber except one
          corrFactorlnconcViable<-df[df[,"replicateID"]==i & df[,"TPNumber"]==1 & df[,"replicateNumber"]==j,"lnconcViable"] - df[df[,"replicateID"]==i & df[,"TPNumber"]==1 & df[,"replicateNumber"]==1,"lnconcViable"] #Difference between ln(VCD) of replicate j and replicate 1 at TP01-->Correction factor
          corrFactorlnconcViableCellVolume<-df[df[,"replicateID"]==i & df[,"TPNumber"]==1 & df[,"replicateNumber"]==j,"lnconcViableCellVolume"] - df[df[,"replicateID"]==i & df[,"TPNumber"]==1 & df[,"replicateNumber"]==1,"lnconcViableCellVolume"] #Difference between ln(VCV) of replicate j and replicate 1 at TP01-->Correction factor
          df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"lnconcViable"]<-df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"lnconcViable"]-corrFactorlnconcViable #Apply correction factor by substraction
          df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"lnconcViableCellVolume"]<-df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"lnconcViableCellVolume"]-corrFactorlnconcViableCellVolume  #Apply correction factor by substraction
        }
      }
    }
    if(max(df[,"replicateNumber"])>1){
     l<-0 #indicates whether a warning has been produced
    } else {
      l<-1
    }
    for(i in 1:length(unique(df[,"replicateID"]))){ #for all replicateIDs
      if(max(df[,"replicateNumber"])>1)
      if(l==1){
        break
      }
      temp<-df[df[,"replicateID"]==i,] #Subsample
      for(j in 1:max(temp[,"TPNumber"])){ #For all timepoints
        if(l==1){
          break
        }
        temp2<-temp[temp[,"TPNumber"]==j,] #Subsample
        if(length(unique(temp2[,"replicateNumber"]))>1){
          for(k in unique(df[,"replicateNumber"]) [1:(length(unique(df[,"replicateNumber"]))-1)]){ #For all replicates
            if(sum(sapply(strsplit(temp2[temp2[,"replicateNumber"]>k,"cilnconcViable"],split="||",fixed=TRUE),'[[',1) < strsplit(temp2[temp2[,"replicateNumber"]==k,"cilnconcViable"],split="||",fixed=TRUE) [[1]] [2] | sapply(strsplit(temp2[temp2[,"replicateNumber"]>k,"cilnconcViable"],split="||",fixed=TRUE),'[[',2) > strsplit(temp2[temp2[,"replicateNumber"]==k,"cilnconcViable"],split="||",fixed=TRUE) [[1]] [1])>0){ # If lower ci limit of at least one replicate >k is larger than higher ci limit of replicate k, or if higher ci limit of at least one replicate >k is lower than lower ci limit of replicate k (No overlap between Confidence intervalls of lnconcViable)--> Warning produced
              warning("Confidence intervals don't overlap between replicates, do seed-correction or calculate growth rates for replicates separately")
              l<-1 #l is set to one, so that all for-loops are determined --> Only one warning produced
              break
            }else if(sum(sapply(strsplit(temp2[temp2[,"replicateNumber"]>k,"cilnconcViableCellVolume"],split="||",fixed=TRUE),'[[',1) < strsplit(temp2[temp2[,"replicateNumber"]==k,"cilnconcViableCellVolume"],split="||",fixed=TRUE) [[1]] [2] | sapply(strsplit(temp2[temp2[,"replicateNumber"]>k,"cilnconcViableCellVolume"],split="||",fixed=TRUE),'[[',2) > strsplit(temp2[temp2[,"replicateNumber"]==k,"cilnconcViableCellVolume"],split="||",fixed=TRUE) [[1]] [1])>0){ # If lower ci limit of at least one replicate >k is larger than higher ci limit of replicate k, or if higher ci limit of at least one replicate >k is lower than lower ci limit of replicate k (No overlap between Confidence intervalls of lnconcViableCellVolume)--> Warning produced
              warning("Confidence intervals don't overlap between replicates, do seed-correction or calculate growth rates for replicates separately")
              l<-1 #l is set to one, so that all for-loops are determined --> Only one warning produced
              break
            }
          }
        }else{
          warning("Only one replicate at a timepoint")
        }
      }
    }
  }
  growthRates<-list()#Initialize variables - will contain growth data
  rconcViable<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]))))  #Initialize variable - Correlation coefficients between lnconcViable and runDate (-1, because it starts with element 1-2)
  rconcViableCellVolume<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"])))) #Initialize variable - Correlation coefficients between lnconcViableCellVolume and runDate (-1, because it starts with element 1-2)
  rMaxconcViable<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Maximum Correlation coefficients between lnconcViable and runDate
  rMaxconcViableCellVolume<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Maximum Correlation coefficients between lnconcViableCellVolume and runDate
  growthRatesconcViable<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"])))) #Initialize variable - Growth rates based on slope between lnconcViable and runDate (-1, because it starts with element 1-2)
  growthRatesconcViableCellVolume<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"])))) #Initialize variable - Growth rates based on slope between lnconcViableCellVolume and runDate (-1, because it starts with element 1-2)
  growthRatesRMaxconcViable<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Growth rates based on slope between lnconcViable and runDate at sampling point with maximum Correlation coefficients
  growthRatesRMaxconcViableCellVolume<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Growth rates based on slope between lnconcViableCellVolume and runDate at sampling point with maximum Correlation coefficients
  growthRatesTRMaxconcViable<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Growth rates based on slope between lnconcViable and runDate at sampling point with median maximum Correlation coefficients
  growthRatesTRMaxconcViableCellVolume<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Growth rates based on slope between lnconcViableCellVolume and runDate at sampling point with median maximum Correlation coefficients
  cigrowthRatesconcViable<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"])))) #Initialize variable - Confidence interval of Growth rates based on slope between lnconcViable and runDate (-1, because it starts with element 1-2)
  cigrowthRatesconcViableCellVolume<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]))))   #Initialize variable - Confidence interval of Growth rates based on slope between lnconcViableCellVolume and runDate (-1, because it starts with element 1-2)
  cigrowthRatesRMaxconcViable<-c(rep(NA,length(unique(df[,"replicateID"]))))   #Initialize variable - Confidence interval of Growth rates based on slope between lnconcViable and runDate at sampling point with maximum Correlation coefficients
  cigrowthRatesRMaxconcViableCellVolume<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Confidence interval of Growth rates based on slope between lnconcViableCellVolume and runDate at sampling point with maximum Correlation coefficients
  cigrowthRatesTRMaxconcViable<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Confidence interval of Growth rates based on slope between lnconcViable and runDate at sampling point with median maximum Correlation coefficients
  cigrowthRatesTRMaxconcViableCellVolume<-c(rep(NA,length(unique(df[,"replicateID"]))))   #Initialize variable - Confidence interval of Growth rates based on slope between lnconcViableCellVolume and runDate at sampling point with median maximum Correlation coefficients
  tRMaxconcViable<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Timepoints of best correlation between lnconcViable and runDate
  tRMaxconcViableCellVolume<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Timepoints of best correlation between lnconcViableCellVolume and runDate
  rownames(rconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(rconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(growthRatesconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(growthRatesconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(cigrowthRatesconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(cigrowthRatesconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(rMaxconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(rMaxconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(growthRatesRMaxconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(growthRatesRMaxconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(growthRatesTRMaxconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(growthRatesTRMaxconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cigrowthRatesRMaxconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cigrowthRatesRMaxconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cigrowthRatesTRMaxconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cigrowthRatesTRMaxconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(tRMaxconcViable)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(tRMaxconcViableCellVolume)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  colnames(rconcViable)<-paste("TP",(1+lagDelay),"-",unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]),sep="") #Add names to matrices/vectors
  colnames(rconcViableCellVolume)<-paste("TP",(1+lagDelay),"-",unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]),sep="") #Add names to matrices/vectors
  colnames(growthRatesconcViable)<-paste("TP",(1+lagDelay),"-",unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]),sep="") #Add names to matrices/vectors
  colnames(growthRatesconcViableCellVolume)<-paste("TP",(1+lagDelay),"-",unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]),sep="") #Add names to matrices/vectors
  colnames(cigrowthRatesconcViable)<-paste("TP",(1+lagDelay),"-",unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]),sep="") #Add names to matrices/vectors
  colnames(cigrowthRatesconcViableCellVolume)<-paste("TP",(1+lagDelay),"-",unique(df[df[,"TPNumber"]>(1+lagDelay),"TPNumber"]),sep="") #Add names to matrices/vectors
  for(n in unique(df[,"replicateID"])){ #For each replicateID
    temp<-df[df[,"replicateID"]==n,] #Subsample
    for(i in unique(temp[temp[,"TPNumber"]>(1+lagDelay),"TPNumber"])){ #Tp1-Tpi
      timepoints<-as.numeric(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"runDate"])/3600 #Calculate time in hours
      rconcViable[n,(i-(1+lagDelay))]<-summary(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViable"]~timepoints))$r.squared #determines correlation coefficients of linear correlations between lnVCD and time of TP1-TPi
      rconcViableCellVolume[n,(i-(1+lagDelay))]<-summary(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViableCellVolume"]~timepoints))$r.squared #determines correlation coefficients of linear correlations between lnVCD and time of TP1-TPi
      #Growth rates
      growthRatesconcViable[n,(i-(1+lagDelay))]<-coef(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViable"]~timepoints))[2] #determines slope (=growth rate/h) of linear correlations between lnVCD and time of TP1-TPi
      growthRatesconcViableCellVolume[n,(i-(1+lagDelay))]<-coef(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViableCellVolume"]~timepoints))[2] #determines slope (=growth rate/h) of linear correlations between lnVCV and time of TP1-TPi
      w<-getOption('warn') #Gets warning-option
      options(warn=-1) #Suppresses warnings, because next cmd can generate a warning, if only two values are available
      cigrowthRatesconcViable[n,(i-(1+lagDelay))]<-confint(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViable"]~timepoints))[2,2]-coef(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViable"]~timepoints))[2] #calculates confidence interval of linear correlation between lnVCD and time of TP1-TPi
      if(is.na(cigrowthRatesconcViable[n,(i-(1+lagDelay))])){ #If NA is generated for a confidence interval, it is replaced by a 0
        cigrowthRatesconcViable[n,(i-(1+lagDelay))]<-0
      }
      cigrowthRatesconcViableCellVolume[n,(i-(1+lagDelay))]<-confint(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViableCellVolume"]~timepoints))[2,2]-coef(lm(temp[temp[,"TPNumber"] %in% c((1+lagDelay):i),"lnconcViableCellVolume"]~timepoints))[2] #calculates confidence interval of linear correlation between lnVCV and time of TP1-TPi
      if(is.na(cigrowthRatesconcViableCellVolume[n,(i-(1+lagDelay))])){ #If NA is generated for a confidence interval, it is replaced by a 0
        cigrowthRatesconcViableCellVolume[n,(i-(1+lagDelay))]<-0
      }
      options(warn=w) #Warnings are turned back on
    }
  }

  for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    if(sum(!is.na(rconcViable[i,c((nmin-1):length(rconcViable[1,]))]))>0){ #If rconcViable not only consists of NA for at least nmin timepoints
      rMaxconcViable[i]<-which(rconcViable[i,]==max(rconcViable[i,c((nmin-1):length(rconcViable[1,]))],na.rm=TRUE)) #Maximum Correlation coefficents between lnVCD and TP are extracted (at least nmin timepoints!)
    }else{ #If no correlation coefficients are stored, an NA is added
      rMaxconcViable[i]<-NA
    }
    if(sum(!is.na(rconcViableCellVolume[i,c((nmin-1):length(rconcViableCellVolume[1,]))]))>0){ #If rconcViableCellVolume not only consists of NA for at least nmin timepoints
      rMaxconcViableCellVolume[i]<-which(rconcViableCellVolume[i,]==max(rconcViableCellVolume[i,c((nmin-1):length(rconcViableCellVolume[1,]))],na.rm=TRUE)) #Maximum Correlation coefficents between lnVCV and TP are extracted (at least nmin timepoints!)
    }else{#If no correlation coefficients are stored, an NA is added
      rMaxconcViableCellVolume[i]<-NA
    }
  }
  for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    tRMaxconcViable[i]<-mean(as.numeric(df[df[,"replicateID"]==unique(df[,"replicateID"]) [i] & df[,"TPNumber"]==(rMaxconcViable[i]+1+lagDelay),c("runDate")])) #Timepoints of best correlations between lnVCD and time are extracted
    tRMaxconcViableCellVolume[i]<-mean(as.numeric(df[df[,"replicateID"]==unique(df[,"replicateID"]) [i] & df[,"TPNumber"]==(rMaxconcViableCellVolume[i]+1+lagDelay),c("runDate")])) #Timepoints of best correlations between lnVCV and time are extracted
  }
  #Determine growth rates based on maximum correlation coefficients
  for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    growthRatesRMaxconcViable[i]<-growthRatesconcViable[i,rMaxconcViable[i]] #lnVCD-based growth rates based on maximum correlation coefficient are extracted
    growthRatesRMaxconcViableCellVolume[i]<-growthRatesconcViableCellVolume[i,rMaxconcViableCellVolume[i]] #lnVCV-based growth rates based on maximum correlation coefficient are extracted
    cigrowthRatesRMaxconcViable[i]<-cigrowthRatesconcViable[i,rMaxconcViable[i]] #lnVCD-based confidence interval based on maximum correlation coefficient are extracted
    cigrowthRatesRMaxconcViableCellVolume[i]<-cigrowthRatesconcViableCellVolume[i,rMaxconcViableCellVolume[i]] #lnVCV-based confidence interval based on maximum correlation coefficient are extracted
  }
  temp2<-c() #Clears subsample
  for(i in unique(df[,"replicateID"])){ #For each replicateID
    for(j in (nmin+lagDelay):length(unique(df[,"TPNumber"]))){ #For each TPNumber
      temp2[j]<-mean(as.numeric(df[df[,"replicateID"]==i & df[,"TPNumber"]==j,"runDate"])) #Extracts average (if more than one sample per timepoint) runDate for each TP
    }
    if(length(temp2[!is.na(temp2)])>0){
      growthRatesTRMaxconcViable[i]<-growthRatesconcViable[i,which(abs(median(tRMaxconcViable,na.rm=TRUE)-temp2)==min(abs(median(tRMaxconcViable,na.rm=TRUE)-temp2),na.rm=TRUE))-1-lagDelay] #Extracts growth rate based on lnVCD vs time of runDate which is closest to the median runDate of max. Correlation Coefficient (-1 because it starts with TP0-1)
      growthRatesTRMaxconcViableCellVolume[i]<-growthRatesconcViableCellVolume[i,which(abs(median(tRMaxconcViableCellVolume,na.rm=TRUE)-temp2)==min(abs(median(tRMaxconcViableCellVolume,na.rm=TRUE)-temp2),na.rm=TRUE))-1-lagDelay] #Extracts growth rate based on lnVCV vs time of runDate which is closest to the median runDate of max. Correlation Coefficient (-1 because it starts with TP0-1)
      cigrowthRatesTRMaxconcViable[i]<-cigrowthRatesconcViable[i,which(abs(median(tRMaxconcViable,na.rm=TRUE)-temp2)==min(abs(median(tRMaxconcViable,na.rm=TRUE)-temp2),na.rm=TRUE))-1-lagDelay] #Extracts confidence interval of growth rate based on lnVCD vs time of runDate which is closest to the median runDate of max. Correlation Coefficient (-1 because it starts with TP0-1)
      cigrowthRatesTRMaxconcViableCellVolume[i]<-cigrowthRatesconcViableCellVolume[i,which(abs(median(tRMaxconcViableCellVolume,na.rm=TRUE)-temp2)==min(abs(median(tRMaxconcViableCellVolume,na.rm=TRUE)-temp2),na.rm=TRUE))-1-lagDelay]#Extracts confidence interval of growth rate based on lnVCD vs time of runDate which is closest to the median runDate of max. Correlation Coefficient (-1 because it starts with TP0-1)
    }else{
      growthRatesTRMaxconcViable[i]<-NA
      growthRatesTRMaxconcViableCellVolume[i]<-NA
      cigrowthRatesTRMaxconcViable[i]<-NA
      cigrowthRatesTRMaxconcViableCellVolume[i]<-NA
    }
  }
  growthRates[["growthRatesconcViable"]]<-growthRatesconcViable #Adds values to growthRates
  growthRates[["growthRatesconcViableCellVolume"]]<-growthRatesconcViableCellVolume #Adds values to growthRates
  growthRates[["cigrowthRatesconcViable"]]<-cigrowthRatesconcViable #Adds values to growthRates
  growthRates[["cigrowthRatesconcViableCellVolume"]]<-cigrowthRatesconcViableCellVolume #Adds values to growthRates
  growthRates[["rconcViable"]]<-rconcViable #Adds values to growthRates
  growthRates[["rconcViableCellVolume"]]<-rconcViableCellVolume #Adds values to growthRates
  growthRates[["rMaxconcViable"]]<-rMaxconcViable #Adds values to growthRates
  growthRates[["rMaxconcViableCellVolume"]]<-rMaxconcViableCellVolume #Adds values to growthRates
  growthRates[["median(rMaxconcViable)"]]<-median(rMaxconcViable,na.rm=TRUE) #Adds values to growthRates
  growthRates[["median(rMaxconcViableCellVolume)"]]<-median(rMaxconcViableCellVolume,na.rm=TRUE) #Adds values to growthRates
  growthRates[["growthRatesRMaxconcViable"]]<-growthRatesRMaxconcViable #Adds values to growthRates
  growthRates[["growthRatesRMaxconcViableCellVolume"]]<-growthRatesRMaxconcViableCellVolume #Adds values to growthRates
  growthRates[["cigrowthRatesRMaxconcViable"]]<-cigrowthRatesRMaxconcViable #Adds values to growthRates
  growthRates[["cigrowthRatesRMaxconcViableCellVolume"]]<-cigrowthRatesRMaxconcViableCellVolume #Adds values to growthRates
  growthRates[["tRMaxconcViable"]]<-tRMaxconcViable#Adds values to growthRates
  growthRates[["tRMaxconcViableCellVolume"]]<-tRMaxconcViableCellVolume  #Adds values to growthRates
  growthRates[["median(tRMaxconcViable)"]]<-median(tRMaxconcViable,na.rm=TRUE) #Adds values to growthRates
  growthRates[["Date of median(tRMaxconcViable)"]]<-as.POSIXct(median(tRMaxconcViable,na.rm=TRUE),origin='1970-01-01 00:00:00 UTC') #Adds values to growthRates
  growthRates[["median(tRMaxconcViableCellVolume)"]]<-median(tRMaxconcViableCellVolume,na.rm=TRUE) #Adds values to growthRates
  growthRates[["Date of median(tRMaxconcViableCellVolume)"]]<-as.POSIXct(median(tRMaxconcViableCellVolume,na.rm=TRUE),origin='1970-01-01 00:00:00 UTC')
  growthRates[["growthRatesTRMaxconcViable"]]<-growthRatesTRMaxconcViable #Adds values to growthRates
  growthRates[["growthRatesTRMaxconcViableCellVolume"]]<-growthRatesTRMaxconcViableCellVolume #Adds values to growthRates
  growthRates[["cigrowthRatesTRMaxconcViable"]]<-cigrowthRatesTRMaxconcViable #Adds values to growthRates
  growthRates[["cigrowthRatesTRMaxconcViableCellVolume"]]<-cigrowthRatesTRMaxconcViableCellVolume #Adds values to growthRates
  return(growthRates) #Returns growthRates-compiled info
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; avgViabDiam; dataTransformation; getDF; calculateGrowthCorr or calculateProdCorr
#Does: Writes data into a text-file
#Use: For later calculations and in-between-save
########################################################################################################################################
outputVariable<-function(variable,parameter,path){
    sink(paste(path,"/",parameter,".txt",sep="")) #Writes Variable info into a txt-file
    print(variable)
    sink()
}
########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; avgViabDiam; dataTransformation; getDF; calculateGrowthCorr
#Does: Creates barcharts of different growth calculations
#Use: For growth representations and comparisons
########################################################################################################################################
drawGrowthCharts<-function(growthRates,path){
  w<-getOption('warn') #Get warning option
  options(warn=-1) #Suppresses warnings (arrows sometimes not drawn, because of 0 length)
  drawBarChart(variable=growthRates,parameter="growthRatesconcViable",path=path,ylab="µ [h^-1]") #Executes function that draws barchart of all Growth rates calculated based on ln(VCD)vstime
  drawBarChart(variable=growthRates,parameter="growthRatesconcViableCellVolume",path=path,ylab="µ [h^-1]") #Executes function that draws barchart of all Growth rates calculated based on ln(VCV)vstime
  drawBarChart(variable=growthRates,parameter="growthRatesTRMaxconcViable",path=path,ylab="µ [h^-1]") #Executes function that draws barchart of Growth rates from median time of RMax based on ln(VCD)vstime
  drawBarChart(variable=growthRates,parameter="growthRatesTRMaxconcViableCellVolume",path=path,ylab="µ [h^-1]") #Executes function that draws barchart of Growth rates from median time of RMax based on ln(VCV)vstime
  drawBarChart(variable=growthRates,parameter="growthRatesRMaxconcViable",path=path,ylab="µ [h^-1]") #Executes function that draws barchart of Growth rates from maximum correlation of each sample based on ln(VCD)vstime
  drawBarChart(variable=growthRates,parameter="growthRatesRMaxconcViableCellVolume",path=path,ylab="µ [h^-1]") #Executes function that draws barchart of Growth rates from maximum correlation of each sample based on ln(VCV)vstime
  options(warn=w) #Put in old warning option
  png(paste(path,"/tRMaxvsRMax.png",sep="")) #Creates a box plot of different methods to select growth rates
    boxplot(list(growthRates$growthRatesTRMaxconcViable,growthRates$growthRatesRMaxconcViable,growthRates$growthRatesTRMaxconcViableCellVolume,growthRates$growthRatesRMaxconcViableCellVolume),names=c("tRMaxVCD","RMaxVCD","tRMaxVCV","RMaxVCV"),main="Growth rate variation",xlab="Calculation method",ylab="µ [h^-1]")
  dev.off()
}


########################################################################################################################################
#Creator: GK
#Dependencies: drawGrowthCharts or drawProdCharts
#Does: Creates barcharts of different growth and specific productivity calculations
#Use: to draw bar charts for growth and productivity-data
#########################################################################################################################################
drawBarChart<-function(variable,parameter,path,ylab){
  if(length(nrow(variable[[paste(parameter)]]))>0){ #Determines which png-file creation to use, based on parameter type (vector or matrix)
    width<-nrow(variable[[paste(parameter)]])*100 #dynamic width of png-file
    if(width>10000){
      width<-10000 #Maximum width of png-file
    }
  } else {
    width<-length(variable[[paste(parameter)]])*100#dynamic width of png-file
    if(width>10000){
      width<-10000 #Maximum width of png-file
    }
  }
  png(paste(path,"/",parameter,".png",sep=""),width=width,height=400) #Creates png-file
  barcenters<-barplot(t(variable[[paste(parameter)]]),ylim=c(0,1.15*max(variable[[paste(parameter)]]+variable[[paste("ci",parameter,sep="")]],na.rm=TRUE)), beside=TRUE,legend.text=TRUE,args.legend=list(x="topright",cex=0.7),ylab=ylab) #Draws barchart
  arrows(barcenters,t(variable[[paste(parameter)]])-t(variable[[paste("ci",parameter,sep="")]]), barcenters,t(variable[[paste(parameter)]])+t(variable[[paste("ci",parameter,sep="")]]), lwd = 1.5, angle = 90,code = 3, length = 0.05) #Adds Confidence interval
  dev.off()
}


########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; getDF;
#Does: Assigns titer to ViCell-samples based on the sampleName
#Use: For titer and productivity-representations and comparisons
########################################################################################################################################
assignTiter<-function(df,inputFile,outputPath,col=NA,seedTime=NA,perRow=NA){
  prod<-read.table(inputFile,header=FALSE,sep="\t") #Reads productivity-File into R
  titer<-matrix(NA,nrow=nrow(df),ncol=ncol(prod)-1)
  for(i in 1:length(df[,"sampleName"])){
    if(length(prod[prod[,1] %in% df[,"sampleName"] [i],2])>0){ #If sampleName is found in productivity-file...
      for(j in 2:ncol(prod)){
        titer[i,(j-1)]<-prod[prod[,1] %in% df[,"sampleName"] [i],j] #...the productivity-value is transferred to the same position in titer
      }
    }else{
      titer[i,(1:(ncol(prod)-1))]<-NA #...otherwise, an NA is added
    }
  }
  if(is.na(seedTime)){ # Calculates time post seed in hours, and also determines maximum x-axis-values for all graphs
    xmaxTime<-(as.numeric(max(df[,"runDate"]))-as.numeric(min(df[,"runDate"])))/3600
  } else { #For offset between seeding and first measurement
    xmaxTime<-(as.numeric(max(df[,"runDate"]))-as.numeric(as.POSIXct(seedTime,origin='1970-01-01 00:00:00 UTC')))/3600
  }
  if(is.na(col)){ #Defines colours for plots
    col<-1:length(unique(df[,"replicateNumber"]))
  }
  temp<-colnames(df)
  df<-cbind(df,titer) #Adds titer(s) to dataframe
  colnames(df)<-c(temp,paste("titer",1:(ncol(prod)-1),sep="-"))
  drawParvsTime(df=df,parameter=paste("titer",1:(ncol(prod)-1),sep="-"),col=col,path=outputPath,xmaxTime=xmaxTime,ylabel="Titer [mg/l]",legend="topleft",suppCI=TRUE,perRow=perRow) #Executes function that draws titervstime
  return(df) #Return dataframe with added Titer
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; getDF;
#Does: Calculates Cumulative cell days, cumulative cell volumes, and the sums of it, and adds it to the dataframe
#Use: For calculation of specific productivities
#Missing: Confidence interval, has to be simulated, according to Michael Melcher
########################################################################################################################################
calculateCumulative<-function(df,outputPath=NA,seedTime=NA,col=NA,perRow=NA){
  ccd<-c(rep(NA,times=length(df[,1]))) #Initialisation - ccd (cumulative cell density)
  ccv<-c(rep(NA,times=length(df[,1]))) #Initialisation - ccv (cumulative cell volume)
  sumccd<-c(rep(NA,times=length(df[,1]))) #Initialisation - sumccd (sum of cumulative cell density)
  sumccv<-c(rep(NA,times=length(df[,1]))) #Initialisation - sumccv (sum of cumulative cell volume)
  for(i in unique(df[,"replicateID"])){
    temp<-df[df[,"replicateID"]==i,] #Subsample by replicateID
    for(j in unique(temp[,"replicateNumber"])){
      temp2<-temp[temp[,"replicateNumber"]==j,] #Subsample by replicateNumber
      timepoints<-c()
      timepoints<-temp2[,"TPNumber"] # Extracts TPNumber of subsample
      if(length(timepoints)>1){
        for(k in 2:length(timepoints)){
          ccd[which(df[,"replicateID"]==i & df[,"replicateNumber"]==j & df[,"TPNumber"]==timepoints[k])]<-(temp2[temp2[,"TPNumber"]==timepoints[k],"concViable"]-temp2[temp2[,"TPNumber"]==timepoints[k-1],"concViable"])*1E6/(temp2[temp2[,"TPNumber"]==timepoints[k],"lnconcViable"]-temp2[temp2[,"TPNumber"]==timepoints[k-1],"lnconcViable"])*((as.numeric(temp2[temp2[,"TPNumber"]==timepoints[k],"runDate"])-as.numeric(temp2[temp2[,"TPNumber"]==timepoints[k-1],"runDate"]))/(3600*24)) #Calculates ccd based on ccd[k]=(vcd[k]-vcd[k-1])/(ln(vcd[k])-ln(vcd[k-1]))*(t[k]-t[k-1])
          ccv[which(df[,"replicateID"]==i & df[,"replicateNumber"]==j & df[,"TPNumber"]==timepoints[k])]<-(temp2[temp2[,"TPNumber"]==timepoints[k],"concViableCellVolume"]-temp2[temp2[,"TPNumber"]==timepoints[k-1],"concViableCellVolume"])/(temp2[temp2[,"TPNumber"]==timepoints[k],"lnconcViableCellVolume"]-temp2[temp2[,"TPNumber"]==timepoints[k-1],"lnconcViableCellVolume"])*((as.numeric(temp2[temp2[,"TPNumber"]==timepoints[k],"runDate"])-as.numeric(temp2[temp2[,"TPNumber"]==timepoints[k-1],"runDate"]))/(3600*24)) #Calculates ccv based ccv[k]=(vcv[k]-vcv[k-1])/(ln(vcv[k])-ln(vcv[k-1]))*(t[k]-t[k-1])
        }
        for(k in 2:length(timepoints)){
          sumccd[which(df[,"replicateID"]==i & df[,"replicateNumber"]==j & df[,"TPNumber"]==timepoints[k])]<-sum(ccd[which(df[,"replicateID"]==i & df[,"replicateNumber"]==j & df[,"TPNumber"] %in% timepoints[2:k])],na.rm=TRUE) #Sums up ccd
          sumccv[which(df[,"replicateID"]==i & df[,"replicateNumber"]==j & df[,"TPNumber"]==timepoints[k])]<-sum(ccv[which(df[,"replicateID"]==i & df[,"replicateNumber"]==j & df[,"TPNumber"] %in% timepoints[2:k])],na.rm=TRUE) #Sums up ccv
        }
      }
    }
  }
  df<-cbind(df,ccd,ccv,sumccd,sumccv) #Adds new parameters to dataframe
  if(is.na(seedTime)){ # Calculates time post seed in hours, and also determines maximum x-axis-values for all graphs
    xmaxTime<-(as.numeric(max(df[,"runDate"]))-as.numeric(min(df[,"runDate"])))/3600
  } else { #For offset between seeding and first measurement
    xmaxTime<-(as.numeric(max(df[,"runDate"]))-as.numeric(as.POSIXct(seedTime)))/3600
  }
  if(is.na(col)){ #Defines colours for plots
    col<-1:length(unique(df[,"replicateNumber"]))
  }
  drawParvsTime(df=df,parameter="ccd",col=col,path=outputPath,xmaxTime=xmaxTime,ylabel="CCD [cells/ml]",legend="topleft",suppCI=TRUE,perRow=perRow) #Executes function that draws CCDvstime
  drawParvsTime(df=df,parameter="ccv",col=col,path=outputPath,xmaxTime=xmaxTime,ylabel="CCV [m³/ml]",legend="topleft",suppCI=TRUE,perRow=perRow) #Executes function that draws CCVvstime
  drawParvsTime(df=df,parameter="sumccd",col=col,path=outputPath,xmaxTime=xmaxTime,ylabel="Sum of CCD [cells/ml]",legend="topleft",suppCI=TRUE,perRow=perRow) #Executes function that draws sum of CCDvstime
  drawParvsTime(df=df,parameter="sumccv",col=col,path=outputPath,xmaxTime=xmaxTime,ylabel="Sum of CCV [m³/ml]",legend="topleft",suppCI=TRUE,perRow=perRow) #Executes function that draws sum of CCVvstime
  return(df) #Returns df + cumulative values
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; getDF; calculateCumulative; assignTiter
#Does: calculates specific productivities based on correlating sum of ccd/ccv and titer
#Use: For calculation of specific productivities
#Missing: Confidence interval-overlap check, because no CI is there for sum of ccd/ccv
########################################################################################################################################
calculateProdCorr<-function(df,nmin=4,corrSeed=FALSE,comb=FALSE,measDelay=0){ #nmin=number of minimum datapoints to be incorporated, measDelay= Number of timepoints post TP01 for first titer measurement
  if(comb==FALSE){ #If replicates should not be combined
    sepreplicateID<-nrow(df[,1]) #New replicateID is initialized
    for(i in 1:nrow(unique(df[,c("replicateID","replicateNumber")]))){
      sepreplicateID[df[,c("replicateID")] %in% unique(df[,c("replicateID","replicateNumber")]) [i,1] & df[,c("replicateNumber")] %in% unique(df[,c("replicateID","replicateNumber")]) [i,2]]<-i #New replicateID is created, based on differences in replicateID and replicateNumber
    }
    for(i in 1:nrow(df)){
      df[i,c("replicateName")]<-paste0(df[i,c("replicateName","replicateNumber")],collapse="-") #New replicateName is created (replicateID+replicateNumber)
    }
    df[,c("replicateID")]<-sepreplicateID #New replicateID replaces old one
  } else{
    #Check whether samples overlap in their confidence interval
    if(corrSeed==TRUE & max(df[,"replicateNumber"])>1){ #if replicates should be combined, and seed-corrected
      corrFactorsumccd<-c() #Initializes correction factors for ln(VCD)
      corrFactorsumccv<-c()  #Initializes correction factors for ln(VCV)
      for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID separately
        for(j in 2:max(df[,"replicateNumber"])){ #For each replicateNumber except one
          corrFactorsumccd<-df[df[,"replicateID"]==i & df[,"TPNumber"]==2 & df[,"replicateNumber"]==j,"sumccd"] - df[df[,"replicateID"]==i & df[,"TPNumber"]==2 & df[,"replicateNumber"]==1,"sumccd"] #Difference between sumccd of replicate j and replicate 1 at TP02 (TP01-->no ccd) -->Correction factor
          corrFactorsumccv<-df[df[,"replicateID"]==i & df[,"TPNumber"]==2 & df[,"replicateNumber"]==j,"sumccv"] - df[df[,"replicateID"]==i & df[,"TPNumber"]==2 & df[,"replicateNumber"]==1,"sumccv"] #Difference between sumccv of replicate j and replicate 1 at TP02 (TP01-->no ccv) -->Correction factor
          df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"sumccd"]<-df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"sumccd"]-corrFactorsumccd #Apply correction factor by substraction
          df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"sumccv"]<-df[df[,"replicateID"]==i & df[,"replicateNumber"]==j,"sumccv"]-corrFactorsumccv  #Apply correction factor by substraction
        }
      }
    }
  }

  df<-df[which(apply(as.matrix(df[,grep("titer",colnames(df))]),1,sum,na.rm=TRUE)>0),] #All data without titer is withdrawn
  specProd<-list()#Initialize variables - will contain growth data
  rsumccd<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"]))))  #Initialize variable - Correlation coefficients between sumccd and titer
  rsumccv<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])))) #Initialize variable - Correlation coefficients between sumccv and titer
  rMaxsumccd<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Maximum Correlation coefficients between sumccd and titer
  rMaxsumccv<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Maximum Correlation coefficients between sumccv and titer
  specProdsumccd<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])))) #Initialize variable - Specific productivities based on slope between sumccd and titer
  specProdsumccv<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])))) #Initialize variable - Specific productivities based on slope between sumccv and titer
  specProdRMaxsumccd<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Specific Productivities based on slope between sumccd and titer at sampling point with maximum Correlation coefficients
  specProdRMaxsumccv<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Specific Productivities based on slope between sumccv and titer at sampling point with maximum Correlation coefficients
  specProdTRMaxsumccd<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Specific Productivities based on slope between sumccd and titer at sampling point with median maximum Correlation coefficients
  specProdTRMaxsumccv<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Specific Productivities based on slope between sumccv and titer at sampling point with median maximum Correlation coefficients
  cispecProdsumccd<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])))) #Initialize variable - Confidence interval of specific productivities based on slope between sumccd and titer
  cispecProdsumccv<-matrix(NA,nrow=length(unique(df[,"replicateID"])),ncol=(length(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"]))))   #Initialize variable - Confidence interval of specific productivities based on slope between sumccv and titer
  cispecProdRMaxsumccd<-c(rep(NA,length(unique(df[,"replicateID"]))))   #Initialize variable - Confidence interval of specific productivities based on slope between sumccd and titer at sampling point with maximum Correlation coefficients
  cispecProdRMaxsumccv<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Confidence interval of specific productivities based on slope between sumccv and titer at sampling point with maximum Correlation coefficients
  cispecProdTRMaxsumccd<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Confidence interval of specific productivities based on slope between sumccd and titer at sampling point with median maximum Correlation coefficients
  cispecProdTRMaxsumccv<-c(rep(NA,length(unique(df[,"replicateID"]))))   #Initialize variable - Confidence interval of specific productivities based on slope between sumccv and titer at sampling point with median maximum Correlation coefficients
  tRMaxsumccd<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Timepoints of best correlation between sumccd and titer
  tRMaxsumccv<-c(rep(NA,length(unique(df[,"replicateID"])))) #Initialize variable - Timepoints of best correlation between sumccv and titer
  rownames(rsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(rsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(specProdsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(specProdsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(cispecProdsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  rownames(cispecProdsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(rMaxsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(rMaxsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(specProdRMaxsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(specProdRMaxsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(specProdTRMaxsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(specProdTRMaxsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cispecProdRMaxsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cispecProdRMaxsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cispecProdTRMaxsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(cispecProdTRMaxsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(tRMaxsumccd)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  names(tRMaxsumccv)<-unique(df[,c("replicateName")]) #Add names to matrices/vectors
  colnames(rsumccd)<-paste("TP",(1+measDelay),"-",sort(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])),sep="") #Add names to matrices/vectors
  colnames(rsumccv)<-paste("TP",(1+measDelay),"-",sort(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])),sep="") #Add names to matrices/vectors
  colnames(specProdsumccd)<-paste("TP",(1+measDelay),"-",sort(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])),sep="") #Add names to matrices/vectors
  colnames(specProdsumccv)<-paste("TP",(1+measDelay),"-",sort(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])),sep="") #Add names to matrices/vectors
  colnames(cispecProdsumccd)<-paste("TP",(1+measDelay),"-",sort(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])),sep="") #Add names to matrices/vectors
  colnames(cispecProdsumccv)<-paste("TP",(1+measDelay),"-",sort(unique(df[df[,"TPNumber"]>(1+measDelay),"TPNumber"])),sep="") #Add names to matrices/vectors
  for(n in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    temp<-df[df[,"replicateID"]==unique(df[,"replicateID"]) [n],] #Subsample
    for(i in unique(temp[temp[,"TPNumber"]>(1+measDelay),"TPNumber"])){ #(Tp1+measDelay)-Tpi (TP1 has no CCD/CCV)
      rsumccd[n,(i-(1+measDelay))]<-summary(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccd"],length(grep("titer",colnames(temp))))))$r.squared #determines correlation coefficients of linear correlations between sum of ccd and titer
      rsumccv[n,(i-(1+measDelay))]<-summary(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccv"],length(grep("titer",colnames(temp))))))$r.squared #determines correlation coefficients of linear correlations between sum of ccv and titer
      #specificProductivites
      specProdsumccd[n,(i-(1+measDelay))]<-coef(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccd"],length(grep("titer",colnames(temp))))))[2] *1e6 #determines slope (=pg/(cell*d)) of linear correlations between sum of ccd and titer
      specProdsumccv[n,(i-(1+measDelay))]<-coef(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccv"],length(grep("titer",colnames(temp))))))[2]/1000 #determines slope (=pg/(cm³*d)) of linear correlations between sum of ccv and titer
      w<-getOption('warn') #Gets warning-option
      options(warn=-1) #Suppresses warnings, because next cmd can generate a warning, if only two values are available
      cispecProdsumccd[n,(i-(1+measDelay))]<-(confint(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccd"],length(grep("titer",colnames(temp))))))[2,2]-coef(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccd"],length(grep("titer",colnames(temp))))))[2])*1e6 #calculates confidence interval of linear correlation between sumccd and titer
      if(is.na(cispecProdsumccd[n,(i-(1+measDelay))])){ #If NA is generated for a confidence interval, it is replaced by a 0
        cispecProdsumccd[n,(i-(1+measDelay))]<-0
      }
      cispecProdsumccv[n,(i-(1+measDelay))]<-(confint(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccv"],length(grep("titer",colnames(temp))))))[2,2]-coef(lm(c(unlist(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),grep("titer",colnames(temp))]))~rep(temp[temp[,"TPNumber"] %in% c((1+measDelay):i),"sumccv"],length(grep("titer",colnames(temp))))))[2])/1000 #calculates confidence interval of linear correlation between sumccv and titer
      if(is.na(cispecProdsumccv[n,(i-(1+measDelay))])){ #If NA is generated for a confidence interval, it is replaced by a 0
        cispecProdsumccv[n,(i-(1+measDelay))]<-0
      }
      options(warn=w) #Warnings are turned back on
    }
  }

  for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    if(sum(!is.na(rsumccd[i,c((nmin-1):length(rsumccd[1,]))]))>0){ #If rsumccd not only consists of NA for at least nmin (-1 because it starts with 2 tp) timepoints
      rMaxsumccd[i]<-which(rsumccd[i,]==max(rsumccd[i,c((nmin-1):length(rsumccd[1,]))],na.rm=TRUE)) #Maximum Correlation coefficents between sum of ccd and titer are extracted (at least nmin timepoints!)
    }else{ #If no correlation coefficients are stored, an NA is added
      rMaxsumccd[i]<-NA
    }
    if(sum(!is.na(rsumccv[i,c((nmin-1):length(rsumccv[1,]))]))>0){ #If rsumccv not only consists of NA for at least nmin (-1 because it starts with 2 tp) timepoints
      rMaxsumccv[i]<-which(rsumccv[i,]==max(rsumccv[i,c((nmin-1):length(rsumccv[1,]))],na.rm=TRUE)) #Maximum Correlation coefficents between sum of ccv and titer are extracted (at least nmin timepoints!)
    }else{#If no correlation coefficients are stored, an NA is added
      rMaxsumccv[i]<-NA
    }
  }
  for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    tRMaxsumccd[i]<-mean(as.numeric(df[df[,"replicateID"]==unique(df[,"replicateID"]) [i] & df[,"TPNumber"]==(rMaxsumccd[i]+(1+measDelay)),c("runDate")])) #Timepoints of best correlations between sum of ccd and titer are extracted
    tRMaxsumccv[i]<-mean(as.numeric(df[df[,"replicateID"]==unique(df[,"replicateID"]) [i] & df[,"TPNumber"]==(rMaxsumccv[i]+(1+measDelay)),c("runDate")])) #Timepoints of best correlations between sum of ccv and titer are extracted
  }
  #Determine specific productivities based on maximum correlation coefficients
  for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    specProdRMaxsumccd[i]<-specProdsumccd[i,rMaxsumccd[i]] #sum of ccd-based specific productivities based on maximum correlation coefficient are extracted
    specProdRMaxsumccv[i]<-specProdsumccv[i,rMaxsumccv[i]] #sum of ccv-based specific productivities based on maximum correlation coefficient are extracted
    cispecProdRMaxsumccd[i]<-cispecProdsumccd[i,rMaxsumccd[i]] #sum of ccd-based confidence interval based on maximum correlation coefficient are extracted
    cispecProdRMaxsumccv[i]<-cispecProdsumccv[i,rMaxsumccv[i]] #sum of ccv-based confidence interval based on maximum correlation coefficient are extracted
  }
  temp2<-c() #Clears subsample
  for(i in 1:length(unique(df[,"replicateID"]))){ #For each replicateID
    for(j in (nmin+measDelay):max(df[,"TPNumber"],na.rm=TRUE)){ #For each TPNumber
      temp2[j]<-mean(as.numeric(df[df[,"replicateID"]==unique(df[,"replicateID"]) [i] & df[,"TPNumber"]==j,"runDate"])) #Extracts average (if more than one sample per timepoint) runDate for each TP
    }
    if(length(temp2[!is.na(temp2)])>0){
      specProdTRMaxsumccd[i]<-specProdsumccd[i,which(abs(median(tRMaxsumccd,na.rm=TRUE)-temp2)==min(abs(median(tRMaxsumccd,na.rm=TRUE)-temp2),na.rm=TRUE))-(1+measDelay)] #Extracts specific productivity based on sum of ccd vs titer of runDate which is closest to the median runDate of max. Correlation Coefficient (-1+measDelay, dep on first measurement of Productivity)
      specProdTRMaxsumccv[i]<-specProdsumccv[i,which(abs(median(tRMaxsumccv,na.rm=TRUE)-temp2)==min(abs(median(tRMaxsumccv,na.rm=TRUE)-temp2),na.rm=TRUE))-(1+measDelay)] #Extracts specific productivity based on sum of ccv vs titer of runDate which is closest to the median runDate of max. Correlation Coefficient (-1+measDelay, dep on first measurement of Productivity)
      cispecProdTRMaxsumccd[i]<-cispecProdsumccd[i,which(abs(median(tRMaxsumccd,na.rm=TRUE)-temp2)==min(abs(median(tRMaxsumccd,na.rm=TRUE)-temp2),na.rm=TRUE))-(1+measDelay)] #Extracts confidence interval of specific productivity based on sum of ccd vs titer of runDate which is closest to the median runDate of max. Correlation Coefficient (-1+measDelay, dep on first measurement of Productivity)
      cispecProdTRMaxsumccv[i]<-cispecProdsumccv[i,which(abs(median(tRMaxsumccv,na.rm=TRUE)-temp2)==min(abs(median(tRMaxsumccv,na.rm=TRUE)-temp2),na.rm=TRUE))-(1+measDelay)]#Extracts confidence interval of specific productivity based on sum of ccv vs titer of runDate which is closest to the median runDate of max. Correlation Coefficient (-1+measDelay, dep on first measurement of Productivity)
    }else{
      specProdTRMaxsumccd[i]<-NA
      specProdTRMaxsumccv[i]<-NA
      cispecProdTRMaxsumccd[i]<-NA
      cispecProdTRMaxsumccv[i]<-NA
    }
  }
  specProd[["specProdsumccd"]]<-specProdsumccd #Adds values to specProd
  specProd[["specProdsumccv"]]<-specProdsumccv #Adds values to specProd
  specProd[["cispecProdsumccd"]]<-cispecProdsumccd #Adds values to specProd
  specProd[["cispecProdsumccv"]]<-cispecProdsumccv #Adds values to specProd
  specProd[["rsumccd"]]<-rsumccd #Adds values to specProd
  specProd[["rsumccv"]]<-rsumccv #Adds values to specProd
  specProd[["rMaxsumccd"]]<-rMaxsumccd #Adds values to specProd
  specProd[["rMaxsumccv"]]<-rMaxsumccv #Adds values to specProd
  specProd[["median(rMaxsumccd)"]]<-median(rMaxsumccd,na.rm=TRUE) #Adds values to specProd
  specProd[["median(rMaxsumccv)"]]<-median(rMaxsumccv,na.rm=TRUE) #Adds values to specProd
  specProd[["specProdRMaxsumccd"]]<-specProdRMaxsumccd #Adds values to specProd
  specProd[["specProdRMaxsumccv"]]<-specProdRMaxsumccv #Adds values to specProd
  specProd[["cispecProdRMaxsumccd"]]<-cispecProdRMaxsumccd #Adds values to specProd
  specProd[["cispecProdRMaxsumccv"]]<-cispecProdRMaxsumccv #Adds values to specProd
  specProd[["tRMaxsumccd"]]<-tRMaxsumccd#Adds values to specProd
  specProd[["tRMaxsumccv"]]<-tRMaxsumccv  #Adds values to specProd
  specProd[["median(tRMaxsumccd)"]]<-median(tRMaxsumccd,na.rm=TRUE) #Adds values to specProd
  specProd[["Date of median(tRMaxsumccd)"]]<-as.POSIXct(median(tRMaxsumccd,na.rm=TRUE),origin='1970-01-01 00:00:00 UTC')
  specProd[["median(tRMaxsumccv)"]]<-median(tRMaxsumccv,na.rm=TRUE) #Adds values to specProd
  specProd[["Date of median(tRMaxsumccv)"]]<-as.POSIXct(median(tRMaxsumccv,na.rm=TRUE),origin='1970-01-01 00:00:00 UTC')
  specProd[["specProdTRMaxsumccd"]]<-specProdTRMaxsumccd #Adds values to specProd
  specProd[["specProdTRMaxsumccv"]]<-specProdTRMaxsumccv #Adds values to specProd
  specProd[["cispecProdTRMaxsumccd"]]<-cispecProdTRMaxsumccd #Adds values to specProd
  specProd[["cispecProdTRMaxsumccv"]]<-cispecProdTRMaxsumccv #Adds values to specProd
  return(specProd) #Returns specProd-compiled info
}

########################################################################################################################################
#Creator: GK
#Dependencies: readTxt; avgViabDiam; dataTransformation; getDF; calculateProdCorr
#Does: Creates barcharts of different specific productivity calculations
#Use: For productivity representations and comparisons
########################################################################################################################################
drawProdCharts<-function(df,specProd,path,col=NA,perRow=NA){
  if(is.na(col)){ #Defines colours for plots
    col<-1:length(unique(df[,"replicateNumber"]))
  }
  w<-getOption('warn') #Get warning option
  options(warn=-1) #Suppresses warnings (arrows sometimes not drawn, because of 0 length)
  drawBarChart(variable=specProd,parameter="specProdsumccd",path=path,ylab="qp [pg/(cell*d)]") #Executes function that draws barchart of all specific productivities calculated based on sum of CCD vs titer
  drawBarChart(variable=specProd,parameter="specProdsumccv",path=path,ylab="qp [pg/(cm³*d)]") #Executes function that draws barchart of all specific productivities calculated based on sum of CCV vs titer
  drawBarChart(variable=specProd,parameter="specProdRMaxsumccd",path=path,ylab="qp [pg/(cell*d)]") #Executes function that draws barchart of specific productivities from median time of RMax based on sum of CCD vs titer
  drawBarChart(variable=specProd,parameter="specProdRMaxsumccv",path=path,ylab="qp [pg/(cm³*d)]") #Executes function that draws barchart of specific productivities from median time of RMax based on sum of CCV vs titer
  drawBarChart(variable=specProd,parameter="specProdTRMaxsumccd",path=path,ylab="qp [pg/(cell*d)]") #Executes function that draws barchart of specific productivities from maximum correlation of each sample based on sum of CCD vs titer
  drawBarChart(variable=specProd,parameter="specProdTRMaxsumccv",path=path,ylab="qp [pg/(cm³*d)]") #Executes function that draws barchart of specific productivities from maximum correlation of each sample based on sum of CCV vs titer
  options(warn=w) #Put in old warning option
  pdf(paste(path,"/tRMaxvsRMax.pdf",sep="")) #Creates a box plot of different methods to select specific productivities
  boxplot(list(specProd$specProdTRMaxsumccd,specProd$specProdRMaxsumccd),names=c("tRMaxCCD","RMaxCCD"),main="Productivity variation",xlab="Calculation method",ylab="qp [pg/(cell*d)]")
  boxplot(list(specProd$specProdTRMaxsumccv,specProd$specProdRMaxsumccv),names=c("tRMaxCCV","RMaxCCV"),main="Productivity variation",xlab="Calculation method",ylab="qp [pg/(cm³*d)]")
  dev.off()
  #TitervsCCD/CCV
  df<-df[which(apply(as.matrix(df[,grep("titer",colnames(df))]),1,sum,na.rm=TRUE)>0),]
  ymax<-max(df[,grep("titer",colnames(df))],na.rm=TRUE) #Determines max y-axis-value
  ymin<-min(df[,grep("titer",colnames(df))],na.rm=TRUE) #Determines min y-axis-value
  xmaxd<-max(df[,"sumccd"],na.rm=TRUE) #Determines max x-axis-value
  xmind<-min(df[,"sumccd"],na.rm=TRUE) #Determines min x-axis-value
  xmaxv<-max(df[,"sumccv"],na.rm=TRUE) #Determines max x-axis-value
  xminv<-min(df[,"sumccv"],na.rm=TRUE) #Determines min x-axis-value
  if(!is.na(perRow)){
    png(paste(path,"/","titervsccd.png",sep=""),width=400*perRow,height=(400*length(unique(df[,"replicateID"]))/perRow))
    par(mfrow=c(ceiling(length(unique(df[,"replicateID"]))/perRow),perRow)) #Determines number of plots
  }else if(length(unique(df[,"replicateID"]))>1){
    png(paste(path,"/","titervsccd.png",sep=""),width=800,height=(400*length(unique(df[,"replicateID"]))/2)) #Creates png-file
    par(mfrow=c(ceiling(length(unique(df[,"replicateID"]))/2),2)) #Determines number of plots
  } else{
    png(paste(path,"/","titervsccd.png",sep=""),width=800,height=400) #Creates png-file
  }
  for(i in unique(df[,"replicateID"])){ #Creates plot for each replicate ID
    par(new=FALSE)
    temp<-df[df[,"replicateID"]==i,] #Subsample
    gtitle<-temp[1,c("replicateName")] # Title selection for each sample
    for(j in unique(temp[,"replicateNumber"])){
      temp2<-temp[temp[,"replicateNumber"]==j,] #Subsample
      for(k in 1:length(grep("titer",colnames(df)))){
        plot(temp2[,c("sumccd")],temp2[,paste("titer",k,sep="-")],xlab="sum of cumulative cell days [cells*d]",col=col[j],ylab="titer [µg/ml]",main=gtitle,ylim=c(ymin,ymax),xlim=c(xmind,xmaxd),pch=20,cex.axis=1.5,cex.lab=1.5,cex.main=2.5)# Plots data
        lines(temp2[,c("sumccd")],temp2[,paste("titer",k,sep="-")],type="l",col=col[j]) #Adds lines between points
        par(new=TRUE)#So that replicate will be added into same plot
      }
    }
    legend("topright",legend=paste(sort(unique(temp[,"replicateNumber"]))),pch=20,cex=1.5,col=col) #Adds replicate legend
  }
  dev.off() #Closes file
  if(!is.na(perRow)){
    png(paste(path,"/","titervsccv.png",sep=""),width=400*perRow,height=(400*length(unique(df[,"replicateID"]))/perRow))
    par(mfrow=c(ceiling(length(unique(df[,"replicateID"]))/perRow),perRow)) #Determines number of plots
  }else if(length(unique(df[,"replicateID"]))>1){
    png(paste(path,"/","titervsccv.png",sep=""),width=800,height=(400*length(unique(df[,"replicateID"]))/2)) #Creates png-file
    par(mfrow=c(ceiling(length(unique(df[,"replicateID"]))/2),2)) #Determines number of plots
  } else{
    png(paste(path,"/","titervsccv.png",sep=""),width=800,height=400) #Creates png-file
  }
  for(i in unique(df[,"replicateID"])){ #Creates plot for each replicate ID
    par(new=FALSE)
    temp<-df[df[,"replicateID"]==i,] #Subsample
    gtitle<-temp[1,c("replicateName")] # Title selection for each sample
    for(j in unique(temp[,"replicateNumber"])){
      temp2<-temp[temp[,"replicateNumber"]==j,] #Subsample
      for(k in 1:length(grep("titer",colnames(df)))){
        plot(temp2[,c("sumccv")],temp2[,paste("titer",k,sep="-")],xlab="sum of cumulative cell volume [m³*d]",col=col[j],ylab="titer [µg/ml]",main=gtitle,ylim=c(ymin,ymax),xlim=c(xminv,xmaxv),pch=20,cex.axis=1.5,cex.lab=1.5,cex.main=2.5)# Plots data
        lines(temp2[,c("sumccv")],temp2[,paste("titer",k,sep="-")],type="l",col=col[j]) #Adds lines between points
        par(new=TRUE)#So that replicate will be added into same plot
      }
    }
    legend("topright",legend=paste(sort(unique(temp[,"replicateNumber"]))),pch=20,cex=1.5,col=col) #Adds replicate legend
  }
  dev.off() #Closes file
}


