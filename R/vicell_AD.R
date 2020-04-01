# give either a filename
# file = "filename.txt"
# or a directory containing only (!) vicell txt files
# dir = "path/to/directory"
# dir is either relative to the working directory, or an absolute path
# timezone is an optional parameter
readTxt <- function(file = F, dir = F, extended=T, timeZone = "Europe/Vienna", purgeDate = F, sampleIDSearch=c(), sampleIDReplace=c()) {
  if(dir != F) {
    files <- list.files(dir)
  } else if(file != F) {
    files <- c(file)
  } else {
    return(FALSE)
  }

  data <- list()
  fileNum <- 0
  for (file in files) {
    fileData <- list()

    conn <- file(paste(dir, "/", file, sep=""), open="r")
    linn <- readLines(conn)

    # check if this file is a proper vicell file by looking for the CS entry in the last row
    if(!grepl("^CS : ", linn[length(linn)])) {
      warning(paste("Skipped file ",file, " - not a proper ViCell file (missing check sum row)!\n", sep=""))
      close(conn)
      next
    }

    for (i in 1:length(linn)){
      if(grepl("^Sample ID : ", linn[i])) {
        # sampleID
        str <- gsub(".* : ", "", linn[i])
        sampleName <- str
        if(purgeDate == T) {
          sampleName <- gsub("^\\d{1,8}[-_]", "", sampleName)
        }
        if(length(sampleIDReplace) > 0 && length(sampleIDSearch) == length(sampleIDReplace)) {
          for (rep in 1:length(sampleIDSearch)) {
            sampleName <- gsub(sampleIDSearch[rep], sampleIDReplace[rep], sampleName)
          }
          #mgsub does NOT behave like gsub - so we do it ourselfs with the for loop
          #sampleName <- mgsub(sampleIDSearch, sampleIDReplace, sampleName)
        }
        fileData[["sampleID"]] <- str
        fileData[["sampleName"]] <- sampleName
      } else if (grepl("^RunDate : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["runDate"]] <- lubridate::parse_date_time(str, "%d %b %Y  %I:%M:%S %p", tz=timeZone)
        fileData[["runDay"]] <- as.character(fileData[["runDate"]], format="%Y-%m-%d")
      } else if (grepl("^Images : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["images"]] <- strtoi(str)
      } else if (grepl("^Total cells : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["totalCells"]] <- strtoi(str)
      } else if (grepl("^Viable cells : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["viableCells"]] <- strtoi(str)
      } else if (grepl("^Nonviable cells : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["nonviableCells"]] <- strtoi(str)
      } else if (grepl("^Viability \\(\\%\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["viability"]] <- as.double(str)
      } else if (grepl("^Total cells \\/ ml \\(x 10\\^6\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["concTotal"]] <- as.double(str) #*10^6
      } else if (grepl("^Total viable cells \\/ ml \\(x 10\\^6\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["concViable"]] <- as.double(str) #*10^6
      } else if (grepl("^Average diameter \\(microns\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["avgDiameter"]] <- as.double(str)
      } else if (grepl("^Average circularity : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["avgCircularity"]] <- as.double(str)
      } else if (grepl("^Average cells \\/ image : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["avgCellsImage"]] <- as.double(str)
      } else if (grepl("^Microns\\/pixel ratio : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["micronsPixelRatio"]] <- as.double(str)
      } else if (grepl("^Total diameter sum : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["totalDiameterSum"]] <- as.double(str)
      } else if (grepl("^Total circularity sum : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["totalCircularitySum"]] <- as.double(str)
      } else if (grepl("^Background intensity sum : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["bgIntensitySum"]] <- as.double(str)
      } else if (grepl("^Cell type : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setCellType"]] <- as.character(str)
      } else if (grepl("^Threshold \\(\\%\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setThreshold"]] <- as.double(str)
      } else if (grepl("^Center threshold \\(\\%\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setCenterThreshold"]] <- as.double(str)
      } else if (grepl("^Minimum diameter \\(microns\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setMinDiameter"]] <- as.double(str)
      } else if (grepl("^Maximum diameter \\(microns\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setMaxDiameter"]] <- as.double(str)
      } else if (grepl("^Minimum center size : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setMinCenterSize"]] <- as.double(str)
      } else if (grepl("^Minimum circularity : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setMinCircularity"]] <- as.double(str)
      } else if (grepl("^Frames : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setFrames"]] <- as.double(str)
      } else if (grepl("^Dilution : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setDilution"]] <- as.double(str)
      } else if (grepl("^Focus parameter : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setFocusParameter"]] <- as.double(str)
      } else if (grepl("^Sample flush cycles : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setSampleFlushCycles"]] <- as.double(str)
      } else if (grepl("^Trypan blue mixing cycles : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setTrypanBlueMixingCycles"]] <- as.double(str)
      } else if (grepl("^Internal Dilution : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setInternalDilution"]] <- as.double(str)
      } else if (grepl("^Decluster degree : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setDeclusterDegree"]] <- as.double(str)
      } else if (grepl("^Number of bins : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setNumberOfBins"]] <- as.double(str)
      } else if (grepl("^Field of view \\(microns\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setFieldOfView"]] <- as.character(str)
      } else if (grepl("^Sample depth \\(microns\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setSampleDepth"]] <- as.double(str)
      } else if (grepl("^Probe volume \\(ml x 10\\^\\-6\\) : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setProbeVolume"]] <- as.double(str)
      } else if (grepl("^Image size : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setImageSize"]] <- as.character(str)
      } else if (grepl("^Comment : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["setComment"]] <- as.character(str)
      } else if (grepl("^SizeData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["sizeData"]] <- vector
      } else if (grepl("^ViableSizeData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["viableSizeData"]] <- vector
      } else if (grepl("^CircData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["circData"]] <- vector
      } else if (grepl("^ViableCircData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["viableCircData"]] <- vector
      } else if (grepl("^ViabilityData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["viabilityData"]] <- vector
      } else if (grepl("^CountData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["countData"]] <- vector
      } else if (grepl("^ViableCellsData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["viableCellsData"]] <- vector
      } else if (grepl("^TotalPerMlData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["totalPerMlData"]] <- vector
      } else if (grepl("^ViablePerMlData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["viablePerMlData"]] <- vector
      } else if (grepl("^AvgDiamData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["avgDiamData"]] <- vector
      } else if (grepl("^AvgCircData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["avgCircData"]] <- vector
      } else if (grepl("^BackgroundData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["backgroundData"]] <- vector
      } else if (grepl("^ClusterSizeData : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        vector <- as.double(strsplit(str, ",")[[1]])
        fileData[["clusterSizeData"]] <- vector
      } else if (grepl("^CS : ", linn[i])) {
        str <- gsub(".* : ", "", linn[i])
        fileData[["CS"]] <- str
      }
    }
    close(conn)
    fileNum <- fileNum+1;
    data[[fileData[["sampleID"]]]] <- fileData
    cat(paste(file, " loaded\n", sep=""))
  }
  if(!extended) {
    data <- .flattenData(data)
  }
  cat(paste("\nLoaded data from ",fileNum, " files\n", sep=""))
  data
}

exportCSV <- function(data, path) {
  # if data is not yet a data.frame
  if(class(data) == "list") {
    df <- getDF(data)
  } else {
    df <- data
  }
  write.table(df, paste(path,"/df.csv",sep=""), sep="\t",row.names=FALSE)
}

exportXLSX <- function(data, path) {
  # if data is not yet a data.frame
  if(class(data) == "list") {
    df <- getDF(data)
  } else {
    df <- data
  }
  write.xlsx(df, paste(path,"/df.xlsx",sep=""),row.names=FALSE)
}

getDF <- function(data) {
  data <- .flattenData(data)
  df <- do.call(rbind, lapply(data, data.frame, stringsAsFactors=FALSE))
  df
}

.flattenData <- function(data){
  newData <- list()
  for (i in data) {
    for (elm in names(i)) {
      if(length(i[[elm]]) > 1) { # if this is a list/vector/etc
        i[[elm]] <- NULL
      }
    }
    newData[[i[[1]]]] <- i
  }
  data <- newData
}

# calculates specific growth d^-1 and culture age (timediff in hours) for rows with same sampleIDs
# replicates can be grouped by using the replicate column
calculateGrowth <- function(data) {
  oldw <- getOption("warn")
  options(warn = -1)
  data[,c("growth")] <- NA # add growth column
  data[,c("growthTotal")] <- NA # add growth column
  data[,c("timediff")] <- 0 # add timediff column

  timeStart <- min(data[,"runDate"])

  for (i in unique(data$sampleName)) { # iterate over every sampleID present in this data set
    sample <- subset(data, sampleName==i) # subset only data of the current sample
    for (j in unique(sample$replicate)) {
      # unique replicates now
      sampleReplicate <- subset(sample, replicate==j) # subset only data of the current sample
      sampleReplicate$growth <- c(NA, (log(sampleReplicate[2:nrow(sampleReplicate), "concViable"]) - log(sampleReplicate[1:(nrow(sampleReplicate)-1), "concViable"])) / ((as.numeric(sampleReplicate[2:nrow(sampleReplicate), "runDate"]) - as.numeric(sampleReplicate[1:(nrow(sampleReplicate)-1), "runDate"]))/3600/24))
      sampleReplicate$growthTotal <- c(NA, (log(sampleReplicate[2:nrow(sampleReplicate), "concViable"]) - log(sampleReplicate[1, "concViable"])) / ((as.numeric(sampleReplicate[2:nrow(sampleReplicate), "runDate"]) - as.numeric(sampleReplicate[1, "runDate"]))/3600/24))

      sampleReplicate$timediff <- c((as.numeric(sampleReplicate[1:nrow(sampleReplicate), "runDate"]) - as.numeric(timeStart))/3600)
      sample[sample$replicate == j,] <- sampleReplicate # reinsert sample vector with growth data into dataset
    }
    data[data$sampleName == i,] <- sample # reinsert sample vector with growth data into dataset
  }
  options(warn = oldw)
  return(data)
}

drawConcVia <- function(data, save=F, height=F, width=F) {
  oldw <- getOption("warn")
  options(warn = -1)
  y_max <- max(data$concViable)*1.4
  y_min <- min(data$concViable)
  x_max <- 24*ceiling(max(data$timediff)/24)

  plot <- ggplot(data=data,aes(x=timediff, y=concViable), group=sampleName) +
    geom_point(aes(colour=viability), size=4) +
    geom_text(aes(label = paste("[", round(concViable, 2), "] ", round(viability, 0), "%", sep="")),hjust=0.5, vjust=-0.2, alpha=0.5, size=2) +
    scale_size_continuous(range=c(1,4), limits=c(0,100))
  if("group3" %in% colnames(data)) {
    plot <- plot +
      facet_grid(group1 ~ group2 + group3, scales = "free_x", space = "free")
  } else if("group2" %in% colnames(data)) {
    plot <- plot +
      facet_grid(group1 ~ group2, scales = "free_x", space = "free")
  } else if ("group1" %in% colnames(data)) {
    plot <- plot +
      facet_grid(group1 ~ replicate, scales = "free_x", space = "free")
  } else {
    plot <- plot +
      facet_grid(sampleName ~ replicate, scales = "free_x", space = "free")
  }
  plot <- plot +
    scale_colour_gradient(low = "red", high="green", limits=c(50,100)) +
    scale_x_continuous(breaks = seq(min(data$timediff), x_max, by=24), minor_breaks = seq(min(data$timediff), x_max, by=12)) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    geom_smooth(method="loess",
                formula=y ~ x,
                colour="darkgray",
                fill=NA) +
    labs(title=paste("Growth and viability", sep=""), x="Culture age [h]", y="Conc. [10^6 cells/mL]")
  #plot
  if(save) {
    dim <- vicellR::getDimensions(data, height, width)
    ggsave(plot=plot, paste("concentration_viability.pdf", sep = ""), width = dim[1], height = dim[2], units="mm", scale=1, dpi = 300)
  } else {
    return(plot)
  }
  options(warn = oldw)
}


drawGrowth <- function(data, save=F, height=F, width=F) {
  oldw <- getOption("warn")
  options(warn = -1)
  y_max <- max(data$growth)*1.2
  y_min <- min(data$growth)
  x_max <- 24*ceiling(max(data$timediff)/24)

  # bar charts for growth
  growthData <- data[with(data, order(sampleName, runDate)), ]
  growthData <- subset(growthData, !growthTotal=="NA")
  growthData$timeDiffHours <- round(growthData$timediff)
  dodge = position_dodge(width=0.9)
  plot <- ggplot(data=growthData, aes(x=timeDiffHours, y=growthTotal, fill=factor(replicate))) +
    geom_bar(stat = "identity", position=dodge) +
    scale_x_continuous(breaks = seq(min(data$timediff), x_max, by=24), minor_breaks = seq(min(data$timediff), x_max, by=12))
  if("group3" %in% colnames(data)) {
    plot <- plot +
      facet_grid(group1 ~ group2 + group3, scales = "free_x", space = "free")
  } else if("group2" %in% colnames(data)) {
      plot <- plot +
        facet_grid(group1 ~ group2, scales = "free_x", space = "free")
  } else if ("group1" %in% colnames(data)) {
      plot <- plot +
        facet_grid(group1 ~ replicate, scales = "free_x", space = "free")
  } else {
      plot <- plot +
        facet_grid(.~sampleName, scales = "free_x", space = "free")
  }
  plot <- plot +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) +
    labs(title=paste("Total specific growth", sep=""), x="Culture age [h]", y="Total specific growth (t0 - tN) [d^-1]")+
    theme(legend.position="none")

 # ggsave(paste(filename, "_specific_growth_barchart.pdf", sep = ""), width = 300, height = 300, units="mm", scale=0.5, dpi = 300)

  if(save) {
    dim <- vicellR::getDimensions(data, height, width)
    ggsave(plot=plot, paste("total_specific_growth.pdf", sep = ""), width = dim[2], height = dim[1], units="mm", scale=1, dpi = 300)
  }


  # draw total growth
  y_max <- max(data$growthTotal)*1.2
  y_min <- min(data$growthTotal)
  x_max <- 24*ceiling(max(data$timediff)/24)

  # bar charts for growth
  growthData <- data[with(data, order(sampleName, runDate)), ]
  growthData <- subset(growthData, !growth=="NA")
  growthData$timeDiffHours <- round(growthData$timediff)
  dodge = position_dodge(width=0.9)
  plot <- ggplot(data=growthData, aes(x=timeDiffHours, y=growth, fill=factor(replicate))) +
    geom_bar(stat = "identity", position=dodge) +
    #  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
    scale_x_continuous(breaks = seq(min(data$timediff), x_max, by=24), minor_breaks = seq(min(data$timediff), x_max, by=12))
    if("group3" %in% colnames(data)) {
      plot <- plot +
        facet_grid(group1 ~ group2 + group3, scales = "free_x", space = "free")
    } else if("group2" %in% colnames(data)) {
      plot <- plot +
        facet_grid(group1 ~ group2, scales = "free_x", space = "free")
    } else if ("group1" %in% colnames(data)) {
      plot <- plot +
        facet_grid(group1 ~ replicate, scales = "free_x", space = "free")
    } else {
      plot <- plot +
        facet_grid(.~sampleName, scales = "free_x", space = "free")
    }
  plot <- plot +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) +
    labs(title=paste("Specific growth", sep=""), x="Culture age [h]", y="Specific growth [d^-1]")+
    theme(legend.position="none")

  # ggsave(paste(filename, "_specific_growth_barchart.pdf", sep = ""), width = 300, height = 300, units="mm", scale=0.5, dpi = 300)

  if(save) {
    dim <- vicellR::getDimensions(data, height, width)
    ggsave(plot=plot, paste("specific_growth.pdf", sep = ""), width = dim[2], height = dim[1], units="mm", scale=1, dpi = 300)
  }

  options(warn = oldw)
  return(plot)
}

getDimensions <- function(data, height=F, width=F) {
  if(height != F) {
    return(c(width, height))
  }
  samples <- length(unique(data$sampleID))
  width <- 210 # a4 width
  if(samples <= 2) {
    height <- 297/2
  } else if (samples > 6) {
    height <- 297*2
  } else {
    height <- 297
  }
  return(c(width, height))
}
