
### THIS SCRIPT CALIBRATES A HEC-HMS MODEL BY ORDER OF BASIN INFLUENCE (PREDETERIMINED IN BASIN SENSITIVITY ANALYSIS) 
#AND SUBSEQUENTYLY BY ORDER OF PARAMETER INFLUENCE (PREDETERMINED BY PARAMETER SENSITIVITY ANALYSIS)
# Make sure a copy of DSSVue is installed to it's default location for your system: https://www.hec.usace.army.mil/software/hec-dssvue/downloads.aspx
# only need to install once, run following line if haven't installed already
# Dependencies:

install.packages("rJava","http://rforge.net")
install.packages("plyr")
install.packages("broom")
install.packages("rjson")
library("reshape2")
library("stringr")
library("xts")
library("rJava")
library("plyr")
library("broom")
library("rjson")
## devtools::install_github("eheisman/DSS-Rip",args="--no-multiarch")
##     -Instructions from: http://eaheisman.blogspot.com/2013/04/hec-dss-files-and-r-and-python.html
##     -Different dss functions: https://github.com/eheisman/dssrip/commit/f4cf4e3c68d67a1f9a1dc934ace8dbc9a9c9c33a
#############################################################

### IMPORT DEPENDENCIES ###
options(dss_override_location="c:\\programs\\HEC-DSSVue-v3.0.00.212\\")
Sys.setenv(JAVA_HOME=paste0(options("dss_override_location"), "java"))
options(dss_config_filename = "C:\\Users\\kflint\\Documents\\for_dssrip\\jar_config.json")
library(dssrip)



### FUNCTION TO RUN HEC-HMS USING COMMAND LINE (NOT GUI) ###
run.HECHMS <- function(dir.runScript,fname.runScript,dir.hechms) {
  #### Define the command to run the script
  run.command = paste0("hec-hms.cmd -s ",dir.runScript,fname.runScript) 
  #### set working directory to directory to directory of HECHMS.exe 
  setwd(dir.hechms)
  #### initiate the run command
  system(run.command)
}

### FUNCTION TO EXPORT .DSS FILES TO TABLE ###
dssToDataframe <- function(dir.dssVue,dir.model,flow.element = NULL) {
  # DEFINE LOCATIONS OF DSS-VUE AMD FILES 
  jars = c("hec", "heclib", "rma", "hecData") 
  jars = paste0(dir.dssVue, "jar\\", jars, ".jar")
  ##libs = "-Djava.library.path=C:\\Program Files (x86)\\HEC\\HEC-DSSVue\\lib\\"
  libs = paste0("-Djava.library.path=",dir.dssVue,"lib\\")
  .jinit(classpath=jars, parameters=libs)
  
  # OUTPUT .DSS FILE 
  ## exract model name
  model.name <- basename(dir.model)
  ## List .DSS files in the model folder
  dss.fnames = list.files(dir.model,"*.dss",full.names=TRUE) 
  ## order .DSS files from most recent to oldest 
  dss.fnames.recent = file.info(dss.fnames)
  dss.fnames.recent.sort = dss.fnames.recent[rev(order(dss.fnames.recent$mtime)),]
  ## select output .DSS file (most recent .DSS file)
  out.dss.fname = paste0(dir.model,"\\",basename(row.names(dss.fnames.recent.sort[1,])))
  out.dss.file = opendss(out.dss.fname) 
  paths.all = dssrip::getCatalogedPathnames(out.dss.file)
  idx.flow = grep("/FLOW/",paths.all)
  paths.flow.all = paths.all[idx.flow]
  if(!is.null(flow.element)) {
    idx.element = grep(paste0("/",flow.element,"/"), paths.flow.all,ignore.case = TRUE)
  } else if(length(idx.element) < 1) {
      idx.element <- 1
  } else {
    idx.element <- 1
  }
  paths.flow <- paths.flow.all[idx.element]
  tsc.full = getFullTSC(out.dss.file,paths.flow[1])
  #### OC rainfall are in Pacific Time Zone but do not account for Daylight Savings time shifts
  #### Need to load time series in GMT timezone and then convert to "America/Los_Angeles" 
  #### Then need to add 8-hours (7-hours if daylight savings time) to account for the time difference between
  #### GMT and America/Los Angeles time zones
  datetime.GMT <- as.POSIXct(index(tsc.full),format = "%Y-%m-%d %H:%M",tz="GMT")
  datetime.PT <- format(datetime.GMT,tz = "America/Los_Angeles", usetz = TRUE)
  datetime.PT <- as.POSIXct(datetime.PT,format = "%Y-%m-%d %H:%M",tz = "America/Los_Angeles")
  date.index.char <- strftime(datetime.PT,tz = "America/Los_Angeles", usetz = TRUE)
  date.index.split <- lapply(date.index.char,strsplit," ")
  date.index.tz <- lapply(date.index.split, "[[",1)
  date.index.tz <- lapply(date.index.tz,"[[",3)
  date.index.PDT <- which(date.index.tz == "PDT")
  date.index.PST <- which(date.index.tz == "PST")
  datetime.PT[date.index.PDT] <- datetime.PT[date.index.PDT] + (60*60*7) 
  datetime.PT[date.index.PST] <- datetime.PT[date.index.PST] + (60*60*8) 
  tsc.full <- xts(coredata(tsc.full),datetime.PT)
  return(tsc.full)
  
}

### FUNCTIONS THAT CALCULATE PERFORMANCE METRICS ###
# function for NSE
calculate.NSE <- function(row.idx,name.streamGauge,flow.observed,flow.modeled) {
  flow.observed <- flow.observed
  if(is.na(nrow(flow.modeled))|is.na(ncol(flow.modeled))) {
    NSE <- 1- (sum((flow.modeled-flow.observed)^2,na.rm = TRUE)/sum((flow.observed-mean(flow.observed, na.rm = TRUE))^2,na.rm = TRUE))
  } else {
    NSE <- 1- (colSums((flow.modeled-flow.observed)^2,na.rm = TRUE)/sum((flow.observed-mean(flow.observed,na.rm = TRUE))^2,na.rm = TRUE))
  }
  resultsDF.NSE <- data.frame(streamGauge = name.streamGauge,NSE = NSE)
  rownames(resultsDF.NSE) <- row.idx
  return(resultsDF.NSE)
}

#function for PBias
calculate.pBias <- function(row.idx,name.streamGauge,flow.observed,flow.modeled) {
  flow.observed <- flow.observed
  if(is.na(nrow(flow.modeled))|is.na(ncol(flow.modeled))) {
    PBias <- (sum(flow.modeled - flow.observed,na.rm = TRUE)/sum(flow.observed,na.rm = TRUE))*100 
  } else {
  PBias <- (colSums(flow.modeled - flow.observed,na.rm = TRUE)/sum(flow.observed,na.rm = TRUE))*100
  }
  resultsDF.pBias <- data.frame(streamGauge = name.streamGauge,PBias = PBias)
  rownames(resultsDF.pBias) <- row.idx
  return(resultsDF.pBias)
}

#function to calculate % error for Peak Flow
calculate.peakFlowPercError <- function(row.idx,name.streamGauge,flow.observed,flow.modeled) {
  flow.observed <- flow.observed
  peakFlow.observed <- max(flow.observed,na.rm = TRUE)
  peakFlow.modeled <- sapply(flow.modeled,max,na.rm = TRUE)
  peakFlowPercError <- 100*(peakFlow.modeled - peakFlow.observed)/peakFlow.observed
  resultsDF.PeakFlow <- data.frame(streamGauge = name.streamGauge,peakFlow.observed = peakFlow.observed, peakFlow.modeled = peakFlow.modeled,percentError.peakFlow = peakFlowPercError)
  rownames(resultsDF.PeakFlow) <- row.idx
  return(resultsDF.PeakFlow)
}

# function to calculate % error for Total Flow Volume (from time series in cfs)
# timestep: the number of time units
# time.unit: the unit of time (i.e. "minutes","hours","days")
calculate.FlowVolumePercError <- function(row.idx,name.streamGauge,flow.observed,flow.modeled,timestep,time.unit) {
  flow.observed <- flow.observed
  if(time.unit == "minutes") {
    timestep.sec <- timestep*60
  } else if(time.unit == "hours") {
    timestep.sec <- timestep * 3600
  } else if(time.unit == "days") {
    timestep.sec <- timestep*3600*24
  } else {
    print("Time Error: time unit is not 'minutes,'hours',or 'days'")
  }
  flowVolume.observed <- sum(flow.observed*timestep.sec,na.rm = TRUE)
  flowVolume.modeled <- sapply(flow.modeled,sum,na.rm = TRUE)*timestep.sec
  FlowVolumePercError <- 100*(flowVolume.modeled - flowVolume.observed)/flowVolume.observed
  resultsDF.flowVolume <- data.frame(streamGauge = name.streamGauge,flowVolume.observed = flowVolume.observed,flowVolume.modeled = flowVolume.modeled,percentError.FlowVolume = FlowVolumePercError)
  rownames(resultsDF.flowVolume) <- row.idx
  return(resultsDF.flowVolume)
}

### FUNCTION TO CALCULATE ALL METRICS AND COMBINE INTO A DATAFRAME ###
calculatePerformanceMetrics <- function(row.idx,name.streamGauge,flow.observed,flow.modeled,timestep,time.unit) {
  resultsDF.NSE <- calculate.NSE(row.idx,name.streamGauge,flow.observed,flow.modeled)
  resultsDF.pBias <- calculate.pBias(row.idx,name.streamGauge,flow.observed,flow.modeled)
  resultsDF.peakFlow <- calculate.peakFlowPercError(row.idx,name.streamGauge,flow.observed,flow.modeled)
  resultsDF.flowVolume <- calculate.FlowVolumePercError(row.idx,name.streamGauge,flow.observed,flow.modeled,timestep,time.unit)
  resultsDF.all_1 <- merge(resultsDF.NSE,resultsDF.pBias, by = "row.names")
  resultsDF.all_2 <- merge(resultsDF.peakFlow,resultsDF.flowVolume, by = "row.names")
  resultsDF.all <- merge(resultsDF.all_1,resultsDF.all_2,by = "Row.names")
  metric.names <- c("NSE","PBias","peakFlow.observed","peakFlow.modeled","percentError.peakFlow","flowVolume.observed","flowVolume.modeled","percentError.FlowVolume")
  resultsDF.all <- resultsDF.all[,c("Row.names","streamGauge.x.x",metric.names)]
  colnames(resultsDF.all)[2] <- "streamGauge"
  
  return(resultsDF.all)
}



### FUNCTION TO REPLACE PARAMETER VALUES IN .BASIN FILE
#### basinFile: The HEC-HMS .basin file that contains all subbasin information
#### basinFile.path: the full path (including file name and extention) to the HEC-HMS .basin file that contains all subbasin information
#### element.type: a string indicating the type of element for which you are replacing parameters (i.e. "Subbasin","Reach", etc.)
#### element.name: that name of the element for which you are replacing parameter values
#### parameter.name: the name of the parameter for which you are replacing the value
#### parameter.value: the value of the parameter to which the current value will be changed
replaceParameterValue <- function(basinFile.path,basinFile.Copy, element.type,element.name,parameter.name,parameter.value,sensitivityAnalysis = FALSE) {
  basinFile <- readLines(basinFile.Copy)
  #### Return the values for rows in Elements.idx from the .basin file
  list.elements <- grep("^[^ \tEnd]",basinFile,value = TRUE)
  #### Identify and return row indices for lines that do not start with a white space and do not contain "End"
  idx.elements.start <- grep("^[^ \tEnd]",basinFile)
  #### Identiy and return row indices for lines in .basin file that contain "End:"
  idx.elements.end <- grep("End:",basinFile)
  idx.modelElements <- list(modelElement = list.elements, element.startLine = idx.elements.start, element.endLine = idx.elements.end)
  element.index <- grep(sprintf("%s: %s",element.type,element.name),idx.modelElements$modelElement)
  element.info <- basinFile[idx.modelElements$element.startLine[element.index]:idx.modelElements$element.endLine[element.index]]
  parameter.LineIndex <- grep(paste0(gsub("_"," ",parameter.name),":"),element.info)
  element.info[parameter.LineIndex] <- sprintf("     %s: %s",gsub("_"," ",parameter.name),parameter.value) ### there is a 5-space indent before element info
  basinFile[idx.modelElements$element.startLine[element.index]:idx.modelElements$element.endLine[element.index]] <- element.info
  if(sensitivityAnalysis == TRUE) {
    writeLines(basinFile,basinFile.path)
  } else {
    writeLines(basinFile,basinFile.Copy)
    writeLines(basinFile,basinFile.path) 
  }
  
}



### FUNCTION TO RANDOMIZE PARAMETER VALUES
randomizeParameterValue <- function(parameter,parameter.value.initial,range.factor,delta,lower.limit,upper.limit,n.samples) {
  if((parameter.value.initial-(delta*range.factor)) <= lower.limit) {
    bound.lower <- lower.limit
  } else {
    bound.lower <- parameter.value.initial-(delta*range.factor)
  }
  if((parameter.value.initial+(delta*range.factor)) >= upper.limit) {
    bound.upper <- upper.limit
  } else {
    bound.upper <- parameter.value.initial+(delta*range.factor)
  }
  parameterSample.randomized <- sample(seq(bound.lower,bound.upper,by=delta),n.samples,replace = T)
  return(parameterSample.randomized)
}



##########################################################################
##########################################################################

# SPECIFY ALL DIRECTORIES THAT WILL BE USED 

#### Directory where HEC-HMS.exe is located
dir.hechms = "C:\\Program Files\\HEC\\HEC-HMS\\4.5\\" 
#### Directory where HEC-DSSVue is located
dir.dssVue <- "C:\\Programs\\HEC-DSSVue-v3.0.00.212\\" 
#### Directory where scripts to run HEC-HMS are located
dir.runScripts <- "C:\\Users\\kflint\\Documents\\SanDiegoCreek\\CalibrationFiles\\"
#### Directory that contains the model folder(s)
dir.modelFolders <- "C:\\Users\\kflint\\Documents\\SanDiegoCreek\\HECHMS\\"
#### Directorys that contains parameter related files 
dir.parameter <- "C:\\Users\\kflint\\Documents\\SanDiegoCreek\\ModelParameters\\Parameters\\20200808\\forCalibration\\"
dir.parameterSpecs <- "C:\\Users\\kflint\\Documents\\SanDiegoCreek\\ModelParameters\\Parameters\\20200808\\forCalibration\\"
#### Directory that contains the observed flow files
dir.observedFlow <- "C:\\Users\\kflint\\Documents\\SanDiegoCreek\\CalibrationFiles\\ObservedFlow\\"
#### Directory that output results will be written to
dir.results <- "C:\\Users\\kflint\\Documents\\SanDiegoCreek\\CalibrationFiles\\Results\\"
#Directories with backup files
dir.backupHMSFiles <-"C:\\Users\\kflint\\Documents\\SanDiegoCreek\\BackupHMSFiles\\"
#### File that contains which storms should be used for calibration for each gauge
run.idx.file <- "C:\\Users\\kflint\\Documents\\SanDiegoCreek\\CalibrationFiles\\Runs_byGage.txt"


# Set directory containing model folders as working directory, and list the model folders within the directory
setwd(dir.modelFolders)
dir.models <- gsub("/","\\",list.dirs(path=dir.modelFolders,full.names = TRUE, recursive = FALSE)) #recursive = FALSE: do not list subdirectories within the subdirectory
model <- dir.models[1]

##load observed flow for all gauges
observedflow.name <- "CFS_hourly_SanDiegoCreek.txt"
flow.observed.all <- read.table(paste0(dir.observedFlow,observedflow.name),sep=",",header = TRUE, check.names = FALSE)
flow.observed.all <- xts(flow.observed.all[ ,2:ncol(flow.observed.all)], order.by = as.POSIXct(strptime(flow.observed.all$Datetime,"%Y-%m-%d %H:%M",tz = "America/Los_Angeles")))

# SPECIFY THE BASIN NAME
basin.name <- "SanDiegoCreek"
# SPECIFY RUNSCRIPT NAME
name.runScript <- paste0(basin.name,".compute.script")
# SPECIFY .BASIN FILE AND .MET FILE LOCATIONS, AND THEIR BACKUP FILE LOCAITONS (HMS ERASES .BASIN AND .MET FILES IF RUN FROM THE COMMAND LINE TOO MANY TIMES)
basinFile.name <- paste0(model,"\\",list.files(model,"*.basin"))
basinFile.Copy <- paste0(dir.backupHMSFiles,basin.name,".basin")
metFile.path <-paste0(model,"\\","Met_1.met")
metFile.Copy <- paste0(dir.backupHMSFiles,"Met_1.met")

# LOAD THE PARAMETER FILE AS A TABLE
#### First Column should indicate the gage that the subbasin immediately drains to. 
#### Second column shoud indicate the name of the subbasin.
#### Remaining columns should indicate the parameters for that basin.
#### The table should be ordered from upstream to downstream.
parameter.table.basin <- read.table(paste0(dir.parameter,"\\","parameters_subbasin_",basin.name,".txt"),header = TRUE,sep = ",")
parameter.table.reach <- read.table(paste0(dir.parameter,"\\","parameters_reach_",basin.name,".txt"),header = TRUE,sep = ",")

# Define elements and parameters
basins <- parameter.table.basin[,3]
reaches <- parameter.table.reach[,3]
elements <- list(Subbasin = basins, Reach = reaches)
row.names(parameter.table.basin) <- basins
row.names(parameter.table.reach) <- reaches
basin.parameters <- colnames(parameter.table.basin)[4:ncol(parameter.table.basin)]
reach.parameters <- colnames(parameter.table.reach)[4:ncol(parameter.table.reach)]
parameters <- list(Subbasin = basin.parameters, Reach = reach.parameters)

# Combine parameter dataframes into a list
parameter.tables <- list(Subbasin = parameter.table.basin, Reach = parameter.table.reach)

#### The table containing parameter specifications (Qualitative Uncertainty Ranking, Qualitative Uncertainty Factor,Min Value, Max Value)
parameter.specs.basin <-read.csv(paste0(dir.parameterSpecs,"parameterSpecifications_subbasin_",basin.name,".txt"),header = TRUE, sep = ",", row.names = 1)
parameter.specs.reach <-read.csv(paste0(dir.parameterSpecs,"parameterSpecifications_reach_",basin.name,".txt"),header = TRUE, sep = ",", row.names = 1)
parameter.specs <- list(Subbasin = parameter.specs.basin, Reach = parameter.specs.reach)

lower.limit.basin <- list()
upper.limit.basin <- list()
lower.limit.reach <- list()
upper.limit.reach <- list()

for(parameter in colnames(parameter.specs.basin)) {
  lower.limit.basin[[parameter]] <- parameter.specs.basin["Min_Value",parameter]
  upper.limit.basin[[parameter]] <- parameter.specs.basin["Max_Value",parameter]
}

for(parameter in colnames(parameter.specs.reach)) {
  lower.limit.reach[[parameter]] <- parameter.specs.reach["Min_Value",parameter]
  upper.limit.reach[[parameter]] <- parameter.specs.reach["Max_Value",parameter]
}

lower.limit <- list(Subbasin = lower.limit.basin, Reach = lower.limit.reach)
upper.limit <- list(Subbasin = upper.limit.basin, Reach = upper.limit.reach)

# DEFINE THE NAMES OF STREAM GAUGES USED FOR CALIBRATION
streamGauge.list <- unique(parameter.table.basin[,1])
junction.list <- unique(parameter.table.basin[,2])
timestep.model <- 15
time.unit.model <- "minutes"
timestep.observedFlow <- 1
time.unit.observedFlow <- "hours"

#DEFINE VARIABLES FOR SENSITIVITY ANALYSIS: 
#percentages = percentage deviation from the initial parameter value; 
#threshold = the minimum percent change in performance metrics to be considered a significant change
percentages <- c(-100,-50,-25,-10,-5,5,10,25,50,100)
threshold <- 5

# DEFINE THE NUMBER OF SAMPLES FROM PARAMETER TO BE SAMPLE FOR EACH PARAMETER
n.samples <- 200




############## SENSITIVTY ANALYSIS #########################################################  
###############################################################################################################################################
#Perform sensitivity analysis and create a dataframe of sensitive parameter intervals (as a percentage of the initital parameter value)

parameterIntervals.subset <- list()

for(element.class in names(elements)) {
#temp_element.class <- "Reach"
#for(element.class in temp_element.class) {
  n.element <- 1
  elements.subset <- elements[[element.class]][sample(1:length(elements[[element.class]]),ceiling(length(elements[[element.class]])*0.1),replace = FALSE)]
  parameterIntervals.subset[[element.class]] <- data.frame(matrix(ncol=length(parameters[[element.class]]),nrow=length(elements.subset)))
  colnames(parameterIntervals.subset[[element.class]]) <- parameters[[element.class]]
  rownames(parameterIntervals.subset[[element.class]]) <- elements.subset
  for(element in elements.subset) {
    n.parameter <- 1
    for(parameter in parameters[[element.class]]) {
      flow.modeled <- data.frame(double())
      parameter.value.initial <- parameter.tables[[element.class]][element,parameter]
      replaceParameterValue(basinFile.name,basinFile.Copy,element.class,element,gsub("_"," ",parameter),parameter.value.initial,sensitivityAnalysis = TRUE)
      metFile <- readLines(metFile.Copy)
      writeLines(metFile,metFile.path)
      run.HECHMS(dir.runScripts,name.runScript,dir.hechms)
      tsc.full <- dssToDataframe(dir.dssVue,model,flow.element=element)
      flow.modeled[1:nrow(tsc.full),1] <- tsc.full
      colnames(flow.modeled)[1] <- "Initial_Value"
      n.run <- 2
      for(percentage in percentages) {
        parameter.value <- parameter.value.initial * (1 + (percentage/100))
        if(parameter.value > upper.limit[[element.class]][[parameter]]|parameter.value < lower.limit[[element.class]][[parameter]]) {
          next
        }
        replaceParameterValue(basinFile.name,basinFile.Copy,element.class,element,gsub("_"," ",parameter),parameter.value,sensitivityAnalysis = TRUE)
        writeLines(metFile,metFile.path)
        run.HECHMS(dir.runScripts,name.runScript,dir.hechms)
        tsc.full <- dssToDataframe(dir.dssVue,model,flow.element=element)
        flow.modeled[1:nrow(flow.modeled),n.run] <- tsc.full
        colnames(flow.modeled)[n.run] <- percentage
        n.run <- n.run + 1
      }
      n.runs <- c(1:n.run-1)
      sensitivityMatrix <- calculatePerformanceMetrics(colnames(flow.modeled)[2:ncol(flow.modeled)],element,flow.modeled$Initial_Value,flow.modeled[,2:ncol(flow.modeled)],timestep.model,time.unit.model)
      sensitivityMatrix <- sensitivityMatrix[,c("Row.names","streamGauge","NSE","PBias","percentError.peakFlow","percentError.FlowVolume")]
      row.names(sensitivityMatrix) <- sensitivityMatrix$Row.names
      names.metrics <- names(sensitivityMatrix)[2:ncol(sensitivityMatrix)]
      sensitivityMatrix$NSE <- (1-(sensitivityMatrix$NSE))*100
      greaterThan.sensitivityThreshold.logical <- abs(sensitivityMatrix[,3:ncol(sensitivityMatrix)]) >= threshold
      rows.greaterThan.sensitivityThreshold <- which(greaterThan.sensitivityThreshold.logical==TRUE,arr.ind = TRUE)[,"row"]
      if(length(rows.greaterThan.sensitivityThreshold)== 0){
        minimumSensitiveDelta <- max(abs(as.numeric(row.names(sensitivityMatrix))))
      } else {
        minimumSensitiveDelta <- min(abs(as.numeric(row.names(sensitivityMatrix)[rows.greaterThan.sensitivityThreshold])))
      }
      parameterIntervals.subset[[element.class]][n.element,n.parameter] <- minimumSensitiveDelta
      n.parameter <- n.parameter + 1
    }
    n.element <- n.element + 1
  }
  write.table(parameterIntervals.subset[[element.class]],paste0(dir.results,"parameterIntervalssubset_",element.class,"_.txt"),sep = ",",row.names = TRUE, col.names = TRUE)
}
################# USE IF NEED TO RELOAD THE PARAMATERINTERVALS.SUBSET FILES ################################
sens.basin <- read.table("C:\\Users\\kflint\\Documents\\SanDiegoCreek\\CalibrationFiles\\Results\\parameterIntervalssubset_Subbasin_.txt",sep = ",",header = TRUE,check.names = FALSE)
sens.reach <- read.table("C:\\Users\\kflint\\Documents\\SanDiegoCreek\\CalibrationFiles\\Results\\parameterIntervalssubset_Reach_.txt",sep = ",",header = TRUE,check.names = FALSE)
parameterIntervals.subset[["Subbasin"]] <- sens.basin
parameterIntervals.subset[["Reach"]] <- sens.reach
################################################################################################################

parameterIntervals <- list(Subbasin = data.frame(matrix(ncol=length(basin.parameters),nrow=1)),Reach = data.frame(matrix(ncol=length(reach.parameters),nrow=1)))
colnames(parameterIntervals[["Subbasin"]]) <- basin.parameters
colnames(parameterIntervals[["Reach"]]) <- reach.parameters
for(element.class in names(elements)) {
  for(parameter in parameters[[element.class]]) {
    n.element <- 1
    parameterIntervals[[element.class]][parameter] <- quantile(parameterIntervals.subset[[element.class]][,parameter],0.1)
  }
}
  
#############################################################################################################################################################
#Create a sample of parameters for each basin and each parameters
parameter.samples <- list()
for(element.class in names(elements)) {
  for(element in elements[[element.class]]) {
    parameter.sample <- list()
    for(parameter in parameters[[element.class]]) {
      parameter.range.factor <- parameter.specs[[element.class]]["Qualitative_Uncertainty_Factor",parameter]
      parameter.sample[[parameter]] <- randomizeParameterValue(parameter,parameter.tables[[element.class]][element,parameter],parameter.range.factor,parameterIntervals[[element.class]][,parameter]*parameter.tables[[element.class]][element,parameter]/100,lower.limit[[element.class]][[parameter]],upper.limit[[element.class]][[parameter]],n.samples)
    }
    parameter.samples[[element.class]][[element]] <- parameter.sample
  }
}

save(parameter.samples,file = paste0(dir.results,"parameter_samples.Rdata"))


#####CALIBRATE SEMI-MANUALLY#############################
#"Initial_Canopy_Storage_Percent","Canopy_Storage_Capacity",
#"Crop_Coefficient","Initial_Surface_Storage_Percent","Surface_Storage_Capacity","Percent_Impervious_Area"
#"Initial_Deficit","Maximum_Deficit","Percolation_Rate","Snyder_Tp","Snyder_Cp","Recession_Factor"

calibration.parameters <- c("Canopy_Storage_Capacity","Surface_Storage_Capacity",
                            "Percent_Impervious_Area","Initial_Deficit","Maximum_Deficit","Percolation_Rate","Snyder_Tp",
                            "Snyder_Cp","Recession_Factor")
parameter.percent.delta <- list(Surface_Storage_Capacity = c(100),
                                Canopy_Storage_Capacity = c(100),
                                Percent_Impervious_Area = c(100),
                                Maximum_Deficit = c(100),
                                Initial_Deficit = c(100),
                                Percolation_Rate = c(20),
                                Snyder_Tp = c(100),
                                Snyder_Cp = c(57.1),
                                Recession_Factor = c(990)
                              )
                              #57.1
                                
# Randomize parameters for basins that drain to a given stream gauge 
run.idx.table <- read.table(run.idx.file,sep = ",",header = TRUE,check.names = FALSE)
n.gauge <- 4
gauge <- streamGauge.list[n.gauge]
run.idx.gauge <- run.idx.table[,gauge]
junction.gauge <- junction.list[grep(gauge,streamGauge.list)]
# identify the basins that drain to the current flow gauge
basins.currentGauge <- parameter.table.basin[which(parameter.table.basin[,1]==gauge),3]
reaches.currentGauge <- parameter.table.reach[which(parameter.table.reach[,1]==gauge),3]
elements.currentGauge <- list(Subbasin = basins.currentGauge, Reach = reaches.currentGauge)

for(run.idx in run.idx.gauge){
  runscript.file <- paste0(dir.runScripts,name.runScript)
  runscript.lines <- readLines(runscript.file)
  runscript.updatedLine <- sub(strsplit(runscript.lines[3],"\"")[[1]][2],paste0("Run ",run.idx),runscript.lines[3])
  runscript.lines[3] <- runscript.updatedLine
  write(runscript.lines,runscript.file)
  
  #### create an empty dataframe to write modeled flow to
  flow.modeled <- data.frame()
  element.class <- "Subbasin"
  # begin calibration
  n.storm <- 1
  for(n.run in c(1)){
    for(parameter in calibration.parameters){
      n.delta <- parameter.percent.delta[[parameter]][n.run] 
      #### assign parameter value from randomized parameter set to parameter
      for(element in elements.currentGauge[[element.class]]) {
        
        parameterValue.initial <- parameter.tables[[element.class]][element,parameter]
        parameterValue.n <- parameterValue.initial*(n.delta/100)
        if(parameterValue.n > upper.limit[[element.class]][[parameter]]) {
          parameterValue.n <- upper.limit[[element.class]][[parameter]]
        }else if(parameterValue.n < lower.limit[[element.class]][[parameter]]){
          parameterValue.n <- lower.limit[[element.class]][[parameter]]
        }
        #### replace existing parameter with new parameter
        replaceParameterValue(basinFile.name,basinFile.Copy,element.class,element,parameter,parameterValue.n,sensitivityAnalysis = FALSE)
      }
    }
  
  
    # RUN HEC-HMS WITHOUT GUI (IN COMMAND WINDOW) 
    #### define filename of script that runs HEC-HMS
    metFile <- readLines(metFile.Copy)
    writeLines(metFile,metFile.path)
    run.HECHMS(dir.runScripts,name.runScript,dir.hechms)
    
    # EXPORT DSS-FILES 
    #### Write most recent output flow .dss file to a dataframe 
    full.tsc <- dssToDataframe(dir.dssVue,model,flow.element=junction.gauge)
    #write column into output dataframe
    flow.modeled[1:nrow(full.tsc),n.storm] <- full.tsc
    run.current <- sprintf("delta_%s",n.delta)
    colnames(flow.modeled)[n.storm] <- run.current
    n.storm <- n.storm + 1
  }
  
  flow.modeled_xts <- as.xts(flow.modeled[,],order.by = as.POSIXct(index(full.tsc),format = "%Y-%m-%d %H:%M",tz = "America/Los_Angeles"))
  colnames(flow.modeled_xts) <- colnames(flow.modeled)
  #identify the end of each hour (for a 15-minute time series, the end of the hour period is the 45-mintue mark)
  ep <- endpoints(index(flow.modeled_xts),on="hours")
  # take the average of the flow for the hour period (00:00, 00:15, 00:30, 00:45)
  flow.modeled_xts_hourly <- period.apply(flow.modeled_xts,ep,mean)
  colnames(flow.modeled_xts_hourly) <- colnames(flow.modeled_xts)
  #shift the time index back by 45 minutes (in unites of seconds), to the 00:00 mark for that hour
  index(flow.modeled_xts_hourly) <- index(flow.modeled_xts_hourly) - (15*60*3)
  last.time.value <- format(index(flow.modeled_xts_hourly)[nrow(flow.modeled_xts_hourly)],"%H:%M")
  if(last.time.value != "00:00") {
    flow.modeled_xts_hourly <- flow.modeled_xts_hourly[1:nrow(flow.modeled_xts_hourly)-1,]
  }
  #write output time series to csv
  dir.junction <-  paste0(dir.results,gauge,"_",junction.gauge,"\\")
  dir.create(dir.junction)
  write.zoo(flow.modeled_xts_hourly, file=paste0(dir.results,junction.gauge,"_storm",run.idx,".txt"), sep = ",")
}

  
######## START HERE 20201004 ################################################ 
#### HERE WE CALCULATE METRICS FOR RESULTS AND CHOOSE PARAMETER VALUE THAT YIELDED THE BEST RESULTS AND INPUT INTO .BASIN FILE
# load observed flow file as a table

################################################################################################################################
############### USE IF RUNNING THIS SECTION INDEPENDTLY OF THE PREVIOUS SECTION OF CODE ########################################
n.gauge <- 4
gauge <- streamGauge.list[n.gauge]
junction.gauge <- junction.list[n.gauge]

flow.modeled.folder <- paste0(dir.results,gauge,"\\")
setwd(flow.modeled.folder)
flow.modeled.files <- list.files(pattern = "*.txt")


#################################################################################################################################
################################################################################################################################
# DEFINE THE PERFORMANCE METRICS, THIER RELATIVE WEIGHTS, THE POINTS ASSIGNED TO EACH RANKED INTERVAL
performanceMetrics.weights <- list(percentError.FlowVolume = 0.35, percentError.peakFlow = 0.40, NSE = 0.25)
### Define the points assigned to each "rank" interval
ranks.percErrorFlowVolume <- c(5,10,25,50,75)
ranks.percErrorPeakFlow <- c(5,10,25,50,75)
#ranks.pBias <- c(5,10,25,50,75)
### NOTE NSE WILL BE MULTIPLIED BY (-1) SO THAT THE DESIRABILITY ORDER (SMALLEST = DESIRABLE VALUE; LARGEST = UNDESIRABLE VALUE) MATCHES THAT OF THE OTHER METRICS
ranks.NSE <- c(-0.90,-0.75,-0.5,-0.25,0)
assignedPoints <- c(6,5,4,3,2,1)

for(file in flow.modeled.files){
  storm <- strsplit(strsplit(file,"_")[[1]][3],"\\.")[[1]][1]
  flow.modeled.table <- read.table(file,sep = ",",header = TRUE,check.names = FALSE)
  date.index <- flow.modeled.table[,1]
  #date.index <- flow.modeled.table$Index
  date.index <- as.POSIXlt(date.index,format = "%Y-%m-%d %H:%M",tz="America/Los_Angeles")
  if(ncol(flow.modeled.table) > 2){
    flow.modeled_xts_hourly <- xts(flow.modeled.table[,2:ncol(flow.modeled.table)],order.by = date.index)
  }else{
    flow.modeled_xts_hourly <- xts(data.frame(flow.modeled.table[,2,drop=FALSE]),order.by = date.index)
  }
  ep.flow.modeled <- endpoints(flow.modeled_xts_hourly,"hours")
  flow.modeled_xts_hourly <- period.apply(flow.modeled_xts_hourly,ep.flow.modeled,mean,na.rm = TRUE)
  date.index <- index(flow.modeled_xts_hourly)
  
  flow.modeled.start <- which(index(flow.observed.all) == index(flow.modeled_xts_hourly)[1])
  flow.modeled.end <- which(index(flow.observed.all) == index(flow.modeled_xts_hourly[length(index(flow.modeled_xts_hourly))]))
  flow.observed <- flow.observed.all[flow.modeled.start:flow.modeled.end,gauge]
  ### confirm that both flow.observed and flow.modeled_xts_hourly are "double" type arrays
  if(typeof(flow.observed) != "double") {
    flow.observed <- type.convert(flow.observed,"double")
  }
  if(typeof(flow.modeled_xts_hourly) != "double") {
    flow.modeled_xts_hourly <- type.convert(flow.modeled_xts_hourly,"double")
  }
  # Caclculated performante metrics
  performanceMetrics <- calculatePerformanceMetrics(colnames(flow.modeled_xts_hourly),gauge,as.vector(coredata(flow.observed)),flow.modeled_xts_hourly,timestep.observedFlow,time.unit.observedFlow)
  #performanceMetrics <- calculatePerformanceMetrics(colnames(flow.modeled_xts),gauge,flow.observed[,1:ncol(flow.observed)],flow.modeled_xts_hourly[,1:ncol(flow.modeled_xts_hourly)],timestep.observedFlow,time.unit.observedFlow)
  row.names(performanceMetrics) <- performanceMetrics$Row.names
  metric.names <- c("NSE","percentError.peakFlow","percentError.FlowVolume")
  outfile.PerformanceMetrics <-paste0(dir.results,junction.gauge,"_",storm,"_PerformanceMetrics.txt")
  write.table(performanceMetrics,outfile.PerformanceMetrics,sep = ",", row.names = FALSE, col.names = TRUE)

  performanceMetrics.forCalc <- performanceMetrics[,metric.names]
  ### multiply NSE by -1 so that the directionality of more desirable values is in the same direction as the other methrics
  performanceMetrics.forCalc$NSE <- -1*performanceMetrics.forCalc$NSE
  ### Set percent error and bias metrics to absolute values
  performanceMetrics.forCalc$percentError.FlowVolume <- abs(performanceMetrics.forCalc$percentError.FlowVolume)
  performanceMetrics.forCalc$percentError.peakFlow <- abs(performanceMetrics.forCalc$percentError.peakFlow)
  #performanceMetrics.forCalc$PBias <- abs(performanceMetrics.forCalc$PBias)
  ### Define weights of performance metrics
  
  performanceGrades <- data.frame(percentError.FlowVolume = ranks.percErrorFlowVolume,percentError.peakFlow = ranks.percErrorPeakFlow,NSE = ranks.NSE)
  ### Assign points for each individual metric value (more points = more desirable metric value)
  for(metric in metric.names) {
    for(i in 1:length(assignedPoints)) {
      if(i==1) {
        performanceMetrics.forCalc[which(performanceMetrics.forCalc[[metric]] <= performanceGrades[i,metric]),metric] <- assignedPoints[i]
      } else if(i==length(assignedPoints)) {
        performanceMetrics.forCalc[which(performanceMetrics.forCalc[[metric]] > performanceGrades[i-1,metric]),metric] <- assignedPoints[i]
      } else {
        performanceMetrics.forCalc[which(performanceMetrics.forCalc[[metric]] <= performanceGrades[i,metric] & performanceMetrics.forCalc[[metric]] > performanceGrades[i-1,metric]),metric] <- assignedPoints[i]
      }
    }
  }
  
  ### Calculate the weighted average of points accumulated for all metrics
  performanceMetric.weighted <- performanceMetrics.forCalc
  for(metric in metric.names){
    performanceMetric.weighted[metric] <- performanceMetric.weighted[[metric]] * performanceMetrics.weights[[metric]]
  }
  performanceMetric.weighted <- data.frame(metric.Weighted = rowSums(performanceMetric.weighted))
  row.names(performanceMetric.weighted) <- row.names(performanceMetrics)
  performanceMetric.weighted <- performanceMetric.weighted[order(performanceMetric.weighted["metric.Weighted"], decreasing = TRUE), ,drop = FALSE]
  #### Identify the top performing 5% of the parameter sets
  #top.five.percent <- performanceMetric.weighted[1:ceiling(nrow(performanceMetric.weighted)*0.05),"metric.Weighted"]
  #top.five.percent_runs <- row.names(performanceMetric.weighted)[1:ceiling(nrow(performanceMetric.weighted)*0.05)]
  #performanceMetrics_top.five.percent <- performanceMetrics[top.five.percent_runs[1:length(top.five.percent_runs)],]
  #performanceMetrics_top.five.percent[["weighted_Metric"]] <- top.five.percent
  outfile.PerformanceMetric.weighted <- paste0(dir.results,junction.gauge,"_",storm,"_PerformanceMetrics_weighted.txt")
  write.table(performanceMetric.weighted,outfile.PerformanceMetric.weighted,sep = ",", row.names = TRUE, col.names = TRUE)
  
  # for(element.class in names(elements)) {
  #   rank <- 1
  #   for(run in top.five.percent_runs) {
  #     n.run <- as.numeric(sub("run_","",run))
  #     parameters_n.run <- data.frame()
  #     for(element in elements[[element.class]]){
  #       for(parameter in names(parameter.samples[[element.class]][[1]])) {
  #         parameters_n.run[element,parameter] <- parameter.samples[[element.class]][[element]][[parameter]][n.run] 
  #       }
  #     }
  #     outfilename.ParameterSets <- paste0(dir.results,junction.gauge,"_",element.class,"_",storm,"_ParameterSet_Rank",rank,"_",sub("run ","set",run),".txt")
  #     write.table(parameters_n.run,outfilename.ParameterSets,sep = ",",row.names = TRUE, col.names = TRUE)
  #     rank <- rank + 1
  #   }
  # }
  
  png(file=paste0(dir.results,gauge,"\\Plot\\",storm,".png"),width=8, height=6,unit="in", pointsize = 6, res = 400)
  ymax <- max(flow.observed,sapply(flow.modeled_xts_hourly,max,na.rm = TRUE),na.rm = TRUE)
  plot(date.index,flow.observed,type = "l",ylim = c(0,ymax))
  axis.POSIXct(1,date.index, format="%Y-%m-%d")
  line.colors <- c("black")
  line.names <- "Observed"
  n.variable <- 2
  for (i in 1:ncol(flow.modeled_xts_hourly)) {
    lines(date.index,flow.modeled_xts_hourly[,i],col=colors()[i*10])
    line.colors[n.variable] <- colors()[i*10]
    line.names[n.variable] <- colnames(flow.modeled_xts_hourly)[i]
    n.variable <- n.variable + 1
  }
  
  legend(min(date.index),ymax,legend=line.names,col=line.colors,lty = 1)  
  dev.off()
}

#### USE IF NEED TO CHANGE INITITAL SUBBASIN PARAMETERS (OTHER THAN MONTHLY LIMIT)

element.class <- "Subbasin"
parameter.table <- read.table("C:\\Users\\kflint\\Documents\\SanDiegoCreek\\ModelParameters\\Parameters\\20200808\\parameters_subbasin.txt",sep=",",header = TRUE)
elements <- parameter.table$HMS_ID
row.names(parameter.table) <- elements

n.gauge <- 2
gauge <- streamGauge.list[n.gauge]
junction.gauge <- junction.list[grep(gauge,streamGauge.list)]
# identify the basins that drain to the current flow gauge
basins.currentGauge <- parameter.table.basin[which(parameter.table.basin[,1]==gauge),3]
reaches.currentGauge <- parameter.table.reach[which(parameter.table.reach[,1]==gauge),3]
elements.currentGauge <- list(Subbasin = basins.currentGauge, Reach = reaches.currentGauge)

elements <- basins.currentGauge

#"Initial_Canopy_Storage_Percent","Canopy_Storage_Capacity",
#"Crop_Coefficient","Initial_Surface_Storage_Percent",
#"Surface_Storage_Capacity","Initial_Deficit","Maximum_Deficit",
#"Percolation_Rate","Snyder_Tp","Snyder_Cp","Recession_Factor"

parameters <- c("Initial_Deficit")
percent.delta <- -40
for(element in elements) {
  for(parameter in parameters) {
    parameter.value.initial <- parameter.table[element,parameter]
    parameter.value.n <- parameter.value.initial
    #parameter.value.n <- parameter.value.initial*(1 + (percent.delta/100))
    #parameter.value.n <- .99
    replaceParameterValue(basinFile.name,basinFile.Copy,element.class,element,gsub("_"," ",parameter),parameter.value.n,sensitivityAnalysis = FALSE)
  }  
}

element.class <- "Reach"
parameter.table <- read.table("C:\\Users\\kflint\\Documents\\SanDiegoCreek\\ModelParameters\\Parameters\\20200808\\parameters_reach.txt",sep=",",header = TRUE)
elements <- parameter.table$HMS_ID
row.names(parameter.table) <- elements
parameters <- c("Index_Flow")
for(element in elements) {
  for(parameter in parameters) {
    parameter.value.initial <- parameter.table[element,parameter]
    replaceParameterValue(basinFile.name,basinFile.Copy,element.class,element,gsub("_"," ",parameter),parameter.value.initial,sensitivityAnalysis = FALSE)
  }  
}


