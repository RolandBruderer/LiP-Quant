#       Copyright (c) 2020 Roland Bruderer, Biognosys AG, Switzerland
#   
#         Permission is hereby granted, free of charge, to any person obtaining a copy
#       of this software and associated documentation files (the "Software"), to deal
#       in the Software without restriction, including without limitation the rights
#       to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#       copies of the Software, and to permit persons to whom the Software is
#       furnished to do so, subject to the following conditions:
#   
#         The above copyright notice and this permission notice shall be included in all
#       copies or substantial portions of the Software.
# 
#       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#       IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#       FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#       AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#       LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#       OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#       SOFTWARE.
# 
#       Title:  Limited proteolysis dose-response based drug deconvolution data analysis algorithm based on machine learning scores
#         
#                                                                         by Roland Bruderer
#
#       Function: This algorithm performs dose-response fit correlation and testing.
#       The ranking of the peptide and protein candidates is achieve based on by multi score machine learning (ML) using LDA.
#       This was performed based on a training set. The ranking of LIP dose response experiments can be executed. 
#
#       Required input: 
#   
#       1.)Spectronaut Report: RnD Report Schema, required column headers column names: EG.ModifiedSequence,  R.Condition, R.FileName, R.Replicate, PG.Coverage, PG.ProteinAccessions, PG.ProteinNames,  
#       EG.PrecursorId,  EG.StrippedSequence, EG.Qvalue, FG.NormalizedMS2PeakArea and PEP.Quantity                  
#       2.)Spectronaut t-test report (unfiltered, remove default fold change and q-value filters), required column headers: (Comparison (group1/group2),	Group,	ProteinGroups,	AVG Log2 Ratio,	Absolute AVG Log2, Qvalue,	UniProtIds,	Genes,	ProteinDescriptions,	ProteinNames,	# Unique Total Peptides)
#       3.)Optionally, a know drug target list file (for training set required) with two columns: ProteinNames (Uniport name, e.g. FPKB1_HUMAN) and Target, for machine learning to generate a weightned scores based on mulitple experiments, 
#       add a third column to the known targets file named: ExpName and a column called priorEC50 for calculation of delta EC50 values
#       This column should contain the file names of the experiments (with file extension) e.g. RapaDoseResp.txt
#       4.)Optionally, a concentration sheet (will be generated based on conditions present in the experiment)
#       ML based on Qvalue and FC filtered list (background drawn from there)
#       5.) for human experiments a LIP-2Protein-Frequency library can be used (by default OFF), columnn headers: PFL.PG and PLFrequency
#       Requirements for training sets:        
#
#       Concentrations (normal, not log scale) either filled by defaults or externally defined via lookup table,  required column headers: R.Condition, Concentrations
#       Automated default condition to concentration annotation or with lookup file ( normal, not log concentrations)
#       Random set drawn on protein level
#       Normalized intensities are optionally log transformed and then averaged per run for modified peptides
#       A dose-response correlation (DRC) test is performed and the pvalues are corrected by qvalue. The pvalue from the correlation test will be used for the ML.
#       
#
#
# Intro End

#_____________________________________________________________________________________________
#
#Part 0: Initilalization
#_____________________________________________________________________________________________
## section0-package loading, clearing ----
#check for required packages
list.of.packages <- c("curl","kernlab","caret", "drc", "dplyr", "ggplot2", "data.table", "plyr", "plotrix", "RColorBrewer", 
                      "stringr", "tools", "reshape2", "zoo", "scales", "devtools", "emdbook", "Hmisc", "zoo", "yaml", "glue", 
                      "rlang", "digest", "stringi", "ModelMetrics", "robustbase", "ipred", "e1071", "rpart", "DataCombine")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load necessary libraries
library(data.table)
library(plyr)
library(ggplot2)
library(plotrix)
library(RColorBrewer)
library(stringr)
library(tools)
library(reshape2)
library(zoo)
library(scales)
library(dplyr)
library(drc)
library(caret)
library(devtools)
library(emdbook)
library(Hmisc)

#BiocManager::install("qvalue") 
#BiocManager::installbiocLite("svm") 

library(qvalue)
library(e1071)
library(rpart)
library(DataCombine)


#clear work space
rm(list = ls(all = TRUE))

cat( "   \n" )
cat( "   Script Title:  Limited proteolysis dose-response based drug deconvolution data analysis algorithm based on machine learning scores.\n
Written and tested in R V.3.5.0 (2018-04-23) / Rstudio V.1.1.423 by R.Bruderer\n
" )
cat( "   \n" )

#manual parameters to set: 

#Qvalue cutoff candidate list
Qcutoff <- 0.01
#foldchange cutoff (log2) candidate list
FCcutoff <- 0.46
#LDA background set
ldaBGsize <- 400
#optional: concentration to be excluded from analysis
#RemConc <- c(200000000)


#_____________________________________________________________________________________________

#setting of pre-calculated score weights based on Rapamycin and Calyculin A and Staurosporine
scoreweights.v <- c(1.317320, 9.191420, 1.332628, 1.452496)



#_____________________________________________________________________________________________
#Question if score weighting should be performed
answeropt <- c("YES", "NO", "_____________________________________________________")

# selection default
cat("Please state, if adjusted score weighting is to be applied based upon machine leaning with a training set (default: NO).\n" )
MLtraining <- select.list( answeropt, multiple=T, title="Please state, if adjusted score weighting is to be applied based upon machine leaning with a training set (default: NO)." )
if ( length(MLtraining) == 0 ) {
  MLtraining <- "NO"
  cat("No adjusted score weighting will be executed, weights from rapamycin (exp 18), calyculin A (exp 20) and staurosporine (exp 36) will be applied.\n" )
} else { cat("Your answer was  ", paste( MLtraining, collapse=", "), "\n", sep="" )
}
if (  MLtraining == "NO" ) { 
  cat("No score weightning will be executed, weights from rapamycin (exp 18), calyculin A (exp 20) and staurosporine (exp 36) will be applied.\n" )
}

answercoropt <- c("pearson", "kendall", "spearman", "_____________________________________________________")
# selection default #any cases when to use spearman or kendall  
cat("Please state, what correlation should be applied (default: pearson).\n" )
answercor <- select.list( answercoropt, multiple=T, title="Please state, what correlation should be applied (default: pearson)." )
if ( length(answercor) == 0 ) {
  answercor <- "pearson"
  cat("Pearson correlation will be executed.\n" )
} else { cat("The selected correlation method was ", paste( answercor, collapse=", "), "\n", sep="" )
}

answerlogopt <- c("Normal", "Log", "_____________________________________________________")
# selection default 
cat("Please select what intensity value transformation should be applied (default: normal).\n" )
answerlog <- select.list( answerlogopt, multiple=T, title="Please select what intensity value transformation should be applied (default: normal)." )
if ( length(answerlog) == 0 ) {
  answerlog <- "Normal"
  cat("Normal intensities will be used.\n" )
} else { cat("The transformation was method was ", paste( answerlog, collapse=", "), "\n", sep="" )
}

#_____________________________________________________________________________________________
#Question if dose response correlation testing should performed
# selection default
cat("Please state, if dose response correlation testing should performed (default: YES).\n" )
DRCtesting <- select.list( answeropt, multiple=T, title="Please state, if dose response correlation testing should performed (default: YES)." )
if ( length(DRCtesting) == 0 ) {
  DRCtesting <- "YES"
  cat("Dose response correlation will be performed.\n" )
} else { cat("Your answer was  ", paste( DRCtesting, collapse=", "), "\n", sep="" )
}
if (  DRCtesting == "NO" ) { 
  cat("No dose response correlation testing will performed.\n" )
}

#_____________________________________________________________________________________________
#load LIP-PFL

#preparation
answeroptPFL <- c("Automatic","YES", "NO", "_____________________________________________________")

cat("Please select, if you select a LIP-Protein-Frequency-library/Crapome.\n" )
PFLselection <- select.list( answeroptPFL, multiple=T, title="Please select, if you select a LIP-Protein-Frequency-library/Crapome..Default is auto)" )
if ( length(PFLselection) == 0 || (PFLselection) == "Automatic") {
  cat("Auto import from file .\n" )
  pfl.df  <- read.delim("/Users/.../LiP-Quant_Crapome.txt", header=T, sep="\t", stringsAsFactors=F, na.strings=c( "n. def." ) ) #add file and path
  setnames(pfl.df,"Frequency","PLFrequency")
  PFLselection <- "Automatic"
} else { cat("Your answer was  ", paste( PFLselection, collapse=", "), " .\n", sep="" )
}

if ( (PFLselection) == "YES" ) {
  #import and merge
  cat( "Select LIP-Protein-Frequency-library/Crapome  \n" )
  temp <- file.choose()
  temp1 <- sub("", "", basename(temp))
  for (i in 1:length(temp)) assign(temp1[i], read.delim(temp[i]))
  
  #list names of data frames to merge
  file_namesPFL=(temp1)
  #preprocessing of loaded data frames
  #Name the files in a new column
  names(file_namesPFL) <- file_namesPFL
  pfl.df = ldply(file_namesPFL, get)
  
}

#if ( (PFLselection) != "NO" ) {
#  #annotation of LIP PFL
#  setnames(pfl.df,"UniProtIds","PFL.PG")}

#generated an empty data frame, if no PFL selected
if ( (PFLselection) == "NO" ) {
  pfl.df <- data.frame(c(1),c(1))
  setnames(pfl.df,"c.1.","PFL.PG")
  setnames(pfl.df,"c.1..1","PLFrequency")
}

#_____________________________________________________________________________________________
#Perform score weighting based on training set(s)
if (  MLtraining == "YES" ) {
  #data loading
  #import and merge
  cat( "Select unfiltered candidate list for ML training set  \n" )
  temp <- file.choose()
  temp1 <- sub("", "", basename(temp))
  for (i in 1:length(temp)) assign(temp1[i], read.delim(temp[i]))
  
  #list names of data frames to merge
  file_names1=(temp1)
  #preprocessing of loaded data frames
  #Name the files in a new column
  names(file_names1) <- file_names1
  t.df = ldply(file_names1, get)
  
  
  #set working directory
  tempwd <- temp[1]
  workingDir <- sub("", "", dirname(tempwd))
  setwd(workingDir)
  
  #comparison filtering
  cat( "   Comparison filtering\n" )
  setnames(t.df,"Comparison..group1.group2.","Filter")
  
  #merge file-name and comparison for multiple candidate lists
  if(length(unique(t.df$Filter)) != 1) {
    t.df$Filter <-
      do.call(paste, c(t.df[c(".id", "Filter")], sep = ""))
  }

  
  #qvalue cutoff filtering
  t.df <- subset(t.df, t.df$Qvalue <= Qcutoff)
  #optional FC cutoff filtering
  t.df <- subset(t.df, t.df$Absolute.AVG.Log2.Ratio >= FCcutoff)
  
  #check if correlation pvalues not out off range for stable scores
  if(sum(-log(t.df$Qvalue,10) > 60) > 0) {
    cat("Qvalues out of range (-log(Qvalue)) > 60, redefine Q-value score range, stopping")
    stop()
  }
  
  #_____________________________________________________________________________________________
  # selection conditions column to retain
  cat( " Selection of the comparison to retain (e.g. 10 to 100x IC50).\n" )
  #get comparison column levels
  filter.v <- unique(t.df$Filter)
  filter.v <- as.character(filter.v)
  filter.v <- append(filter.v, "_____________") 
  
  cond.true.v <- select.list( filter.v, multiple=T, title="Selection of the comparison to retain (e.g. 10 to 100x IC50)" )
  if ( length( cond.true.v ) == 0 ) {
    cat( " No conditions selected Stopping.\n" )
    stop()
  }
  cat( " conditions selected: ", paste( cond.true.v, collapse=", " ), "\n", sep="" )
  cat( "   \n" )
  
  #Filtering based on condition
  #true
  t1.df <- t.df[(t.df$Filter %in% cond.true.v), ]
  t1.df$Qvalue <- as.numeric(t1.df$Qvalue)
  t1.df <- subset(t1.df,!is.na(t1.df$Qvalue))
  t1.df <- t1.df[with(t1.df, order(Qvalue)), ]
  

  #check if protein names exist, otherwise take accessions
  if(length(unique(t1.df$ProteinNames)) == 1) {
    t1.df$ProteinNames <- NULL
    setnames(t1.df,"ProteinGroups","ProteinNames")
  }
  
  #Protein annotation LIP_PFL
  t1.df$PFL.PG <- gsub("\\;.*","",t1.df$ProteinGroups)
  
  #comparison table protein name
  t1.df$Comp.PG <-
    do.call(paste, c(t1.df[c(".id", "ProteinNames")], sep = ""))
  
  #saving candidate list file name
  filename <- str_sub(file_names1, end=-5)
  
  
  #_____________________________________________________________________________________________
  #load lookup table with know drug candidates
  cat( "  Select lookup table with known drug targets for ML training set (two columns called: ProteinNames and Target) \n" )
  #data loading
  #import and merge
  temp <- file.choose()
  temp2 <- sub("", "", basename(temp))
  for (i in 1:length(temp)) assign(temp2[i], read.delim(temp[i]))
  
  #list names of data frames to merge
  file_names2=(temp2)
  #preprocessing of loaded data frames
  #Name the files in a new column
  names(file_names2) <- file_names2
  lookup.df = ldply(file_names2, get)
  
  lookup.df$.id <- NULL
  lookup.df$Target <- as.character(lookup.df$Target)
  
  #merge file-name and comparison for multiple candidate lists
  if("ExpNameComp" %in% colnames(lookup.df)) {
    #comparison table protein name
    lookup.df$Comp.PG <-
      do.call(paste, c(lookup.df[c("ExpNameComp", "ProteinNames")], sep = ""))
    
    lookup.df$SNRep.PG <-
     do.call(paste, c(lookup.df[c("ExpNameSN", "ProteinNames")], sep = ""))
  }
  lookup.df$ProteinNames <- NULL
  lookup.df$ExpNameComp <- NULL
  lookup.df$ExpNameSN <- NULL
  
  #________________________________________________________________________________________________
  #annotate protein frequency library in comparison table
  t1.df <- merge(t1.df, pfl.df, by=c("PFL.PG"), all.x=TRUE, all.y=FALSE)
  t1.df$PFL.PG <- NULL
  t1.df$PLFrequency[is.na(t1.df$PLFrequency)] <- 0
  
  #_____________________________________________________________________________________________
  #load SN report
  cat( "  Select Spectronaut quantitative Report (SN12 - long format) for ML training set \n" )
  
  temp <- file.choose()
  temp3 <- sub("", "", basename(temp))
  for (i in 1:length(temp)) assign(temp3[i], read.delim(temp[i]))
  
  #list names of data frames to merge
  file_names3=(temp3)
  #preprocessing of loaded data frames
  #Name the files in a new column
  names(file_names3) <- file_names3
  r.df = ldply(file_names3, get)
  
  #saving SN report file name
  filename1 <- str_sub(file_names3, end=-5)
  
  #check if external format loaded
  if ("EG.ModifiedPeptide" %in% colnames(r.df)) {
    setnames(r.df,"EG.ModifiedPeptide","EG.ModifiedSequence") }
  
  #check if protein names exist, otherwise take accessions
  if(length(unique(t1.df$PG.ProteinNames))== 1) {
    setnames(r.df,"PG.ProteinAccessions","ProteinNames")
  } else {
    setnames(r.df,"PG.ProteinNames","ProteinNames")
  }
  #merge experiment and protein accession for multiple candidate lists
  # if(length(unique(r.df$.id)) != 1) {
    #sn report protein
    r.df$SNRep.PG <-
      do.call(paste, c(r.df[c(".id", "ProteinNames")], sep = ""))
  # }
  

  #getting protein coverages and proteotypicity
  if ("PG.Coverage" %in% colnames(r.df) && "PEP.IsProteotypic" %in% colnames(r.df)) {
    rc.df <- subset(r.df, !duplicated(r.df$EG.ModifiedSequence))
    rc.df <- subset(rc.df, select= c(EG.ModifiedSequence, PG.Coverage, PEP.IsProteotypic))}
  
  if ("PG.Coverage" %in% colnames(r.df) && !("PEP.IsProteotypic" %in% colnames(r.df))) {
    rc.df <- subset(r.df, !duplicated(r.df$EG.ModifiedSequence))
    rc.df <- subset(rc.df, select= c(EG.ModifiedSequence, PG.Coverage))}
  
  if (!("PG.Coverage" %in% colnames(r.df)) && "PEP.IsProteotypic" %in% colnames(r.df)) {
    rc.df <- subset(r.df, !duplicated(r.df$EG.ModifiedSequence))
    rc.df <- subset(rc.df, select= c(EG.ModifiedSequence, PEP.IsProteotypic))}
  
  
  #candidate list processing (ranking, annotation)
  ## section2 ----
  #_____________________________________________________________________________________________
  #candidate list processing (ranking, annotation)
  cat("Annotating sequences of unfilterd candidate list (target/protein coverage).\n" )
  ####t1.df$Target <- NULL 
  #annotating drug targets and background experiment wise in comparison candidate list
  t1.df <- merge(t1.df, lookup.df, by=c("Comp.PG"), all.x=TRUE, all.y=FALSE)
  t1.df$Target[is.na(t1.df$Target)] <- "BG"
  
  #add protein coverage (and peptides per protein, Proteotypic?) (in case of multiple not solved)
  t1.df$EG.ModifiedSequence <- t1.df$Group
  t1.df <- merge(t1.df, rc.df, by=c("EG.ModifiedSequence"), all.x=TRUE, all.y=FALSE)
  
  #write.table(t1.df, "CompReport.txt", sep="\t", row.names = FALSE, quote = FALSE)
  
  # Preparing SN report and performing dose response analysis
  ## section3-DR anlysis training set----
  #_____________________________________________________________________________________________
  #Dose response analysis training set
  # Preparing SN report and performing dose response analysis
  cat("Preparing SN report and performing dose response analysis.\n" )
  #combine experiments, runs, conditions and proteins
  r.df$RunModseq <-
    do.call(paste, c(r.df[c(".id","R.FileName", "EG.ModifiedSequence")], sep = ""))
  r.df$CondModseq <-
    do.call(paste, c(r.df[c(".id","R.Condition", "EG.ModifiedSequence")], sep = ""))
  r.df$ProtModSeq <-
    do.call(paste, c(r.df[c("SNRep.PG", "EG.ModifiedSequence")], sep = ""))
  r.df$ExpCond <-
    do.call(paste, c(r.df[c(".id","R.Condition")], sep = ""))
  
  
  #annotate known candidates for training set true and false separation
  r2a.dt <- merge(r.df,lookup.df, by=c("SNRep.PG"), all.x= T)
  r2a.dt$Target[is.na(r2a.dt$Target)] <- "BG"
  r2aTARG.dt <- subset(r2a.dt, r2a.dt$Target == "YES")
  r2aBG.dt <- subset(r2a.dt, r2a.dt$Target == "BG")
  
  #get training set
  #get know targets significantly enriched
  ttrg.df <- subset(t1.df, t1.df$Target == "YES")
  ttrg.df <- subset(ttrg.df, ttrg.df$Qvalue < (Qcutoff/100))
  ttrg.df <- data.frame(ttrg.df$Group)
  setnames(ttrg.df, "ttrg.df.Group", "EG.ModifiedSequence")
  #get known candidates for DRC
  r2at.dt <- merge(r2aTARG.dt,ttrg.df, by="EG.ModifiedSequence", all.x=FALSE, all.y= T)
  
  #get random background set form candidate list for SN report extraction
  r2ab.dt <- subset(t1.df, t1.df$Target == "BG")
  
  setnames(r2ab.dt, "Group", "EG.ModifiedSequence")
  #sampling X random modified sequences from background spectronaut report
  r2arpart.dt <-
    data.frame(sample(unique(
      r2ab.dt$EG.ModifiedSequence),
      ldaBGsize,
      replace = FALSE,
      prob = NULL
    ))
  setnames(r2arpart.dt, "sample.unique.r2ab.dt.EG.ModifiedSequence...ldaBGsize..replace...FALSE..","EG.ModifiedSequence")
  r2arpart.dt <- data.frame(subset(r2arpart.dt, select=("EG.ModifiedSequence")))
  #select samples mod sequences from bg in spectronaut report
  r3atb.dt <- merge(r2aBG.dt, r2arpart.dt, by="EG.ModifiedSequence",all.x=F,all.y=T)
  #combine targets and background
  r3a.dt <- rbind(r3atb.dt, r2at.dt)
  #reconvert to data.frame containing only targets and random background
  r3.df <- data.frame(r3a.dt)
  #remove missing values
  r3.df <- subset(r3.df,!is.na(r3.df$ProtModSeq))

  
  #calculate run means for modified sequences and remove sparse modified sequences from data set
  # set as numeric
  r3.dt <- data.table(r3.df)
  
  # SN quantity fix
  if ("PEP.Quantity" %in% colnames(r3.dt)) {
    nSN1.dt$FG.NormalizedMS2PeakArea <- as.numeric(nSN1.dt$PEP.Quantity)
  } 
  if ("FG.MS2Quantity" %in% colnames(r3.dt)) {
    nSN1.dt$FG.NormalizedMS2PeakArea <- as.numeric(nSN1.dt$FG.MS2Quantity) #fix for MS2 quantity in SN15
  }
  if ("FG.NormalizedMS2PeakArea" %in% colnames(r3.dt)) {
    nSN1.dt$FG.NormalizedMS2PeakArea <- as.numeric(nSN1.dt$FG.NormalizedMS2PeakArea) #fix for MS2 quantity in older SN reports
  }
  

  if (  answerlog == "Log" ) {
    #log transform intensities
    r3.dt$log.FG.NormalizedMS2PeakArea <- log10(r3.dt$FG.NormalizedMS2PeakArea)
    #Modified sequences averaged per run
    r3.dt[, Runmean := mean(log.FG.NormalizedMS2PeakArea, na.rm = TRUE), by =
            list(RunModseq)]
  } else {
    #Modified sequences averaged per run
    r3.dt[, Runmean := mean(FG.NormalizedMS2PeakArea, na.rm = TRUE), by =
            list(RunModseq)]
  }
  r3.dt <- subset(r3.dt, !duplicated(r3.dt$RunModseq))
  
  #________________________________________________________________________________________________
  #standardize data from 0 to 1
  r3.dt[, Runmeanstd := rescale(Runmean,c(0,1)), by =
          list(EG.ModifiedSequence)]
  
  r3.dt <- data.frame(r3.dt)
  #rename for plots
  setnames(r3.dt, "Runmean", "Intensity")
  #r3.dt$R.Condition[r2.dt$R.Condition == 0] <- 0.0001
  
  #remove scarse modified sequences
  r3a.dt <- data.table(r3.dt)
  #first peptides with not all enough replicates
  r3a.dt[, Repoccurence := length(unique(R.FileName)), by =
           list(CondModseq)]
  r3a.dt <- subset(r3a.dt,r3a.dt$Repoccurence >= (max(r3a.dt$Repoccurence))/2)
  #second peptides with not all conditions
  r3a.dt[, Occurence := length(unique(ExpCond)), by =
           list(ProtModSeq)]
  r3a.dt <- subset(r3a.dt,r3a.dt$Occurence >= (max(r3a.dt$Occurence))*0.66)
  
  r3.df <- data.frame(r3a.dt)
  
  
  #_____________________________________________________________________________________________
  # Dose response fitting and correlation measurement
  #get your runs/samples/proteins of interest in a vector
  #prepare correlation data frame
  correlation.df <- subset(r3.df, select=c("EG.ModifiedSequence", "ProtModSeq"))
  correlation.df <- subset(correlation.df, !duplicated(correlation.df$ProtModSeq))
  ProtModSeq.v <- unique(correlation.df$ProtModSeq)
  setnames(correlation.df, "EG.ModifiedSequence", "Group")
  correlation.df$DRC <- 0
  correlation.df$ED50 <- 0
  correlation.df$corpvalue <- 0
  #reset state counter
  state <- c(0)
  
  #perform dose response analysis and correlation for selected modified sequences
  for (k in 1:(length(ProtModSeq.v))) {
    #get current variable
    modSeq <-  ProtModSeq.v[k]
    
    #get current iteration
    state <- state + 1
    
    
    #take valid subset of the whole data set
    r1.dt.temp <- subset(r3.df, r3.df$ProtModSeq == modSeq)
    
    tryCatch(r1fit.dt.temp <-
               drm(
                 Intensity ~ R.Condition,
                 ProtModSeq,
                 data = r1.dt.temp,
                 fct = LL.4(),
                 #logDose=10,
                 na.action = na.omit
               ), error=function(e) print("Correlation LL.4 failed, LL3.4 backup performed"))
    if(exists("r1fit.dt.temp") == TRUE) {
      cat(round((100*state/length(ProtModSeq.v)), digits = 0), "%\n")} else {
        tryCatch(r1fit.dt.temp <-
                   drm(
                     Intensity ~ R.Condition,
                     ProtModSeq,
                     data = r1.dt.temp,
                     fct = LL.3(),
                     #logDose=10,
                     na.action = na.omit
                   ), error=function(e) print("Correlation LL.3 failed, L.4 backup performed"))
      }
    if(exists("r1fit.dt.temp") == TRUE) {
      cat(" ")} else {
        tryCatch(r1fit.dt.temp <-
                   drm(
                     Intensity ~ R.Condition,
                     ProtModSeq,
                     data = r1.dt.temp,
                     fct = L.4(),
                     #logDose=10,
                     na.action = na.omit
                   ), error=function(e) print("No correlation performed"))
      }
    
    correlations <- cor(r1.dt.temp$Intensity, fitted(r1fit.dt.temp), method = answercor) # , use = "complete.obs"
    corrtest <- cor.test(r1.dt.temp$Intensity, fitted(r1fit.dt.temp), method = answercor)$p.value
    correlation.df[state, 3] <- correlations
    ED50est <- data.frame(ED(r1fit.dt.temp, 50, display = FALSE))
    correlation.df[state, 4] <- ED50est$Estimate
    correlation.df[state, 5] <- corrtest
    try(rm(r1fit.dt.temp))
    try(rm(correlations))
    try(rm(corrtest))
    try(rm(ED50est))  
  }
  
  
  
  
  ##add dose response correlations to candidate list
  t1c.dt <-
    merge(
      t1.df,
      correlation.df,
      by = c("Group"),
      all.x = F,
      all.y = T
    )
  
  #check if correlation pvalues not out off range for stable scores
  if(sum(-log(t1c.dt$corpvalue,10) > 32) > 0) {
    cat("Correlation test p-values out of range (-log(pvalue)) > 32, stopping")
    stop()
  }
  
  #_____________________________________________________________________________________________
  #Preparation for Machine learning
  t1c.dt$corpvalue[is.na(t1c.dt$corpvalue)] <- 1
  t1c.dt$PG.Coverage <- as.numeric(str_sub(t1c.dt$PG.Coverage, end=-2))
  t1c.dt$PG.Coverage[is.na(t1c.dt$PG.Coverage)] <- 0.01
  t2c.dt <- t1c.dt
  t2c.dt$DRCNorml <- -log10(t2c.dt$corpvalue) / 32
  t2c.dt <- subset(t2c.dt,!is.na(t2c.dt$Target))
  t2c.dt$PG.Coverage <- as.numeric(round(t2c.dt$PG.Coverage))
  
  #calculate normalized Qvalue
  t2c.dt$logQvalue <- -log(t2c.dt$Qvalue,10)
  #divide trough defined max 50
  t2c.dt$logQvalueNorml <- t2c.dt$logQvalue / 50
  #calculate normalized FC
  t2c.dt$FCNorml <- t2c.dt$Absolute.AVG.Log2.Ratio / max(t2c.dt$Absolute.AVG.Log2.Ratio)
  
  #Prepare LIP frequency library
  t2c.dt$normPFL <-   1-t2c.dt$PLFrequency
  t2c.dt$normPFL[is.na(t2c.dt$normPFL)] <- 1
  
  #rank all candidates per comparison
  t2c.dt <- data.table(t2c.dt)
  #t2c.dt[,Rank:=rank(Qvalue),by=list(Filter)]
  t2c.dt[,Rank:=rank(corpvalue),by=list(Filter)] #more relevant score to rank
    t2c.dt[,Freq:=length(Filter),by=list(ProteinNames)]
  #remove potential off targets from ranking
  ##subset first bg data set
  #t2c.dt <- filter(t2c.dt, Rank > 20 | DRC < 0.75 & Target == "BG")  #or use , instead of |
  #t2c.dt$Rank <- "NULL"
  
  #calculate top10% occurrence
  t1cTop10perc.dt <- subset(t2c.dt, t2c.dt$Rank <= (0.1*nrow(t2c.dt))) 
  t1cRest.dt <- subset(t2c.dt, t2c.dt$Rank > (0.1*nrow(t2c.dt)))
  t1cTop10perc.dt[,Freq:=length(Filter),by=list(ProteinNames)]
  t1cTop10perc.dt$FreqNorml <- t1cTop10perc.dt$Freq / max(t1cTop10perc.dt$Freq)
  #t1cRest.dt$Freq <- 0.001
  t1cRest.dt$FreqNorml <- 1/nrow(t2c.dt) #why do we do this??? 
  t2c.dt <- rbind(t1cTop10perc.dt,t1cRest.dt)
  t2c.dt <- data.frame(t2c.dt)
  
  #select relevant columns

  t2c.dt <- subset(t2c.dt, select=c("Target", "logQvalueNorml", "DRCNorml", "FreqNorml", "normPFL"))


  #Machine learning
  ## section4-ML training set ----
  
  #LDA analysis of the training set
  ml.df <- lda(formula = Target ~ ., 
               data = t2c.dt)
  #SVM
  #mlsvm.df <- svm(Target ~ ., data = t2c.dt)
  
  #get score weights
  scoreweights.v <- ml.df$scaling
  
  #summary ML analysis
  cat("Calculated score weights: \n")
  print(ml.df$scaling)
  cat("Score weight trainnig executed.")
  
}
#_____________________________________________________________________________________________
#predict true candidates in new data set
## section5-ML-based-ranking ----
#second drug to apply weighted scores on
#data loading
#import and merge
cat( "Selection unfiltered candidate list   \n" )
tp <- file.choose()
tp1 <- sub("", "", basename(tp))
for (i in 1:length(tp)) assign(tp1[i], read.delim(tp[i]))

#list names of data frames to merge
file_names6=(tp1)
#preprocessing of loaded data frames
#Name the files in a new column
names(file_names6) <- file_names6
n.df = ldply(file_names6, get)
#saving candidate list file name
filenameCL <- str_sub(file_names6, end=-5)
n.df$.id <- NULL

#set working directory
tempwd <- tp[1]
workingDir <- sub("", "", dirname(tempwd))
setwd(workingDir)


#comparison filtering
cat( "   Comparison filtering\n" )
#check if already renamed
if("Comparison..group1.group2." %in% colnames(n.df)) {
  setnames(n.df,"Comparison..group1.group2.","Filter")}
#get comparison column levels
filter.v <- unique(n.df$Filter)
filter.v <- as.character(filter.v)
filter.v <- c(filter.v, c("_____________"))


# selection conditions column to retain
cat( " Selection of the comparison to retain (e.g. 10 to 100x IC50).\n" )
cond1.true.v <- select.list( filter.v, multiple=T, title="Selection of a statistical comparison to be used ineh analysis (e.g. 10 to 100x IC50)" )
if ( length( cond1.true.v ) == 0 ) {
  cat( " No conditions selected Stopping.\n" )
  stop()
}
cat( " conditions selected: ", paste( cond1.true.v, collapse=", " ), "\n", sep="" )
cat( "   \n" )

#Filtering based on condition
#true
n.df <- n.df[(n.df$Filter %in% cond1.true.v), ]

#qvalue cutoff filtering
n.df <- subset(n.df, n.df$Qvalue <= Qcutoff)
#optional FC cutoff filtering
n.df <- subset(n.df, n.df$Absolute.AVG.Log2.Ratio >= FCcutoff) 

#check if correlation pvalues not out off range for stable scores
if(max(-log(n.df$Qvalue,10)) > 60) {
  cat("Qvalues out of range (-log(Qvalue)) > 60, stopping") ### check consistency with ML part # added by 
  stop() # should not be commented out 
}

#generating sub data frame
n1.df <- n.df

#ensure Qvalues are numbers and there is no NA
n1.df$Qvalue <- as.numeric(n1.df$Qvalue)
n1.df <- subset(n1.df,!is.na(n1.df$Qvalue))
n1.df <- n1.df[with(n1.df, order(Qvalue)), ]

#remove duplicated mod sequences from candidates keeping the best qvalue
n1.df <- subset(n1.df,!duplicated(n1.df$Group))

#check if protein names exist, otherwise take accessions
if(length(unique(n1.df$ProteinNames)) == 1) {
  n1.df$ProteinNames <- NULL
  n1.df$ProteinNames <- n1.df$ProteinGroups
}

#check if protein names exist, otherwise take accessions
if (!("ProteinGroups" %in% colnames(n1.df))) {
  n1.df$ProteinGroups <- n1.df$ProteinNames
}

#Protein annotation LIP_PFL
if(length(unique(n1.df$ProteinGroups)) != 1) {
  n1.df$PFL.PG <- gsub("\\;.*","",n1.df$ProteinGroups)}

#________________________________________________________________________________________________
#annotate protein frequency library in comparison table
n1.df <- merge(n1.df, pfl.df, by=c("PFL.PG"), all.x=TRUE, all.y=FALSE)
n1.df$PFL.PG <- NULL
n1.df$PLFrequency[is.na(n1.df$PLFrequency)] <- 0

#_____________________________________________________________________________________________
# selection 
cat("Please select, if you want to use a list of known drug targets.\n" )
knowntargets <- select.list( answeropt, multiple=T, title="Use of a list with know drug targets (Uniprot Protein Names), default NO" )
if ( ( length(knowntargets) == 0 ) ) {
  cat("No table with konwn drug targets will be imported.\n" )
  knowntargets <- "NO"
} else { cat("Your answer was  ", paste( knowntargets, collapse=", "), "\n", sep="" )
}

if ( knowntargets == "YES" ) {
  #load known targets
  #import and merge
  cat( "Select list of known drug targets   \n" )
  tempr <- file.choose()
  tempr1 <- sub("", "", basename(tempr))
  for (i in 1:length(tempr)) assign(tempr1[i], read.delim(tempr[i]))
  
  #list names of data frames to merge
  file_names5=(tempr1)
  #preprocessing of loaded data frames
  #Name the files in a new column
  names(file_names5) <- file_names5
  lookupntargets.df = ldply(file_names5, get)
  
  lookupntargets.df$.id <- NULL
  #annotate know targets in candidates
  n1.df <- merge(n1.df, lookupntargets.df, by="ProteinNames", all.x= TRUE)
  n1.df$Target <- as.character(n1.df$Target)
  n1.df$Target[is.na(n1.df$Target)] <- "BG"
}

#prefilter candidate list by qvalue and FC cutoff user defined
n1s.df <- subset(n1.df, n1.df$Qvalue <= Qcutoff)
n1s.df <- subset(n1s.df, n1s.df$Absolute.AVG.Log2.Ratio >= FCcutoff)
n1s.df$Prefilter <- "Pass"
n1s.df <- subset(n1s.df, select=c("Group", "Prefilter"))
setnames(n1s.df,"Group","EG.ModifiedSequence")



#________________________________________________________________________________________________
#annotate tryptic, semitryptic peptides in comparison table 
#change folder to directory where PFL file is located or select manually, during running of the script
# try(ISdigest.df  <- read.delim("S:/Big_Temp/Roland/LIP Resource/uniprot_sprot_2018-07-01_HUMAN_ISOFORMS_digest.tsv", header=T, sep="\t", stringsAsFactors=F, na.strings=c( "n. def." ) ))
# if (exists("ISdigest.df") == TRUE) {
#   ISdigest.df$Semitryptic <- "NO"
#   ISdigest.df$ProteinId <- NULL
#   ISdigest.df$PeptideLength <- NULL
#   ISdigest.df$NrDigestPeptides <- NULL
#   setnames(ISdigest.df,"PeptideSequence","EG.StrippedSequence")
#   #generate stripped sequence in comparison table
#   #saving candidate list file name
#   n1.df$EG.StrippedSequence <- n1.df$Group
#   #remove underscores
#   n1.df$EG.StrippedSequence <- str_sub(n1.df$EG.StrippedSequence, end=-2)
#   n1.df$EG.StrippedSequence <- str_sub(n1.df$EG.StrippedSequence, start=2)
#   #remove cmids
#   n1.df$EG.StrippedSequence <-str_replace_all(n1.df$EG.StrippedSequence, "\\[\\+C2\\+H3\\+N\\+O\\]", "")
#   n1.df$EG.StrippedSequence <-str_replace_all(n1.df$EG.StrippedSequence, "\\[\\+C2\\+H2\\+O\\]", "")
#   n1.df$EG.StrippedSequence <-str_replace_all(n1.df$EG.StrippedSequence, "\\[\\+O\\]", "")
#   n1.df$EG.StrippedSequence <-str_replace_all(n1.df$EG.StrippedSequence, "\\[\\+O\\-H\\-N\\]", "")
#   n1.df$EG.StrippedSequence <-str_replace_all(n1.df$EG.StrippedSequence, "\\[\\+C6\\+H10\\+O5\\]", "")
#   n1.df$EG.StrippedSequence <-str_replace_all(n1.df$EG.StrippedSequence, "\\[\\+H\\+O3\\+P\\]", "")
#   #[+C6+H10+O5]
#   #[+H+O3+P]
#   #convert precursors into stripped sequences alternative approch
#   #t.df$Sequence <-str_replace_all(t.df$Precursor, "\\[.*?\\]", "")
#   #t.df$Sequence <-str_replace_all(t.df$Sequence, "\\_", "")
#   #t.df$Sequence <-str_replace_all(t.df$Sequence, "\\..*", "")
#   
#   
#   #annotate candidate list
#   n1.df <- merge(n1.df, ISdigest.df, by = c("EG.StrippedSequence"), all.x=TRUE, all.y=FALSE)
#   n1.df$Semitryptic[is.na(n1.df$Semitryptic)] <- "YES"
# }

#_____________________________________________________________________________________________
#load SN report for data set where the ML should be applied to
#data loading
#import and merge
cat( "Select SN report for ML based ranking  \n" )
tmp <- file.choose()
tmp1 <- sub("", "", basename(tmp))
for (i in 1:length(tmp)) assign(tmp1[i], read.delim(tmp[i]))

#list names of data frames to merge
file_names4=(tmp1)
#preprocessing of loaded data frames
#Name the files in a new column
names(file_names4) <- file_names4
nSN.df = ldply(file_names4, get)
#remove unnecessary column
nSN.df$.id <- NULL

nSN.df <- subset(nSN.df, !grepl("iRT", nSN.df$PG.ProteinAccessions))

#check if external format
if ("EG.ModifiedPeptide" %in% colnames(nSN.df)) {
  setnames(nSN.df,"EG.ModifiedPeptide","EG.ModifiedSequence") }

#check if protein names exist, otherwise take accessions
if(length(unique(nSN.df$PG.ProteinNames)) == 3) {   ### why 3 now? 
  nSN.df$PG.ProteinNames <- NULL
  setnames(nSN.df,"PG.ProteinAccessions","ProteinNames")
} else {
  setnames(nSN.df,"PG.ProteinNames","ProteinNames")
}

#________________________________________________________________________________________________
#annotate SN report semitryptic peptides #cannot be done without the previous annotation
if (exists("ISdigest.df") == TRUE) {
  nSN.df <- merge(nSN.df, ISdigest.df, by = c("EG.StrippedSequence"), all.x=TRUE, all.y=FALSE)
  nSN.df$Semitryptic[is.na(nSN.df$Semitryptic)] <- "YES"
  nSN.df$Count <- 1
  
  tryp.df <- table(nSN.df$Semitryptic)
  NO.df <- as.vector(tryp.df[names(tryp.df)=="NO"])
  
  #test if more than 5% tryptic peptides, otherwise do not annotate
  if ((nrow(nSN.df)*0.05 < NO.df) == TRUE) {
    
    #tryptic semitryptic plots
    nSNplot.dt<- data.table(nSN.df)
    nSNplot.dt$Count <- 1
    #count semi/tryptic per run
    STrypticplot.df <- nSNplot.dt[ , .(Totalcount = sum(Count)), by = .(R.FileName, Semitryptic)]
    #shorten RAW file names
    STrypticplot.df$R.FileName <- str_sub(STrypticplot.df$R.FileName, start=20)
    STrypticplot.df$R.FileName <- str_sub(STrypticplot.df$R.FileName, end=-13)
    
    #order truns
    #STrypticplot.df$R.FileName <- as.numeric(STrypticplot.df$R.FileName)
    STrypticplot.df <- STrypticplot.df[with(STrypticplot.df, order(R.FileName)), ]
    
    
    #write semi trypric table
    write.table(STrypticplot.df, "Semitryptic-analysis-summary.txt", sep="\t", row.names = FALSE, quote = FALSE)
    
    #plot semi/trypric
    #with colorbrewer
    s.tryp.plot <- ggplot(STrypticplot.df, aes(x = R.FileName, fill= Semitryptic, y = Totalcount)) +
      geom_bar(colour="black", position="stack", stat = "identity") +
      theme(text = element_text(size=20)) + 
      theme(axis.text.x = element_text(colour="black",size=16)) +
      theme(axis.text.y = element_text(colour="black",size=16)) +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
    
    pdf(paste(filenameCL, "-semitrypic-ana.pdf", sep=""))
    plot(s.tryp.plot)
    dev.off()
    rm(s.tryp.plot)
    
  }} else {n1.df$Semitryptic <- NULL}
#_________________________________________________________________________________________________

#________________________________________________________________________________________________
#getting protein coverages from Spectronaut report file
if ("PG.Coverage" %in% colnames(nSN.df) && "PEP.IsProteotypic" %in% colnames(nSN.df)) {
  nc.df <- subset(nSN.df, !duplicated(nSN.df$EG.ModifiedSequence))
  nc.df <- subset(nc.df, select= c(EG.ModifiedSequence, PG.Coverage, PEP.IsProteotypic))
  #add protein coverage to 2nd candidate file
  n1.df$EG.ModifiedSequence <- n1.df$Group
  n1a.df <- merge(n1.df, nc.df, by=c("EG.ModifiedSequence"), all.x=TRUE, all.y=FALSE)
  n1.df <- subset(n1a.df, !duplicated(n1a.df$Group))}


if ("PG.Coverage" %in% colnames(nSN.df) && !("PEP.IsProteotypic" %in% colnames(nSN.df))) {
  nc.df <- subset(nSN.df, !duplicated(nSN.df$EG.ModifiedSequence))
  nc.df <- subset(nc.df, select= c(EG.ModifiedSequence, PG.Coverage))
  #add protein coverage to 2nd candidate file
  n1.df$EG.ModifiedSequence <- n1.df$Group
  n1a.df <- merge(n1.df, nc.df, by=c("EG.ModifiedSequence"), all.x=TRUE, all.y=FALSE)
  n1.df <- subset(n1a.df, !duplicated(n1a.df$Group))}

if (!("PG.Coverage" %in% colnames(nSN.df)) && "PEP.IsProteotypic" %in% colnames(nSN.df)) {
  nc.df <- subset(nSN.df, !duplicated(nSN.df$EG.ModifiedSequence))
  nc.df <- subset(nc.df, select= c(EG.ModifiedSequence, PEP.IsProteotypic))
  #add protein coverage to 2nd candidate file
  n1.df$EG.ModifiedSequence <- n1.df$Group
  n1a.df <- merge(n1.df, nc.df, by=c("EG.ModifiedSequence"), all.x=TRUE, all.y=FALSE)
  n1.df <- subset(n1a.df, !duplicated(n1a.df$Group))}

DRCtesting <- "YES"
#____________________________________________________________________________________________

if (  DRCtesting == "YES" ) {
#extract protein candidates with qvalue (and FC) prefiltered from SN report for DRC analysis
nSN.df <- merge(nSN.df, n1s.df, by="EG.ModifiedSequence", all.x=FALSE, all.y= TRUE)
#prepare combined columns

nSN.df$ProteinNamestrm <- strtrim(nSN.df$ProteinNames, 25)

nSN.df$ProtModSeqtrm <-
  do.call(paste, c(nSN.df[c("ProteinNamestrm", "EG.ModifiedSequence")], sep = ""))
nSN.df$ProtModSeq <-
  do.call(paste, c(nSN.df[c("ProteinNames", "EG.ModifiedSequence")], sep = ""))
nSN.df$RunModseq <-
  do.call(paste, c(nSN.df[c("R.FileName", "EG.ModifiedSequence")], sep = ""))
nSN.df$CondModseq <-
  do.call(paste, c(nSN.df[c("R.Condition", "EG.ModifiedSequence")], sep = ""))


#________________________________________________________________________________________________
# selection annotation of conditions as concentrations
#conditions
#remove runs with no condition
nSN.df <- nSN.df[complete.cases(nSN.df$R.Condition), ]

answeropt1 <- c("YES", "NO","Retain from SN report", "_____________________________________________________")

cat("Please select, if you want to use a concentration table.\n" )
knownconc <- select.list( answeropt1, multiple=T, title="Please select, if you want to import a concentration table (YES) or generate standard concentrations (NO).Default is retain)" )
if ( length(knownconc) == 0 ) {
  cat("No concentration file import, Spectronaut report conditions/concentrations will be taken .\n" )
  knownconc <- "Retain from SN report"
} else { cat("Your answer was  ", paste( knownconc, collapse=", "), ". A concentration table as txt was generated in the analysis directory. Please fill in the generated table.\n", sep="" )
  condsimport.df <- data.frame(unique(nSN.df$R.Condition))
  setnames(condsimport.df, "unique.nSN.df.R.Condition." , "R.Condition")
  condsimport.df <- data.frame(condsimport.df[with(condsimport.df, order(R.Condition)), ])
  condsimport.df$Concentrations <- "fill-in"
  setnames(condsimport.df, "condsimport.df.with.condsimport.df..order.R.Condition....." , "R.Condition")
  write.table(condsimport.df, "Concentration-table.txt", sep="\t", row.names = FALSE, quote = FALSE)
}

if ( knownconc == "YES" ) {
  #lookup table (columns R.Condition, Concentrations)
  cat( "Please provide a conditions concentration table\n" )
  tempc <- file.choose()
  tempc1 <- sub("", "", basename(tempc))
  for (i in 1:length(tempc)) assign(tempc1[i], read.delim(tempc[i]))
  
  #list names of data frames to merge
  file_names7=(tempc1)
  #preprocessing of loaded data frames
  #Name the files in a new column
  names(file_names7) <- file_names7
  conds1.v = ldply(file_names7, get)
  conds1.v$.id <- NULL
  conds1.v$Concentrations <- as.numeric(conds1.v$Concentrations)
  
}

if ( knownconc == "NO" ) {
  ##specific for experiments without conc annotation in conditions but plain numbering /letters
  if (exists("PosContrl") == TRUE) {
    nSN.df <- subset(nSN.df, nSN.df$R.Condition != PosContrl )}
  conds1.v <- data.frame(unique(nSN.df$R.Condition))
  conds1.v$a <- "a"
  setnames(conds1.v, "unique.nSN.df.R.Condition." , "R.Condition")
  conds1.v <- conds1.v[with(conds1.v, order(R.Condition)), ]
  ncond <- length(unique(nSN.df$R.Condition)) -1
  conds1.v$Concentrations <- c(0, lseq(1,10^(ncond-1), length = ncond))
  conds1.v$a <- NULL
}


#replace condition annotation
if ( knownconc != "Retain from SN report" ) {
  nSN.df <- merge(nSN.df, conds1.v, by=c("R.Condition"), all.x=FALSE, all.y=TRUE)
  nSN.df$R.Condition <- NULL
  setnames(nSN.df,"Concentrations","R.Condition")
}


#________________________________________________________________________________________________
answeropt2 <- c("nM", "pM","fM", "aM", "_____________________________________________________")

cat("Please select the concentration unit (default fM).\n" )
concUnit <- select.list( answeropt2, multiple=T, title="Please select the concentration unit (default nM)" )
if ("_____________________________________________________" == concUnit) {
  cat("No concentration file import, Spectronaut report conditions/concentrations will be taken .\n" )
  concUnit <- "nM"
} else { cat("Your answer was  ", paste( concUnit, collapse=", "), ".\n", sep="" )
}
#answerlog <- "Log"

#________________________________________________________________________________________________
# calculate means
nSN1.dt <- data.table(nSN.df)
# set as numeric
#check if PEP.Quantity exists, otherwise use alternative column names
if ("PEP.Quantity" %in% colnames(nSN.df)) {
  nSN1.dt$PEP.Quantity <- as.numeric(nSN1.dt$PEP.Quantity)
} 
if ("FG.MS2Quantity" %in% colnames(nSN.df)) {
  nSN1.dt$PEP.Quantity <- as.numeric(nSN1.dt$FG.MS2Quantity) #fix for MS2 quantity in SN15
}
if ("FG.NormalizedMS2PeakArea" %in% colnames(nSN.df)) {
  nSN1.dt$PEP.Quantity <- as.numeric(nSN1.dt$FG.NormalizedMS2PeakArea) #fix for MS2 quantity in older SN reports
}

if (  answerlog == "Log" ) {
  #log transform intensities
  nSN1.dt$log.PEP.Quantity <- log10(nSN1.dt$PEP.Quantity)
  #calculate means based on the log10 intensities
  nSN1.dt[, Runmean := mean(log.PEP.Quantity, na.rm = TRUE), by =
            list(RunModseq)]
} else {
  #calculate means for modified peptides per run based  (charge states)
  nSN1.dt[, Runmean := mean(PEP.Quantity, na.rm = TRUE), by =
            list(RunModseq)]
}
# aggregation of charge states? 

#remove duplicated mod sequences per run
nSN2.dt <- subset(nSN1.dt,!duplicated(nSN1.dt$RunModseq))


#standardize data from 0 to 1
nSN2.dt[, Runmeanstd := rescale(Runmean,c(0,1)), by =
          list(EG.ModifiedSequence)]

setnames(nSN2.dt, "Runmeanstd", "Intensity")


#__________________________________________________________________________________________________
#remove scarse modified sequences
#first peptides with not all enough replicates
nSN2.dt[, Repoccurence := length(unique(R.FileName)), by =
          list(CondModseq)]
nSN2.dt <- subset(nSN2.dt,nSN2.dt$Repoccurence >= (max(nSN2.dt$Repoccurence)*0.25)) 
#second peptides with not all conditions
nSN2.dt[, Occurence := length(unique(R.Condition)), by =
          list(EG.ModifiedSequence)]
nSN2.dt <- subset(nSN2.dt,nSN2.dt$Occurence >= (max(nSN2.dt$Occurence)-2))

#remove excluded concentrations
if (exists("RemConc") == TRUE) {
  nSN2.dt <- subset(nSN2.dt, nSN2.dt$R.Condition != RemConc )
  }

#__________________________________________________________________________________________________
#Correlation analysis
#get your runs/samples/proteins of interest in a vector
modSeq1.v <- unique(nSN2.dt$EG.ModifiedSequence)
ProtModSeq2.v <- unique(nSN2.dt$ProtModSeq)
ProtModSeq1.v <- unique(nSN2.dt$ProtModSeqtrm)
NumbCond <- length(unique(nSN2.dt$R.Condition))-1

#prepare correlation data frame
correlation1.df <- data.frame(modSeq1.v, ProtModSeq2.v)
setnames(correlation1.df, "modSeq1.v", "Group")
correlation1.df$DRC <- 0
correlation1.df$ED50 <- 0
correlation1.df$corpvalue <- 0
correlation1.df$model <- "" #model info added - to check which models are fitted

#reset correlation iteration counter
state <- c(0)

#generate subdirectory for plots
dir.create("DRCplots")


#perform dose response analysis and correlation for selected modified sequences
for (k in 1:(length(ProtModSeq1.v))) {
  #get current variable
  modSeq <-  ProtModSeq1.v[k]
  #get current iteration
  state <- state + 1
  
  
  #take valid subset of the whole data set
  r1.dt.temp <- subset(nSN2.dt, nSN2.dt$ProtModSeqtrm == modSeq)
  
  #perform dose response fit
  tryCatch(r1fit.dt.temp <-
             drm(
               Intensity ~ R.Condition,
               ProtModSeq,
               data = r1.dt.temp,
               fct = LL.4(),
               #logDose=10,
               na.action = na.omit
             ), error=function(e) print("Correlation LL.4 failed, LL.3 backup performed"))
  if(exists("r1fit.dt.temp") == TRUE) {
    cat("Progres: ")
    cat(round((100*state/length(ProtModSeq1.v)), digits = 0), "%\n")} else {
      tryCatch(r1fit.dt.temp <-
                 drm(
                   Intensity ~ R.Condition,
                   ProtModSeq,
                   data = r1.dt.temp,
                   fct = LL.3(),
                   #logDose=10,
                   na.action = na.omit
                 ), error=function(e) print("Correlation LL.3 failed, L.4 backup performed"))
    }
  if(exists("r1fit.dt.temp") == TRUE) {
    cat(" ")} else {
      tryCatch(r1fit.dt.temp <-
                 drm(
                   Intensity ~ R.Condition,
                   ProtModSeq,
                   data = r1.dt.temp,
                   fct = L.4(),
                   #logDose=10,
                   na.action = na.omit
                 ), error=function(e) print("No correlation performed"))
    }
  #perform correlation analysis data to fit (pearson)
  correlations <- cor(r1.dt.temp$Intensity, fitted(r1fit.dt.temp), method = answercor) # , use = "complete.obs"
  corrtest <- cor.test(r1.dt.temp$Intensity, fitted(r1fit.dt.temp), method = answercor)$p.value
  
  correlation1.df[state, 3] <- correlations
  
  # Calculate ED50 (IC50)
  ED50est <- data.frame(ED(r1fit.dt.temp, 50, display = FALSE))
  correlation1.df[state, 4] <- ED50est$Estimate
  correlation1.df[state, 5] <- corrtest
  correlation1.df[state, 6] <- if(is.character(r1fit.dt.temp)){
    r1fit.dt.temp
  } else {
    r1fit.dt.temp[["fct"]][["name"]]
  } #ADDED
  
  #Plot the data
  pdf(paste("DRCplots/Cent-Plot-",modSeq,".pdf", sep=""))
  par(mar = c(7, 7, 7, 5))
  plot <- plot( r1fit.dt.temp,  broken = TRUE, ylab= " ", xlab= paste("Concentration (", concUnit, ")", sep = ""))
  title(ylab=paste(answerlog, "-response (scaled)", sep = ""), line=5)
  title(main = paste(modSeq," corr=", round(correlations,digits = 3), sep=" "))
  dev.off()
  
  #Plot the data
  pdf(paste("DRCplots/Reps-Plot-",modSeq,".pdf", sep=""))
  par(mar = c(7, 7, 7, 5))
  plot <- plot( r1fit.dt.temp,  broken = TRUE, type = "all", ylab= " ", xlab= paste("Concentration (", concUnit, ")", sep = ""))
  title(ylab=paste(answerlog, "-response (scaled)", sep = ""), line=5)
  title(main = paste(modSeq," corr=", round(correlations,digits = 3), sep=" "))
  dev.off()
  
  #remove temporary objects
  try(rm(ED50est))
  try(rm(r1fit.dt.temp))
  try(rm(correlations))
  try(rm(corrtest))
  try(rm(plot))
}


#all_models <- correlation1.df
#l.4_model <- correlation1.df

#correlation1.df <- all_models

#add dose response correlations to candidate list
#prepare for application of weights to scores

n1c.df <-
  merge(
    n1.df,
    correlation1.df,
    by = c("Group"),
    all.x = F,
    all.y = T
  )

n1c.df$corpvalue[is.na(n1c.df$corpvalue)] <- 1
#check if correlation pvalues not out off range for stable scores
if(max(-log(n1c.df$corpvalue,10)) > 30) {
  cat("Correlation p-values out of range (-log(pvalue)) > 30, stopping")
  #stop()
}

#qvalue of correlation test
qvalues <- qvalue(p =  n1c.df$corpvalue)
n1c.df <- cbind( n1c.df, qvalues$qvalues)
setnames(n1c.df, "qvalues$qvalues", "DRCqvalue")
#prepare for application of weights to scores
n1c.df <- n1c.df[with(n1c.df, order(corpvalue)), ]

#n1c.df$qval <- p.adjust(n1c.df$corpvalue, method = "BH")

#normalize DRC
n1c.df$logcorpval <- -log(n1c.df$corpvalue,10)
#stabilization of LIP score to max of 30 for logcorpval (never observed higher pvalue)
n1c.df$DRCNorml <- n1c.df$logcorpval / 30
#ensure not larger than 1
n1c.df$DRCNorml[n1c.df$DRCNorml > 1] <- 1

} else{
  #subset data frame
  n1c.df <- n1.df
  #generate empty DRC score
  n1c.df$DRCNorml <- 0
    #order by qvalue
  n1c.df <- n1c.df[with(n1c.df, order(Qvalue)), ]  
}

#qvalue of t-test
n1c.df$logQvalue <- -log(n1c.df$Qvalue,10)
#pvalue stabilization to max 40 (above highest detected)
n1c.df$logQvalueNorml <- n1c.df$logQvalue / 60 
#ensure not larger than 1
n1c.df$logQvalueNorml[n1c.df$logQvalueNorml > 1] <- 1



n1c.dt <- data.table(n1c.df)
n1c.dt[,QvalueRank:=rank(Qvalue),by=list(Filter)]

#calculate Top10% occurence
n1cTop10perc.dt <- subset(n1c.dt, n1c.dt$QvalueRank <= (0.1*nrow(n1c.dt)))
n1cRest.dt <- subset(n1c.dt, n1c.dt$QvalueRank > (0.1*nrow(n1c.dt)))
n1cTop10perc.dt[,Freq:=length(Filter),by=list(ProteinNames)]
n1cTop10perc.dt$FreqNorml <- n1cTop10perc.dt$Freq / max(n1cTop10perc.dt$Freq)
n1cTop10perc.dt$Freq <- NULL
#n1cRest.dt$Freq <- 0
n1cRest.dt$FreqNorml <- 1/nrow(n1c.dt)
n1c.dt <- rbind(n1cTop10perc.dt,n1cRest.dt)
n1c.dt <- data.frame(n1c.dt)

#Prepare LIP frequency library score
if (  (PFLselection) != "NO" ) { 
n1c.dt$normPFL <-   1-n1c.dt$PLFrequency
n1c.dt$normPFL[is.na(n1c.dt$normPFL)] <- 1
} else {
  n1c.dt$normPFL <- 1
}

if ( (PFLselection) == "NO" ) {
  scoreweights.v[4] <- 0}
if (  DRCtesting == "NO" ) {
  scoreweights.v[2] <- 0}

#normalize result of total LIP-score to six
#initoal total lip score sum
ini.sum.LIPscore.v <- sum(scoreweights.v)
#adapt score to a combined max of 6
scoreweightsmax6.v <- scoreweights.v/ini.sum.LIPscore.v*6
#check back if 6 is the max now
final.sum.LIPscore.v <- sum(scoreweightsmax6.v)

#calculate combined score
n1c.dt$LIPscore <- n1c.dt$logQvalueNorml*scoreweightsmax6.v[1] + n1c.dt$DRCNorml*scoreweightsmax6.v[2] + n1c.dt$FreqNorml*scoreweightsmax6.v[3] + n1c.dt$normPFL*scoreweightsmax6.v[4]
#prevent NAs
n1c.dt$LIPscore[is.na(n1c.dt$LIPscore)] <- 0
#order by LipScore
n1c.dt <- n1c.dt[with(n1c.dt, order(LIPscore)), ]

#ranking by combined score
#convert to data.table
n1c.dt <- data.table(n1c.dt)
#generate constant column
n1c.dt$c <- "a"
#rank all candidates
n1c.dt[,LIPRank:=(length(LIPscore)-rank(LIPscore)+1),by=list(c)]
n1c.dt[,QvalueRank:=rank(Qvalue),by=list(c)]
#remove constant column again
n1c.dt$c <- NULL


#for protein level ranking calculate max LIPscore
n1c.dt[,ProteinLipScore:=max(LIPscore),by=list(ProteinNames)]

n1c.df <- data.frame(n1c.dt)
n1c.df <- n1c.df[with(n1c.df, order(-LIPscore)), ]

#write final table
write.table(n1c.df, paste(filenameCL , "-ModSeq-ranked.txt", sep = ""), sep="\t", row.names = FALSE, quote = FALSE)

#generate constant column
n1c.df$c <- "a"
#retain best LIP scoring peptide per protein
n2c.df <- subset(n1c.df,!duplicated(n1c.df$ProteinNames))
#rank all candidates on protein level
n2c.dt <- data.table(n2c.df)
n2c.dt[,LIPRank:=(length(LIPscore)-rank(LIPscore)+1),by=list(c)]
n2c.df <- data.frame(n2c.dt)
#remove constant column again
n2c.df$c <- NULL

#write final table colapsed to protein
write.table(n2c.df, paste(filenameCL , "-Protein-ranked.txt", sep = ""), sep="\t", row.names = FALSE, quote = FALSE)

if ( knowntargets == "YES" ) {
  #plot LIP-score histogram modified peptides
  LIPScore.plot <- ggplot(n1c.df, aes(x = LIPscore, fill= Target)) +
    scale_fill_brewer(palette = "Paired")+
    xlim(c(0), 6) +
    geom_rect(aes(xmin=0, xmax=2, ymin=-10, ymax=Inf), fill="Red") +
    geom_rect(aes(xmin=2, xmax=2.5, ymin=-10, ymax=Inf), fill="Orange") +
    geom_rect(aes(xmin=2.5, xmax=4, ymin=-10, ymax=Inf), fill="Yellow") +
    geom_rect(aes(xmin=4, xmax=6, ymin=-10, ymax=Inf), fill="Green") +
    geom_histogram(aes(fill=Target)) +
    theme(text = element_text(size=20)) + 
    theme(axis.text.x = element_text(colour="black",size=16)) +
    theme(axis.text.y = element_text(colour="black",size=16)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  pdf(paste(filenameCL, "-ModSeq-Hist.pdf", sep=""))
  plot(LIPScore.plot)
  dev.off()
  rm(LIPScore.plot)
  

  #plot LIP-score histogram protein
  LIPScore.plot <- ggplot(n2c.df, aes(x = LIPscore, fill= Target)) +
    scale_fill_brewer(palette = "Paired")+
    xlim(c(0), 6) +
    geom_rect(aes(xmin=0, xmax=2, ymin=-10, ymax=Inf), fill="Red") +
    geom_rect(aes(xmin=2, xmax=2.5, ymin=-10, ymax=Inf), fill="Orange") +
    geom_rect(aes(xmin=2.5, xmax=4, ymin=-10, ymax=Inf), fill="Yellow") +
    geom_rect(aes(xmin=4, xmax=6, ymin=-10, ymax=Inf), fill="Green") +
    geom_histogram() +
    theme(text = element_text(size=20)) + 
    theme(axis.text.x = element_text(colour="black",size=16)) +
    theme(axis.text.y = element_text(colour="black",size=16)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  pdf(paste(filenameCL, "-Protein-Hist.pdf", sep=""))
  plot(LIPScore.plot)
  dev.off()
  rm(LIPScore.plot) 
} else {
  #plot LIP-score histogram modified peptides
  LIPScore.plot <- ggplot(n1c.df, aes(x = LIPscore)) +
    scale_fill_brewer(palette = "Paired")+
    xlim(c(0), 6) +
    geom_rect(aes(xmin=0, xmax=2, ymin=-10, ymax=Inf), fill="Red") +
    geom_rect(aes(xmin=2, xmax=2.5, ymin=-10, ymax=Inf), fill="Orange") +
    geom_rect(aes(xmin=2.5, xmax=4, ymin=-10, ymax=Inf), fill="Yellow") +
    geom_rect(aes(xmin=4, xmax=6, ymin=-10, ymax=Inf), fill="Green") +
    geom_histogram() +
    theme(text = element_text(size=20)) + 
    theme(axis.text.x = element_text(colour="black",size=16)) +
    theme(axis.text.y = element_text(colour="black",size=16)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  pdf(paste(filenameCL, "-ModSeq-Hist.pdf", sep=""))
  plot(LIPScore.plot)
  dev.off()
  rm(LIPScore.plot)
  
  

  #plot LIP-score histogram protein
  LIPScore.plot <- ggplot(n2c.df, aes(x = LIPscore)) +
    scale_fill_brewer(palette = "Paired")+
    xlim(c(0), 6) +
    geom_rect(aes(xmin=0, xmax=2, ymin=-10, ymax=Inf), fill="Red") +
    geom_rect(aes(xmin=2, xmax=2.5, ymin=-10, ymax=Inf), fill="Orange") +
    geom_rect(aes(xmin=2.5, xmax=4, ymin=-10, ymax=Inf), fill="Yellow") +
    geom_rect(aes(xmin=4, xmax=6, ymin=-10, ymax=Inf), fill="Green") +
    geom_histogram() +
    theme(text = element_text(size=20)) + 
    theme(axis.text.x = element_text(colour="black",size=16)) +
    theme(axis.text.y = element_text(colour="black",size=16)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  pdf(paste(filenameCL, "-Protein-Hist.pdf", sep=""))
  plot(LIPScore.plot)
  dev.off()
  rm(LIPScore.plot)
}

##clear work space
#rm(list = ls(all = TRUE))


cat(" Dose response plots saved in sub-folder.   \n")
cat(" Script has concluded.   \n")
cat("   \n")

