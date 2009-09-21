`InitialProcessing` <-
function( StructIn){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# INITIALPROCESSING Inital Prep of Data from Signal Pr0cessing 
# 
# DESCRIPTION
#       Takes peaklists that have been imported into MATLAB and prepares 
#       them for Bayesian Analysis.
# 
# USAGE
#       StructOut = InitialProcessing( StructIn)
# 
# INPUTS
#       Structure with the following double-typed arrays
#       Intensities: n x m real-valued array with variables (peaks) in
#           columns, cases (samples) in rows.
#       MZ: List of the labels (m/z value) for each of the variables.
#           Must be the same size as the number of variables in Intensities
#       Class: Classification of each sample (disease state)-- 1 or 2--must
#       be the same size as the number of cases in Intensities
#       ID: Case or patient ID number, same size as class.  May have second
#           column, so each row is [ID1 ID2} where ID2 is replicate number.
#       Options (logical):  Array of processing options with elements:
#           1. Normalize
#           2. Clip Data (remove negatives)
#           3. Replicate Average
#           4. Auto threshold MI
#           5. Use Log of Data
#           6. Remove Low Signal cases(removed in R package by Qian)
#           NOT DONE: 3 Bin (2 Bin if False)
# 
# OUTPUTS
#
#       DataStructure: MATLAB data structure with the following components:
#           RawData: Intensities as input
#           ClipData: RawData where all values less than 1 are set to 1
#           NormData: ClipData normalized by total ion count, i.e.
#               divided by the sum of all variables for each case
#           LogData: Natural logarithm of NormData
#           Class, MZ: Same as input
#           ID: SIngle column. If replicates are not averaged, the entries
#               are now ID1.ID2. If replicates averaged, then just ID1 
#           DeltaMZ: difference in peak m/z values to look for adducts
#           RatioMZ: ratios of m/z values ot look for satellites
# 
# CALLED FUNCTIONS
# 
#       None. (cdfplot is MATLAB "stat" toolbox)
 
 
## Initialize  Data
#  find the size, create the output structure,and transfer info

 
rows<-nrow(StructIn$Intensities)
cols<-ncol(StructIn$Intensities)
StructOut<-list(Intensities=StructIn$Intensities,MZ=StructIn$MZ,ID=StructIn$ID,Class=StructIn$Class,Options=StructIn$Options,n=StructIn$n,Repeats=StructIn$Repeats,Drop=StructIn$Drop,Threshold=StructIn$Threshold,RawData=StructIn$Intensities)
## Option 2: Clip Negatives from data
#  set values below 0 to be 1 because negative
#   molecule counts are not physically reasonable
# 1 is chosen rather than 0 in case log(data) is used
# Note: the decision to do this before normalization was based on
# discussions with Dr. William Cooke, who created the data set.
 
if (StructOut$Options[2]){
    StructOut$Intensities[which(StructOut$Intensities<1)]<-1 }##ok<FNDSB>

 
##  Option 6: Removal of Cases with Low Signal
#   find the sum of all values for eah row, then normalize each row to
#   account for the effects of signal strenght over time and other
#   instrumental variations in total strength of the signal
 
# Find the total ion count for each case, then the global average.
# Determine a correction factor for each case (NormFactor)
NormFactor<-NULL  ###To make a record for the items in the StructOut (by Qian)
if (StructOut$Options[1]){
    RowTotalIonCount<- rowSums(StructOut$Intensities)
    AvgTotalIonCount<-mean(RowTotalIonCount) #Population average
    NormFactor<-AvgTotalIonCount/RowTotalIonCount #Vector of norm factors
StructOut<-list(Intensities=StructOut$Intensities,MZ=StructOut$MZ,ID=StructOut$ID,Class=StructOut$Class,Options=StructOut$Options,n=StructOut$n,Repeats=StructOut$Repeats,Drop=StructOut$Drop,Threshold=StructOut$Threshold,RawData=StructOut$RawData, NormFactor=NormFactor)
                           }
# If Remove Low Signal is desired, interact with user to determine
# threshold, then remove all cases that are below the threshold.
 
 
## Option 3: Replicate Average
# This option causes cases with same ID numbers to be averaged, peak by
# peak.
 
if (StructOut$Options[3]){ #Replicate Average
                           # Collapse to unique IDs only, throw out replicate ID column
             if(length(NormFactor)!=0){ StructOut<-list(Intensities=StructOut$Intensities,MZ=StructOut$MZ,ID=StructOut$ID,Class=StructOut$Class,Options=StructOut$Options,n=StructOut$n,Repeats=StructOut$Repeats,Drop=StructOut$Drop,Threshold=StructOut$Threshold,RawData=StructOut$RawData, NormFactor= StructOut$NormFactor, Replicate_ID=StructOut$ID, Replicate_Class=StructOut$Class)
                                      }else{  
                                            StructOut<-list(Intensities=StructOut$Intensities,MZ=StructOut$MZ,ID=StructOut$ID,Class=StructOut$Class,Options=StructOut$Options,n=StructOut$n,Repeats=StructOut$Repeats,Drop=StructOut$Drop,Threshold=StructOut$Threshold,RawData=StructOut$RawData, Replicate_ID=StructOut$ID,Replicate_Class=StructOut$Class)}
                                            newID<-sort(unique(StructOut$ID[,1])) # List of unique IDs
                                            num<-length(newID) #how many are there?
                                            newClass<-matrix(0,num,1) # Holders for extracted class, data
                                            newData<-matrix(0,num,cols)
                                            for (i in 1:num){ # for each unique ID
                                                             id<-newID[i] # work on this one
                                                             cases<-which(StructOut$ID[,1]==id) # Get a list of cases with this ID
                                                             newClass[i]<-StructOut$Class[cases[1]] # save their class
                                                             casedata<-StructOut$Intensities[cases,] # get their data
                                                             if(length(cases)>1){
                                                                                newData[i,]<-colMeans(casedata) # and save the average
                                                                                } else{ newData[i,]<-mean(casedata)}
                                             }
              StructOut$Intensities<-newData
              StructOut$Class<-newClass
              StructOut$ID<-newID
              rm(newID,newClass,newData)
}else{ # If replicates exist, combine the 2 column ID into a single ID
    ID<- StructOut$ID
    if (min(dim(ID))==2){
        shortID<-as.double(ID[,1])+(as.double(ID[,2])*.001) # Now single entry is ID1.ID2
        StructOut$OldID<-StructOut$ID
        StructOut$ID<-shortID
        rm(ID,shortID)
                         }
    }
 
## Option 1: Normalize total ion count
# Apply the normalization factor to each row to normalize total ion count.
# We'll recalc norm factors in case data was replicate averaged.
if (StructOut$Options[1]){
    RowTotalIonCount<-rowSums(StructOut$Intensities)
    AvgTotalIonCount<-mean(RowTotalIonCount) #Population average
    NormFactor<-AvgTotalIonCount/RowTotalIonCount #Vector of norm factors
    StructOut$NormFactor<-NormFactor  #save this in the structure

m<-1
n<-cols
mx<-1
nx<- length(NormFactor)
NFmat<-matrix(t(matrix(NormFactor,mx,nx*n)),mx*m,nx*n,byrow=T)

 # match size of Intensities
    StructOut$Intensities<-StructOut$Intensities*as.vector(NFmat)
rm(NFmat,RowTotalIonCount,AvgTotalIonCount,NormFactor,m,n,mx,nx)
}
  ##  Option 5: Work with log (data)
 
if (StructOut$Options[5]){
    StructOut$Intensities<-log10(abs(StructOut$Intensities))}

return(StructOut)
}

