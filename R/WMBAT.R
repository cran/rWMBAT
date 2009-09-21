`WMBAT` <-
function (tofListMetaData , alignedPeakList ,Options, nfold, repeats, threshold){
# The William and Mary Bayesian Analysis Tool
#           (c) 2009 Karl Kuschner, Qian Si and William Cooke, College of William and Mary
# 
# DESCRIPTION
#       WMBAT takes an array of mass spec peak intensities, a vector
#       describing which of two classes each sample belongs to, 
#       and other information and builds and assesses a Bayesian network
#       after selecting features (peaks) from within the data array that
#       are diagnostic of the class. The primary output is an adjacency
#       matrix describing the resulting Bayesian network.
# 
# USAGE
#       [IntOut IDOut PredClass SumAdj SumMV TrialErr] = WMBAT (Intensities,
#                       Class, ID, MZ, Options, nfold, repeats, threshold, drop)
# 
# INPUTS
#       Intensities: Double array of intensity values of size #cases x #variables.
#           Each row is a case (spectrum), each column a specific global 
#           m/z (mass) position. Each entry in the row is the intensity
#           value of that case's spectrum at that specific m/z position.
#       Class: Integer vector of length "#cases", values 1 or 2 identifying the 
#           class of each case, such as "disease, non-disease"
#       ID: Double one or two column array of length #cases containg the sample ID
#           for each case. Second column is optional and would identify
#           replicates of the same sample. 
#       MZ: Double vector of length "#variables" holding m/z labels for peaks
#       Options: Logical 6x1 array. Options are:
#          1. Normalize on population total ion count (sum across rows)
#          2. Remove negative data values by setting them to zero
#          3. After normalizing, before binning, average cases with same ID
#          4. NOT USED - SET TO FALSE
#          5. Take log(data) prior to binning.  Negative values set to 1.
#          6. NOT USED - SET TO FALSE
#       nfold: the "n" in n-fold cross validation (integer 4-10). 10 is
#           recommended.
#       repeats: Integer, times to repeat the whole process (e.g.
#           re-crossvalidate). 100 is recommended.
#       threshold: Factor by which the maximum "random" MI is multiplied to
#           find the minimum "significant" MI (double, 1.0-5.0). We
#           recommend starting with 1 and increasing until a "reasonable"
#           number of diagnostic peaks is reached and error rates are
#           minimized. This setting is dependant on the data and the
#           correlations between variables.
#       drop: MI loss pecentage threshold for testing independance. Set to
#           .75 (75#) and adjust to filter too few/too many variable-to-
#           variable connections.
# 
# OUTPUTS
#       IntOut: The Intensities input array, after processing by the
#           various options selected by the logical Options above.
#       IDOut: The ID number of each row in the IntOut array. With no replicate
#           averaging, each ID will be preserved (but reformatted) from the
#           input.  With replicate averaging, only the primary ID number
#           remains.
#       PredClass: The predicted class of each case, during each of the
#           trials (from input "repeats") 
#       Class2Vars: A vector whose ith value is the fraction of times peak i
#           (from the vector MZ) was selected as being connected to the
#           class. The maximum times it could have been selected was
#           nfold*repeats.
#       Var2Vars: An integer array whose (i,j) entry is the fraction of times a
#           second level link was found from peak i to peak j, when peak i
#           was connected to the class, as found in SumLvl1. 
#        MetaVars: An integer array whose (i,j) entry is the fraction of times a
#           metavariable was created using peak i and peak j and stored in
#           the level 1 variable peak i, once peak i was found connected to
#           the class.
#       TrialErr: The error rate for each of the "repeats" possible
#           trials. Records the percentage of cases where PredClass was not
#           equal to the input Class.
# 
# CALLED FUNCTIONS
# 
#       DoTheMath: Learns a Bayesian Network from the data
 
## Initialize
# Package the inputs into a data structure, as needed by DoTheMath
 
### Convert to m/z
### U*((_A+Kscale*_A*(_K-1))*(C2*Tdelta-_T0)^2+_B)
Tdelta <- 4.00E-09
Kscale <- 1
U <- 20000
A <- 257380326.6
B <- 0.000121327
T0 <- 1.66949E-07
K <- 0.999927616
MZ <- U*((A + Kscale*A*(K-1))*(alignedPeakList$peaks*Tdelta-T0)^2+B)

spectraName <- names(alignedPeakList$data)
###Intensities,ID,Class
Intensities<- array(0,dim=c(length(alignedPeakList$data),length(alignedPeakList$peaks)))
ID<-array(0,dim=c(length(alignedPeakList$data),2))
Class<-array(0,length(alignedPeakList$data))
for (k in 1:length(alignedPeakList$data)) {
                             Intensities[k,] <-as.double(alignedPeakList$data[[spectraName[k]]]$Intensities)
                             if(tofListMetaData[k,"sampleInfo.groupName"]=="Normal") Class[k]<-1
                             if(tofListMetaData[k,"sampleInfo.groupName"]== "Leukemia")  Class[k]<-2      
                             ID[k,1]<-tofListMetaData[k,"sampleInfo.sampleName"]
                             ID[k,2]<-tofListMetaData[k,"replicateNumber"]
                             }

In<-list(Intensities=Intensities,MZ=MZ,ID=ID,Class=Class,Options=Options,nfold=nfold,Repeats=repeats,Drop=0.75,Threshold=threshold)
 
## Call the main function
 
Out<-DoTheMath(In)
 
## Format the results
numvars<-max(dim(as.matrix(MZ)))
numtrials<-nfold*repeats
IntOut <- Out$Intensities
IDOut <- Out$ID
PredClass <- Out$PredictedClass
FirstLevel<-Out$SumAdj[numvars+1,]
Class2Vars<-FirstLevel/numtrials
 
Var2Vars<-matrix(0,numvars,numvars)
MetaVars<-matrix(0,numvars,numvars)
 
Chances<-Repmat(t(FirstLevel),1, numvars) # How often a variable could be linked
VtoV<-Out$SumAdj[1:numvars,1:numvars] # Variable to variable connections
#t<-matrix(1,dim(VtoV)[1],dim(VtoV)[2])
#for(i in 1:dim(VtoV)[1]){ index<-which(VtoV[i,]!=0)
#                          if(length(index)!=0) t[i,index]=0
#                         } # Places where there are no V to V conn
#Var2Vars[-t]<-VtoV[-t]/Chances[-t] # scale the V to V
for(i in 1:dim(VtoV)[1]){ index<-which(VtoV[i,]!=0)
                          if(length(index)!=0) {
                                               if(Chances[i,1]!=0) {                                                                                              
                                               Var2Vars[i,index]<-VtoV[i,index]/Chances[i,index]}
                                               }
                         } 
MetaVars[Out$MetaVars!=0]<-Out$MetaVars[Out$MetaVars!=0]/Chances[Out$MetaVars!=0]
TrialErr<-Out$ErrorRate
 

return (list(IntOut= IntOut, IDOut=IDOut, PredClass= PredClass, Class2Vars= Class2Vars, Var2Vars= Var2Vars, MetaVars= MetaVars, TrialErr= TrialErr))

}

