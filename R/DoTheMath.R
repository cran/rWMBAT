`DoTheMath` <-
function(InputStructure){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# DoTheMath takes a data set and performs feature selection 
# 
# DESCRIPTION
#       DoTheMath takes a data array, class vector, and other information
#       and builds and assesses a Bayesian network after selecting features
#       from within the data array.  It is called from the user interface
#       "orca.m." 
#
#       This is the umbrella script that loops a specified number of times
#       (see "repeats" below), each time doing a full n-fold cross
#       validation and recording the results.  All input and output data
#       are stored in a single data structure, described below.
# 
# USAGE
#       OutputDataStructure <- DoTheMath (InputStructure)
# 
# INPUTS
#       InputStructure: Data repository with fields: 
#       Intensities: Array of intensity values of size #cases x #variables
#       Class: Vector of length "#cases", with discrete values identifying
#          class of each case (may be integer)
#       ID: Patient ID array of length #cases, with one or more cols
#       MZ: Vector of length "#variables" holding labels for variables
#       Options: Logical 6x1 array. Options are:
#          1. Normalize on population total ion count (sum across rows)
#          2. Remove negative data values by setting them to zero
#          3. After normalizing, before binning, average cases with same ID
#          4. Find the MI threshold by randomization
#          5. Take log(data) prior to binning.  Negative values set to 1.
#          6. Remove Low Signal cases
#              NOT DONE: 3 Bin (2 Bin if False)
#       n: the "n" in n-fold cross validation
#       repeats: Times to repeat the whole process (e.g. re-crossvalidate)
#       threshold: Factor by which the maximum "random" MI us multiplied to
#           find the minimum "significant" MI (double, 1.0-5.0).
# 
# OUTPUTS
#       OutputDataStructure: all the fields of InputStructure, plus: 
#       ErrorRate: Vector containing misclassification rate for each repeat
#       KeyFeatures: Index to vector MZ that identifies features selected
# 
# CALLED FUNCTIONS
# 
#       InitialProcessing: Applies the options listed above
#       BuildBayesNet: Learns a Bayesian Network from the training data
#       ChooseMetaVars: Combines variables that may not be physically
#           separate molecules.
#       TestCases: Given the BayesNet, tests the "test group" to determine
#           the probability of being in each class.
#       opt3bin: Discretizes continuous data into 3 bins, optimizing MI
#       FindProbTables: Learns the values P(C,V) for each variable
#       cvpartition and training are MATLAB Statistics toolbox functions.
 
 
## Initialize

                
## Initial Processing
# According to options, remove negative values, normalize and/or take
# logarithm of data, replicate average. Store in output data structure.
 
print('Starting Initial Processing of Data')
OutputStructure <- InitialProcessing(InputStructure)
print("Initial processing complete.")
print ("")

# Get values out of Data structure to be used later 
drop<-InputStructure$Drop  # MI loss pecentage threshold for testing
                           # independance, see clipclassconnections
ff<-OutputStructure$Threshold
n<- OutputStructure$n # for n-fold cross validation default is 10
repeats<-OutputStructure$Repeats # Number of times to repeat CV, default 30
numtrials<-repeats*n
cverrorrate<-array(0,c(numtrials,1))
errorrate<-array(0,c(repeats,1))
data<-OutputStructure$Intensities
class<-OutputStructure$Class

# Find some sizes and initialize variables
rows<-nrow(as.matrix(data))
cols<-ncol(as.matrix(data))
class_predict<-matrix(0,rows,repeats)
class_prob<-matrix(0,rows,repeats)
trial<-0  # counter of how many times we perform Bayes Analysis (n*repeats)

ProbTables<-list()

## "Repeat Entire Process" Loop
 
# Repeat all processes the number of times requested
for (r in 1:repeats){

        print("")
        print(c("Working on repetition number",r,"at"))
        print(Sys.time())

        ## Cross Validation Loop
        # This section selects a training and testing group out of the data by
        # dividing it into n groups, and using n-1 of those for training and 1
        # for testing. 
        cvGroups <- array(0,dim=length(class))
        ### If class is multidimensional, then use cvGroups<-array(0,dim=length(class[,1]))
	gA <- sample(which(class==1))
        ### If class is multidimensional, then use gA <- sample(which(class[,1]==1))
	gB <- sample(which(class==2))
        ### If class is multidimensional, then use gB <- sample(which(class[,1]==1))
	test <- rep(seq(1,n),1+length(gA)/n)
	cvGroups[gA] <- test[1:length(gA)]
	test <- rep(seq(1,n),1+length(gB)/n)
	cvGroups[gB] <- test[1:length(gB)]
        ### cvGroups is now a vector that identifies the 'nfold' that each sample will be in the test group
        for (cv in 1:n) {# for each of n test groups, together spanning all cases
              
                trial<-trial+1 # Keep track of each trial
                print(c("Working on cross-validation number",cv,"of",n))
                # The test cases are identified by having a value of cv in cvGroups
                testgrpindex<-which(cvGroups == cv)
                testgrp <- data[which(cvGroups == cv),] 
                ### If class is multidimensional, then use testgrpclass <-class[which(cvGroups == cv),]
# add as.matrix on Aug_24 by Qian
                testgrpclass <- as.matrix(class[which(cvGroups == cv)])
                #testgrpclass <- class[which(cvGroups == cv)]
                # The training set is the rest
                traingrpindex <-which(cvGroups != cv)
                traingrp <- data[-which(cvGroups == cv),]
                traingrpclass <- class[-which(cvGroups == cv)]
                ### If class is multidimensional, then use traingrpclass <-class[-which(cvGroups == cv),]
                
                ## Discretize the groups into hi-med-low
                # by optimizing MI(V,C) for each V (feature) in the training data
                yarray<- opt3bin(traingrp,traingrpclass)
                leftbndry<-yarray$l
                rightbndry<-yarray$r
                traingrpbin<-yarray$binned
                maxMI<-yarray$mi
                rm(yarray)

                ## Build an augmented Naive Bayesian Network with the training data
                # The adjacency matrix is a logical with true values meaning "there
                # is an arc from row index to column index." The last row
                # represents the class variable.

                adjmat <- BuildBayesNet( traingrpbin, traingrpclass, ff, drop )
                n<- OutputStructure$n
                
                #store the  adjmat values for every trail 
                if(trial==1){
                       Adjacency<-array(NA,c(repeats*n,dim(as.matrix(adjmat))))}


                ## Find MetaVariables, rebuild data                    
                # Depending on the option set, reduce the V->V links by removing
                # them, or combining them into a single variable. The result is a
                # naive Bayesian network with only connections C->V

                meta_option<-1 # Hard coded for now
                classrow<-cols+1
                listvec<-1:cols # just a list of numbers
                               
                varlist<-unique(listvec[as.logical(adjmat[classrow,])]) # top level vars
                
                if (meta_option==1){
                   yarray<- ChooseMetaVars (traingrp, traingrpclass, adjmat)
                   finaldata<-as.matrix(yarray$finaldata)
                   metas<-yarray$metamatrix
                   leftbndry<-yarray$leftbound
                   rightbndry<-yarray$rightbound
save(yarray,file="cmv.rdat")
                   rm(yarray)      }

                   n<- OutputStructure$n

                   # store metas value for every trial
                   if(trial==1){
                         MetaVariablesFound<-array(NA,c(repeats*n,dim(as.matrix(metas))))}
                   
                   # Bin up the test group using these final results, combining
                   # variables per the instructions encoded in the "metas" logical
                   # matrix

                   testdata<-array(0,dim(testgrp))


                    if (length(varlist)==0){# in case no links are found
                                     print ("Not finding any links yet...")
                                     errorrate[trial]<- 1}
                    if (length(varlist)!=0){# if we do find links
                               for (var in  varlist){ # each of the parents of metavariables
                                             metavar<-c(var ,listvec[as.logical(metas[var,])])  # concatenate children
                                             metavar<-metavar[!is.na(metavar)]  
                                             # sum parent/child
                                             if (length(metavar)==1){ 
                                                        testdata[,var]<-testgrp[,metavar]}
                                             if(length(metavar)>1){ 
                                                        X<- testgrp[,metavar]
                                                              for(i in 1:dim(testdata)[1]){
                                                                     testdata[i,var]<-sum(X[i,])}
                                                        rm(X)} 
                                                   } 
                    # Now remove empty rows                  
                    finaltestdata<-testdata[,varlist]
                    

                    # And bin the result
                    testgrpbin<-array(0,dim(finaltestdata)) #will be stored here
                    
                     # Build boundary arrays to test against
                    testcases<-dim(testgrp)[1]


                    mx<-1
                    nx<-length(leftbndry)
                    m<-testcases
                    n<-1            
                    lb<- matrix(t(matrix(leftbndry,mx,nx*n)),mx*m,nx*n,byrow=T)
                    rm(mx,nx,m,n)
                    mx<-1
                    nx<-length(rightbndry)
                    m<-testcases
                    n<-1            
                    rb<- matrix(t(matrix(rightbndry,mx,nx*n)),mx*m,nx*n,byrow=T)
                    rm(mx,nx,m,n)
                    n<- OutputStructure$n

                    #  test each value and record the bin
                    testgrpbin[finaltestdata<lb]<-1
                    testgrpbin[finaltestdata>=lb]<-2
                    testgrpbin[finaltestdata>rb]<-3                       
                    
                    ## Populate Bayesian Network
            
                    # With the final set of data and the adjacency matrix, build the
                    # probability tables and test each of the test group cases, to see
                    # if we can determine the class.
            
                    # Build the probability tables empirically with the training group
                    # results
                    save(finaldata,traingrpclass,file="fpt.rdat")
                    ptable<-FindProbTables(finaldata, traingrpclass)
                    ProbTables[[trial]]<-array(NA,c(dim(ptable)))
                    Nclass<-sort(unique(traingrpclass))
                    lengthNclass<-length(Nclass)
                    prior<-array(0,c(lengthNclass,1))
                    Ntraingrpclass<- max(dim(as.matrix(traingrpclass)))
                           for(i in 1: lengthNclass){
                                prior[i]<-length(which(class== Nclass[i]))/ Ntraingrpclass }
                     rm(i, lengthNclass,Nclass, Ntraingrpclass)
                     
                     # find out the probability of each cases bing in class 1,2,etc.
                     # Cases are in rows, class in columns.
                     classprobtable <- TestCases (ptable, prior, testgrpbin)
                     rowClassprobtable<- dim(as.matrix(classprobtable))[1]
                     P_C<-array(NA,c(rowClassprobtable,1))
                     predclass<-array(NA,c(rowClassprobtable,1))
                     for(i in 1:rowClassprobtable){
                                     P_C[i]<-max(classprobtable[i,])
                                     predclass[i]<-which.max(classprobtable[i,])}
                     rm(rowClassprobtable)

                     class_prob[-traingrpindex,r]<-P_C
                     class_predict[-traingrpindex,r]<-predclass
                     
                     #Get the per trial error rate
                     cverrorrate[trial]<- sum(predclass==as.matrix(testgrpclass))/testcases
                     
                     #Store some "per trial" data
                     Adjacency[trial,,]<-adjmat
                     MetaVariablesFound[trial,1:cols,1:cols]<-metas
                     ProbTables[[trial]]<-ptable
        } # end of finding metavariables

    }# end of Cross Validation loop

    wrong<-sum(!(class==class_predict[,r]))
    errorrate[r]<-wrong/rows
    
} # of repeating entire process loop

MetaVariablesFound[is.na(MetaVariablesFound)]<-0

ErrorRate<-errorrate # one for each repeat
CvErrorRate<-cverrorrate # one for each of n*repeats trials
PredictedClass<-class_predict

# Find out the error for each case
mx<-length(class)
nx<-1
m<-1
n<-r
classrep<- matrix(t(matrix(class,mx,nx*n)),mx*m,nx*n,byrow=T)
rm(mx,nx,m,n)

# Record the results in the output structure
WasIright<-classrep==PredictedClass
CasePredictionRate<-rowSums(WasIright)/r
ClassProbability<-class_prob
SumAdj<- colSums(Adjacency)
MetaVars<-colSums(MetaVariablesFound)


OutputStructure<-list(Intensities=OutputStructure$Intensities,MZ=OutputStructure$MZ,ID=OutputStructure$ID,Class=OutputStructure$Class,Options=OutputStructure$Options,n=OutputStructure$n,Repeats=OutputStructure$Repeats,Drop=OutputStructure$Drop,Threshold=OutputStructure$Threshold,RawData=OutputStructure$RawData,NormFactor=OutputStructure$NormFactor,Replicate_ID=OutputStructure$ Replicate_ID,Replicate_Class=OutputStructure$Replicate_Class, Adjacency= Adjacency, MetaVariablesFound= MetaVariablesFound, ErrorRate= ErrorRate, CvErrorRate= CvErrorRate, PredictedClass= PredictedClass, CasePredictionRate= CasePredictionRate, ClassProbability= ClassProbability, ProbTables= ProbTables, SumAdj= SumAdj, MetaVars= MetaVars)

return(OutputStructure)

}

