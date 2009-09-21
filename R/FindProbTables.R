`FindProbTables` <-
function(data, class){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# FindProbTables estimates the probabilities P(class=c|data=D)  
# 
# DESCRIPTION
#       Input a training group of data arranged with cases in rows and 
#       variables in columns, as well as the class value c for that vector. 
#       Each case represents a data vector V.  For each possible data value 
#       vi, and each variable Vi, it calculates P(C=c|Vi=vi) and stores 
#       that result in a 3-D table.  The table is arranged with the 
#       dimensions (class value, data value, variable number).
# 
# USAGE
#       probtable = FindProbTables(data, class)
# 
# INPUTS
#       data: double array of discrete integer (1:n) values, cases in rows 
#           and variables in columns.
#       class: double column vector, also 1:n. Classification of each case.
# 
# OUTPUTS
#
#       probtable: 3-D array whose (c,d,v) value is P(class=c|data=p) for
#           variable v.
# 
# CALLED FUNCTIONS
# 
#       None.
 
## Intialize
# Find the sizes of the inputs and the number of possible values
cases<-nrow(data)
numvars<-ncol(data)
datavals<-length(unique(as.vector(data)))
classvals<-length(unique(as.vector(class)))
# Build some placeholders and loop indices
p<-array(0,dim=c(classvals, datavals, numvars )) # triplet: (class, value, variable#) 
databins<-1:datavals
classbins<-1:classvals
 
## Find Probabilities
# For each classification value, extract the data with that class
for (c in classbins){
save(data,class,c,file="r.rdat")
    datainthatclass<-as.matrix(data)[class==c,] # array of just cases with class=c
    datainthatclass<-as.matrix(datainthatclass)
# find the percentage of data with each possible data value
for( i in 1:numvars){
                     p[c,1,i]= sum(datainthatclass[,i]==1)/cases
                     p[c,2,i]= sum(datainthatclass[,i]==2)/cases
                     p[c,3,i]= sum(datainthatclass[,i]==3)/cases
                     }
                       }
return(p)
}

