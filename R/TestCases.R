`TestCases` <-
function ( p, prior, data){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# classprobs uses Bayes rule to classify a case
# 
# DESCRIPTION
#       Tests each of a set of data vectors by looking up P(data|class) in
#       a probability table, then finding P(case|class) by multiplying each
#       of those values in a product.  Then uses Bayes' rule to calculate
#       P(class|data) for each possible value of class.  Reports this as an
#       array of class probabilities for each case.
# 
# USAGE
#       classprobs = TestCases( p, prior, data)
# 
# INPUTS
#       data: double array of discrete integer (1:n) values, cases in rows 
#           and variables in columns.
#       p: 3-D double array of probabilities (c,d,v).  The first dimension 
#           is the class, the second is the data value, the third is the 
#           variable number. The entry is P(var v=value d | class=value c).
# 
# OUTPUTS
#
#       classprobs: 2-D double array whose value is P(class=c|data) for
#           each case. Cases are in rows, class in cols.
# 
# CALLED FUNCTIONS
# 
#       None.

## Intialize

# Find the sizes of the inputs and the number of possible values
data<-as.matrix(data)
cases<-nrow(data)
numvars<-ncol(data)
classvals<-dim(p)[1]
pvec<-matrix(0,classvals,numvars)
classprobs<-matrix(0,cases, classvals) # holds the classification results

## Find the probabilities

# Create pvec, an array whose first row is P(V=v|c=1) for each V
for (casenum in 1:cases){
    casedata<-data[casenum,]
    for (c in 1:classvals){
        for (v in 1:numvars){
            pvec[c,v]<-max(p[c,casedata[v],v],.01)
                            }
                          }
        Pdc<-matrix(NA,dim(pvec)[1],1) 
        for(i in 1:dim(pvec)[1]){Pdc[i]<-prod(pvec[i,]) }                           
    classprobs[casenum,]<-(Pdc*prior)/sum(Pdc*prior) 
                         }

return(classprobs)

}

