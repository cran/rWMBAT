`automi` <-
function ( data, class, repeats ){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# automi finds a threshold for randomized MI(V C) 
# 
# DESCRIPTION
#       Finds the threshold of a data set's mutual information with a class
#       vector, above which a variable's MI(class, variable) can be
#       expected  to be significant. The threshold for mi (significance 
#       level) is found by taking the data set and randoomizing the class
#       vector, then calculating MI(CV) for all the variables. This is
#       repeated a number of times. The resulting list of length (#repeats
#       * #variables) is sorted,  and the 99th percentile max MI is taken
#       as the threshold.
 
# USAGE
#       threshold = automi( data, class )
# 
# INPUTS
#       data: double array of discrete integer (1:n) values, cases in rows 
#           and variables in columns.
#       class: double column vector, also 1:n. Classification of each case.
#       repeats: the number of times to repeat the randomization
# 
# OUTPUTS
#
#       threshold: the significance level for MI(CV)
# 
# CALLED FUNCTIONS
# 
#       MIarray(data,class): returns a vector with MI(ViClass) for each V
#           in the data set
 
## Intialize
 
# Find the size of the data (cases x variables) and check against class

rows<-nrow(data)
cols<-ncol(data)
cases<-max(length(class))

if (rows!=cases){
                 print("# of rows in the data and class must be equal.")
                }

if (rows==cases){
                rm(cases)
                }
    
## Repeat a number of times
 
mifound<-matrix(0,cols,repeats) # stores the results of the randomized MI
for (i in 1:repeats){   
                     c<-class[sample(rows)]  # creates a randomized class vector
                     rbin<-nrow(data)
                     cbin<-ncol(data)
                     mivec<-array(0,cbin)
                     for (v in 1:cbin){
                                       mivec[v]<-MutualInfo(data[,v],c) 
                                      }
                     mifound[,i]<-mivec
                     }
 
# pull off the 99th percentile highest MI
mi_in_a_vector<-matrix(mifound,nrow=repeats*cols) 
threshold<-quantile(mi_in_a_vector,prob=.99) 
return(threshold)
}

