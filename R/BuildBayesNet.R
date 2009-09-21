`BuildBayesNet` <-
function(data, class, ffactor, drop ) {
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# BuildBayesNet selects features and metafeatures based on mutual info.
#  
# 
# DESCRIPTION
#       This function takes a set of training data and an additional
#       variable called "class" and tries to learn a Bayesian Network
#       Structure by examining Mutual Information.  The class variable C is
#       assumed to be the ancestor of all other variables V.  Arcs from C
#       to V are declared if MI(CV)>>z, where z is a maximum expected MI
#       of similar, but random data...multiplied by a "fudge factor."  Arcs
#       from Vi to Vj are similarly declared. Then various tests are
#       performed to prune the network structure and combine variables that
#       exhibit high correlations. Finally the network is pruned to be a
#       Naive Bayesian Classifier, with only C->V arcs remaining.
# 
# USAGE
#       network_structure = BuildBayesNet( training_data, class )
######network_structure not use here  by Qian
# 
# INPUTS
#       training_data: cases in rows, variables in cols, integer array
#               containing the data used to build the Bayes net
#       class: the known class variable for each case (1:c col vector)
#       ffactor: multiple of auto MI to use to threshold C->V connections
#       drop: 
# 
# OUTPUTS
#
#       adjmatrix: a matrix of zeros and ones, where one in row i, column j
#               denotes a directed link in a Bayesian network between 
#               variable i and variable j. The class variable is the last
#               row/column.
# 
# CALLED FUNCTIONS
# 
#       automi: finds an MI threshold based on data
#       findmutualinfos: finds all values MI(VC), MI(VV) and MI(VC|V)
 
## Initialize
 
# Initialize the network object and some constants

automireps<-10 #times to repeat the auto MI thresholding to find avg.
 
# Check the sizes of various things
rows<-nrow(data)
cols<-ncol(data)
cases<-max(dim(as.matrix(class))) 
if(rows!=cases) print ("# of rows in the data and class must be equal.")
if (rows==cases) rm(cases)

 
# network.adjmat=zeros(cols+1) # all variables plus class as last row/col
dataalphabet<-length(unique(as.vector(data)))
# number of possible values of data
classalphabet<-length(unique(class)) # Number of values of class
 
## Step 0: Find all the necessary mutual information values, thresholds
# The function below finds all values MI(VC|V) and other combos needed and
# stores them in the network structure.
 

yarray<- findmutualinfos( data, class )
mi_vc<-yarray$mi_vc 
mi_vv<-yarray$mi_vv
mi_vc_v<-yarray$mi_vc_v
rm(yarray)



                            
# Find a threshold MI by examining MI under randomization
# ******************************
# Come back to the next line
# ****************************
#scalar MI threshold, 10 repetitions
vcthreshold<-automi(data, class, automireps)
vcthreshold<-vcthreshold*ffactor
vvthreshold <-vcthreshold * log(dataalphabet)/log(classalphabet)
 
 
## Step 1: Find all the possible arcs.
# Find the variables with high MI with the class, i.e. MI(V,C)>>0 and
# connect a link in the adjacency matrix C->V.  Also connect variable Vi,Vj
# if MI(ViVj)>>0
 

adjmat1<-getarcs(mi_vc, vcthreshold, mi_vv, vvthreshold)

## Step 2: Prune the variable set by clearing irrelevant features
# If there is no path from V to the class, clear all entries V<->Vi (all i) 
adjmat2 <-clearirrarcs(adjmat1)
 
## Step 3: Cut connections to class
# Where two variables are connected to each other and also to the class,
# attempt to select one as the child of the other amd disconnect it from
# the class. Use MI(ViC|Vj)<<MI(ViC) as a test.
 
temp <-clipclassconnections (adjmat2,t(mi_vc),mi_vc_v,drop)
 
# and once again clear features no longer near class and end function
adjacency<-clearirrarcs( temp )

 return(adjacency)
}

