`getarcs` <-
function( mvc, vcthreshold, mvv, vvthreshold ){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# GETARCS builds the adjacency matrix for a set of variables 
# 
# DESCRIPTION
#       By comparing mutual information between two variables to thresholds
#       determined seperately, this function declares there to be an arc in
#       a Bayesian network. Arcs are stored in an adjacency matrix,
#       described below.
# 
#       The primary tests are:
#       MI(Vi;Cj)>>vcthreshold : tests for links between Vi and the class
#       MI(Vi;Vj)>>vvthreshold : tests the links between variables
# 
# USAGE
#       adjacency = getarcs( mvc, vcthreshold, mvv, vvthreshold )
# 
# INPUTS
#       mvc [mvv]: double vector [array] with mutual information between
#           variables and the class [variables and other variables]. The
#           (i,j) entries of mvv are MI(Vi,Vj).
#       vc/vvthreshold: scalar threshold used to test for existence linkz
# 
# OUTPUTS
#
#       adjacency: logical matrix whose entries "1" at (i,j) mean "an arc
#            exists from the Bayesian network node Vi to Vj." The class 
#            variable C is added at row (number of V's + 1). "0" values
#            mean no arc.
# 
# CALLED FUNCTIONS
# 
#       None.
# 
# For more information on the tests and the links, see my dissertation.


## Initialize

numvars<-max(dim(as.matrix(mvc))) #the number of variables
classrow<-numvars+1; #row to store links C->V
adjacency<-matrix(NA,classrow,numvars) #the blank adjacency matrix

## Test for adjacency to class
adjacency[classrow,]<- mvc > as.numeric(vcthreshold)


## Test for links between variables
# This test results in a symmetric logical matrix since MI (X;Y) is
# symetric. To create a directed graph, these arcs will need to be pruned.
adjacency[1:numvars, 1:numvars]<- mvv > as.numeric(vvthreshold)
return(adjacency)
}

