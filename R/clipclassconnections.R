`clipclassconnections` <-
function( adj, mivc_vec, mivcv, dropthreshold ){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# clipclassconnections delinks variables from class 
# 
# DESCRIPTION
#        Where two variables are connected to each other and also 
#        to the class, attempt to select one as the child of the other and
#        disconnect it from the class. Use MI(Vi;C|Vj)<<MI(Vi;C) as a test.
# 
# USAGE
#       probtable = FindProbTables(data, class)
# 
# INPUTS
#       adj: (logical) matrix where "true" entries at (i,j) mean "an arc
#            exists from the Bayesian network node Vi to Vj." The class 
#            variable C is added at row (number of V's + 1). "0" values
#            mean no arc.
#       mivc_vec: (double) row vector containing MI(C;Vi) for each variable
#       mivcv: (double) array whose (i,j) entry is MI(Vi,C|Vj).
#       dropthreshold: percentage drop from MI(Vj;C) to MI(Vj;C|Vi) before
#           declaring that Vi is between C and Vj. 
# 
# OUTPUTS
#
#       adjout: copy of adj with the appropriate arcs removed.
# 
# CALLED FUNCTIONS
# 
#       None.


## Intialize


classrow<-nrow(adj)
numvars<-ncol(adj)


classconnect<-adj[classrow,] # the last row of adj stores arcs C->V
adjout<-matrix(NA,classrow, numvars)# placeholder for output array

## Identify triply connected arcs

# First look for pairs that are connected to each other and connected to
# the class. 

# Connected to each other: build logical array with (i,j) true if Vi<->Vj
vv_conn<-adj[1:numvars, 1:numvars]

# Connected to the class: logical array with (i,j) true if C->Vi and C->Vj

mx<-1
nx<-length(classconnect)
m<-numvars
n<-1

A<- matrix(t(matrix(classconnect,mx,nx)),mx*m,nx,byrow=T)
rm(mx,nx,m,n)

mx<-length(t(classconnect))
nx<-1
m<-1
n<-numvars

B<- matrix(t(matrix(t(classconnect),mx,nx*n)),mx*m,nx*n,byrow=T)
rm(mx,nx,m,n)

vcv_conn<-A&B 
rm(A,B)

# Find all (i,j) with both true
triple_conn<-vv_conn & vcv_conn

## Determine preferred direction on V<->V arcs

# Determine the Vi<->Vj direction by finding the greater of MI(C;i|j) or
# (C;j|i).  Greater MI means less effect of the instantiation of i or j.
arcdirection<-mivcv > t(mivcv) #Only the larger survive
dag_triple_conn<-arcdirection & triple_conn # Wipes out the smaller ->

# find links should NOT be kept under the test above,
linkstoremove<-(!arcdirection) & triple_conn
# and if they are in the connection list, remove them
adjout[1:numvars, 1:numvars]<-xor(vv_conn,linkstoremove)

# Now we need to test whether we can remove the link between C and which
# ever V (i or j) is the child of the other. We look for a "significant"
# drop in MI(Vj;C) when instantiating Vi, e.g. MI(Vj;C|Vi)<<MI(Vj;C).
#
# dropthreshold of .7, for example, means link breaks if 1st term is less
# than 30# of the second term.
#
# If there is a big drop in MI(C;Vj) when Vi is given, and Vi->Vj exists in
# the DAG, then we can remove the link C->Vj and leave C->Vi->Vj.

# Build an array out of the mivc_vec vectorX


mx<-nrow(t(mivc_vec))
nx<-ncol(t(mivc_vec))
m<-1
n<-numvars
mivc<- matrix(t(matrix(t(mivc_vec),mx,nx*n)),mx*m,nx*n,byrow=T)
rm(mx,nx,m,n)


# Test for the large drop described above
bigdrop<-((mivc-mivcv)/mivc) > as.numeric(dropthreshold)
# Test for the big drop and the V-V connection
breakconn<-t(bigdrop) & dag_triple_conn
# If any of the elements in a column of the result are true, remove that
# variable's C->V link, since it is a child.

linkstokeep<-matrix(NA,1,ncol(breakconn))
for(i in 1:nrow(breakconn)) {linkstokeep[i]<-!any(breakconn[,i])}


adjout[classrow,]<-adj[classrow,] & linkstokeep
return(adjout)
}

