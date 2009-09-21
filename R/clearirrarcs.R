`clearirrarcs` <-
function(adjin){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# CLEARIRRARCS clears arcs that are not C->V or C->V<->V 
# 
# DESCRIPTION
#       Given an adjacency matrix with V<->V arcs in a square matrix and an
#       additional row representing C->V (class to variable), this function
#       clears out all V1->V2 arcs where V1 is not a member of the set of
#       V's that are class-connected, i.e. have arcs in the final row.
# 
# USAGE
#       adjout = clearirrarcs( adjin )
# 
# INPUTS
#       adjin: a logical array where a true value at position (i,j) means
#           that there is an arc in a directed acyclic graph between
#           (variable) i and variable j.
# 
# OUTPUTS
#       adjout: copy of adjin with unneeded arcs cleared
# 
# CALLED FUNCTIONS
#       None.

## Intialize
# Find the sizes of the input
classrow<-nrow(adjin)
numvars<-ncol(adjin)

## Main processing
# Find out which variables are connected to class
conntocls<-adjin[classrow,]

# Remove all arcs that don't have at least one variable in this list,
# e.g. all Vi<->Vj such that ~(Vi->C or Vj->C). These are all the entries 
# in the adjacency matrix whose i and j are NOT in the list above.

# Make a matrix with ones where neither variable is in the list above
A<-!conntocls
mx<- 1                 
nx<- length(A)
B<-matrix(t(matrix(A,mx,nx)),mx*numvars,nx,byrow=T)
A<-t(!conntocls)
mx<-length(A)
nx<-1
C<-matrix(t(matrix(A,mx,nx*numvars)),mx,nx*numvars,byrow=T)
noconnmat<-B&C

# Use that to erase all the irrelevant entries in the square adj matrix, at
# the same time remove the diagonal (arcs Vi<->Vi)
adjout<-adjin[1:numvars,1:numvars]&!noconnmat&!diag(numvars)

# Bidirectional arcs are temporarily permitted between nodes connected
# directly to the class, but not between nodes where only one is connected
# to the class- those are assumed to flow C->V1->V2 only.  Remove V2->V1.

# Get a matrix of ones in rows that are class connected. V->V arcs are only
# allowed to be in these rows:
A<-t(conntocls)
mx<-length(A)
nx<-1
parents<-matrix(t(matrix(A,mx,nx*numvars)),mx,nx*numvars,byrow=T)
# Remove anything else
adjout<-adjout & parents

# Now add back in the class row at the bottom of the square matrix

ADJOUT<-adjin
ADJOUT[1:(classrow-1),]<-adjout
adjout<-ADJOUT

return(adjout)
}

