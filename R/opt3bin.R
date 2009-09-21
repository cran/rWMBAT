`opt3bin` <-
function (data, class){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# FunctionName short description 
# 
# DESCRIPTION
#       This function takes an array of continuous sample data of size
#       cases (rows) by variables (columns), along with a class vector of
#       integers 1:c, each integer specifying the class. The class vector 
#       has the same number of cases as the data.  The function outputs the
#       position of the 2 bin boundaries (3 bins) that optimize the mutual
#       information of each variable's data vector with the class vector.    
# 
# USAGE
#       [l,r,binned, mi]=opt3bin(data,class)
# 
# INPUTS
#       data: double array of continuous values, cases in rows and 
#           variables in columns. Distribution is unknown.
#       class: double column vector, values 1:c representing classification
#           of each case. 
# 
# OUTPUTS
#
#       l     - row vector of left boundary position for each var.
#       r     - row vector of right boundary position for each var.
#       binned- data array discretized using boundaries in l and r
#       mi    - row vector of mutual info between each discr. variable 
#                  and class 
# 
# CALLED FUNCTIONS
# 
#       opt2bin: Similar function that finds a single boundary. This is
#           used as a seed for the 3 bin optimization.
#       looklr: See below.
 
 
## Intialize
# 
#  Variable Prep : find sizes of arrays and create placeholders for locals
 
steps<-150
rows<-nrow(as.matrix(data))
cols<-ncol(as.matrix(data))
boundary<-matrix(0,2,cols)
 
## Method
# Find starting point by finding the maximum value of a 2 bin mi. Next, go
# left and right from that position, finding the position of the
# next boundary that maximizes MI.

yarray<- opt2bin (data, class, steps, 2)
mi<-yarray$mi
boundary[1,]<-yarray$ boundary
rm(yarray)

#[mi boundary[1,]] = opt2bin (data, class, steps, 2)
 
# We've located a good starting (center) bin boundary.  Search L/R for a
# second boundary to do a 3 bin discretization.
#[mi boundary[2,]] = looklr (data, class, boundary[1,], steps)
#return(list(mi=mi,boundary=boundary, binneddata=binneddata))
yarray<-looklr(data, class, boundary[1,], steps)
mi<-yarray$miout
boundary[2,]<-yarray$nextboundary
rm(yarray)

# We've now found the optimum SECOND boundary position given the best 2 bin
# center boundary.  Now re-search using that SECOND boaundary position,
# dropping the original (2 bin).  The result should be at, or near, the
# optimal 3 bin position.
#[mi boundary[1,] binned] = looklr (data, class, boundary[2,], steps)
yarray<- looklr (data, class, boundary[2,], steps)
mi<- yarray$miout
boundary[1,]<- yarray$nextboundary
binned<- yarray$binned
rm(yarray)

# from the two boundaries found above, sort the left and right
r<-array(NA,cols)
l<-array(NA,cols)
for(i in 1:cols){
r[i]<-max(boundary[,i])
l[i]<-min(boundary[,i])} ##max min need test 
 
# Now retutn the vector of left and right boundaries, the disc. data, and
# max MI found.
return(list(l=l, r=r, binned=binned, mi=mi))
}

