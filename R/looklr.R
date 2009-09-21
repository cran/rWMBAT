`looklr` <-
function (data, class, startbd, steps){
# given a start position, finds another boundary (to create 3 bins) that
# maximizes MI with the class

rows<-nrow(as.matrix(data))
cols<-ncol(as.matrix(data))
farleft<-array(NA,cols)
farright<-array(NA,cols)
for(i in 1:cols){
                 farleft[i]<-min(as.matrix(data)[,i])
                 farright[i]<-max(as.matrix(data)[,i])
                 }
miout<-array(0,cols)
binned<-matrix(0,rows,cols)
nextboundary<-array(0,cols)
 
for (peak in 1:cols){ 
    # for each peak/variable separately...
    # discretize this variables' values. Sweep through the possible
    # bin boundaries from the startbd to the furthest value of the
    # data, creating 2 boundaries for 3 bins. Record the binned values in
    # a "cases x steps" array, where "steps" is the granularity of the
    # sweep. The data vector starts off as a column...

X<-as.matrix(data)[,peak]
mx<-rows
nx<-1
m<-1
n<-steps
testmat<- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
rm(X,mx,nx,m,n)


# Create same size array of bin boundaries. Each row is the same.
X<-seq(farleft[peak],startbd[peak],length=steps)
mx<-1
nx<- length(X)
m<-rows
n<-1
checkptsL<- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
rm(X,mx,nx,m,n)
X<- seq(startbd[peak],farright[peak],length=steps)
mx<-1
nx<-length(X)
n<-1
m<-rows
checkptsR<- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
rm(X,mx,nx,m,n)




    # Create a place to hold the discrete info, starting with all ones. The
    # "left" array will represent data binned holding the center boundary
    # fixed and sweeping out a second boundary to the left similarly the
    # right boundary starts at "startbd" and sweeps higher.
    binarrayL<-matrix(1,rows,steps) 
    binarrayR<-matrix(1,rows,steps) 
   
    # Those in the L test array that are higher than the left boundary -> 2
    binarrayL[testmat>checkptsL]<-2
    binarrayL[testmat>startbd[peak]]<-3 # >center boundary -> 3
    
    # Similarly using center and right boundaries
    binarrayR[testmat>startbd[peak]]<-2
    binarrayR[testmat>checkptsR]<-3
 
    # Now at each of those step positions, check MI (varclass).
    miout[peak]<- 0
    for (j in 1:steps){
        miL <- MutualInfo(binarrayL[,j],class)# MI(VC) using left/center
        miR <- MutualInfo(binarrayR[,j],class)# MI(VC) using center/right
        if (miL>miout[peak]){ # check if that steps MI is the highest yet
            miout[peak]<-miL # if so, record it
            newboundary<-checkptsL[1,j] # and record the boundary
            binned[,peak]<-binarrayL[,j]} # and record the discrete data
        
        if (miR>miout[peak]){ # and check the center/right combo similarly
            miout[peak]<-miR
            newboundary<-checkptsR[1,j]
            binned[,peak]<-binarrayR[,j]}
            
    } # checking each possible boundary position.
    
    # we should now know the best possible place to put a second boundary,
    # either left or right of a given starting position. Record that.
    nextboundary[peak]<-newboundary
    
}# of that variable's search.  Go to next variable.
 
return(list(miout= miout, nextboundary= nextboundary, binned= binned))
 }

