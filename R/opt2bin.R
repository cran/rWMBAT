`opt2bin` <-
function(rawdata, class, steps,typesearch, minint=NA, maxint=NA){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# opt2bin finds the best single boundary for each variable to maximize MI 
# 
# DESCRIPTION
#       This function takes an array of continuous data, with cases in rows
#       and variables in columns, along with a vector "class" which holds
#       the known class of each of the cases, and returns an array
#       "binneddata" that holds the 2 bin discretized data.  The
#       discretization bin boundary is found by maximizing the mutual
#       information with the class the resulting MI and boundary are also
#       returned. The starting boundaries for the search can be given in
#       the vectors min and max, or either one, or neither, in which case
#       the data values determine the search boundaries.# 
#
# USAGE
#       [mi boundary binneddata] = maxMIbin(rawdata, class, typesearch [,
#           min, max])
# 
# INPUTS
#       rawdata: double array of continuous values, cases in rows and 
#           variables in columns. Distribution is unknown.
#       class: double column vector, values 1:c representing classification
#           of each case. 
#       steps: Number of steps to test at while finding maximum MI
#       typesearch =0: starting bndry based on data's actual max/min values
#                  =1: use the value passed in max as maximum (right) value
#                  =-1: use the value passed in min as minimum (left) value
#                  =2: used values passed via max, min
#       the two optional arguments are vectors whose values limit the range
#       of search for each variables boundaries.
# 
# OUTPUTS
#
#       mi: row vector holding the maximum values of MI(CVi) found
#       boundary: The location used to bin the data to get max MI
#       binneddata: The resulting data binned into "1" (low) or "2" (hi)
# 
# CALLED FUNCTIONS
# 
#       MIarray: Finds the MI of each col in an array with a separate
#           vector (the class in this case)
 
## Intialize


rows<-nrow(as.matrix(rawdata))
cols<-ncol(as.matrix(rawdata))
mi<-array(0,cols)
boundary<-array(0,cols)
binneddata<-matrix(0,rows,cols)
currentmi<-matrix(0,steps,cols)
 
# if not passed, find the left and rightmost possible bin boundaries from data


if((!is.na(minint))&& (!is.na(maxint))&& typesearch!=2&&typesearch!=1&& typesearch!=1){
    print("typesearch must = 0,1,-1,2")
      }
if((!is.na(minint))&&(!is.na(maxint))&& typesearch==2){
                                                       print("using passed values")}
if((!is.na(minint))&&(!is.na(maxint))&&typesearch==-1){ minint<-array(NA,cols)
                                                         for(i in 1:cols){
                                                             maxint[i]<-max(as.matrix(rawdata)[,i])}}
if ((!is.na(minint))&&(!is.na(maxint))&&typesearch==1){ minint<-array(NA,cols)
                                                         for(i in 1:cols){
                                                             minint[i]<-min(as.matrix(rawdata)[,i])}}

if (is.na(minint)){minint<-array(NA,cols)
                       for(i in 1:cols){
                           minint[i]<-min(as.matrix(rawdata)[,i])}}


if(is.na(maxint)){maxint<-array(NA,cols)
                      for(i in 1:cols){
                          maxint[i]<-max(as.matrix(rawdata)[,i])}}
 
## Find best boundary
 
for (peak in 1:cols){ #look at each variable separately
 
    # Create an array of bin boundary's possible locations min->max
    
                           X<-seq(minint[peak],maxint[peak],length=steps)
                           mx<-1
                           nx<-length(X)
                           m<-rows
                           n<-1
checkpoints<- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
rm(X,mx,nx,m,n)

   # discretize the variable's values at each of these possible
    # boundaries, putting 2's everywhere (value > boundary), 1 elsewhere 
    #binarray=(repmat(rawdata[,peak], 1, steps)>checkpoints)+1 

mx<-rows
nx<-1
m<-1
n<-steps
X<- matrix(t(matrix(as.matrix(rawdata)[,peak],mx,nx*n)),mx*m,nx*n,byrow=T)
binarray<-(X>checkpoints)+1
rm(X,mx,nx,m,n)


    # find the MI(C,V) for each possible binning  
rbin<-nrow(binarray)
cbin<-ncol(binarray)
    mivec<-array(0,cbin)
    for (v in 1:cbin){
        mivec[v]<- MutualInfo(binarray[,v],class) # Fast using Pengs DLLS
                         }

    currentmi[1:steps,peak]<-mivec
    
    # Now pick out the highest MI, i.e. best bin boundary 
mi[peak]<-max(currentmi[,peak])
atstep<-which.max(currentmi[,peak])
    boundary[peak]<-checkpoints[1,atstep]
    
    # and record the binned data using that boundary.
    binneddata[,peak]<-binarray[,atstep]
}
return(list(mi=mi,boundary=boundary, binneddata=binneddata)) 
}

