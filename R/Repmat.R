`Repmat` <-
function(X,m,n){
if((dim(as.matrix(X))[1]==1)||(dim(as.matrix(X))[2]==1)) {X<-t(X)}
mx<-dim(as.matrix(X))[1]
nx<-dim(as.matrix(X))[2]
repmat<-matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
return(repmat)}

