`CondMutualInfo` <-
function(V1,V2,condV){
###value in condition
N<-length(condV)
v3Val <- unique(condV)
n3 <- length(v3Val)

##for all the case in condv

Nxyz<-array(0,n3)
Nz<-array(0,n3)
Npart<-array(0,n3)

for(k in v3Val){
               test <- which(condV == v3Val[k])
               Nz[k]<-length(test)  
               Nxyz[k]<-MutualInfo(V1[test],V2[test])
               Npart[k]<-Nz[k]*Nxyz[k]/N              
               }
MIxyz<-sum(Npart)
return(MIxyz)               
}

