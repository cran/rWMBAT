`MutualInfo` <-
function(v1, v2) {
 MI <-0
 test1 <- miCalc(v1,v2,length(v1),MI)
 return(test1[[4]])
}

