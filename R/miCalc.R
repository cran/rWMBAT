`miCalc` <-
function(v1,v2,n,MI) .C("calcMI", as.integer(v1), as.integer(v2),as.integer(n), MI)

