`ChooseMetaVars` <-
function( data, class, adj){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# ChooseMetaVars attempts to combine variables into better variables 
# 
# DESCRIPTION
#       Finds the V-V pairs in the adjacency matrix, and attempts
#       to combine them into a metavariable with a higher mutual
#       information than either variable alone. If it is possible to do
#       this, it returns a new data matrix with the variables combined. 
# 
# USAGE
#       [finaldata metamatrix leftbound rightbound] =
#                        ChooseMetaVars ( data, class, adj)
# 
# INPUTS
#       data: double array of discrete integer (1:n) values, cases in rows 
#           and variables in columns.
#       class: double column vector, also 1:n. Classification of each case.
#       adj: Adjacency matrix, #variables+1 by #variables. Last row is
#           class node. Logical meaning "there is an arc from i to j."
# 
# OUTPUTS
#       metamatrix: logical whose (i,j) means "variable j was combined into
#           variable i (and erased)"
#       finaldata: The data matrix with the variable combined and rebinned
#       leftbound: The new left boundary (vector) for binning.
#       rightbound: The new right boundary (vector) for binning.    
# 
# CALLED FUNCTIONS
#       opt3bin: rebins combined variables to determine highest MI.
 
## Intialize

rows<-nrow(data)
cols<-ncol(data)
classrow<-nrow(adj)
numvars<-ncol(adj)
bindata<-matrix(0,rows,cols)
metamatrix<-matrix(NA,cols,cols)
 
# Create a list of all the variables V to check by examining the adjacency
# matrix's last row, i.e. those with C->V connections
listvec=1:numvars
varstocheck<-unique(listvec[as.logical(adj[classrow,])])
l<-matrix(0,1,numvars)
r<-matrix(0,1,numvars)
 
# Now go through that list, testing each V->W connection to see if adding V
# and W creates a new variable Z that has a higher MI with the class than V
# alone.  V is the list above, W is the list of variables connected to a V.
 
for (v in varstocheck){ # Pull out the W variables connected to V and test
wlist<-unique(listvec[as.logical(adj[v,])])
yarray<- opt3bin(data[,v], class)
l[v]<-yarray$l
r[v]<-yarray$r
binned<-yarray$binned
mitobeat<-yarray$mi
rm(yarray)
    #[l(v), r(v), binned, mitobeat] = opt3bin(data(:,v), class)
    bindata[,v]<-binned
    if (length(wlist)!=0){
        for (w in wlist){
            newdata<-data[,v]+data[,w]
            #[left, right, binned, newmi] = opt3bin(newdata, class)
            yarray<- opt3bin(newdata, class)
            left<-yarray$l
            right<-yarray$r
            binned<-yarray$binned
            newmi<-yarray$mi
            rm(yarray)
            if (newmi>mitobeat){
                mitobeat<-newmi
                data[,v]<-data[,v]+data[,w]
                metamatrix[v,w]<-TRUE # record the combination
                bindata[,v]<-binned 
                l[v]<-left
                r[v]<-right}
                       }
      }     
    
                
                # ********************************************************
                #              Too simple - should check all combos
                #    May be that a later variable is better
                # ***************************************************
}
 
#pull out just the V->C columns from the data matrix.
finaldata<-bindata[,as.logical(adj[classrow,])]
leftbound<-l[as.logical(adj[classrow,])]
rightbound<-r[as.logical(adj[classrow,])]

return(list(finaldata= finaldata, metamatrix= metamatrix, leftbound= leftbound, rightbound= rightbound))
}

