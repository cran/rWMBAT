`findmutualinfos` <-
function ( data, class ){
# (c) Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.
#
# FunctionName short description 
# 
# DESCRIPTION
#       Input a training group of data arranged with cases in rows and 
# 
# USAGE
#       probtable = FindProbTables(data, class)
# 
# INPUTS
#       data: double array of discrete integer (1:n) values, cases in rows 
#           and variables in columns.
#       class: double column vector, also 1:n. Classification of each case.
# 
# OUTPUTS
#
#       probtable: 3-D array whose (c,d,v) value is P(class=c|data=p) for
#           variable v.
# 
# CALLED FUNCTIONS
# 
#       MIarray - finds the MI of each column in an array with a class vector
 
## Intialize
# Find the sizes of the inputs and the number of possible values# 
# THIS VERSION USES SOMEONE ELSE'S CODE, BUT ITS 10x FASTER THAN MINE.
# 
# FINDMUTUALINFOS finds the various mutual info combos among variables.
# Given a set of data (many cases, each with values for many variables) and
# an additional value stored in the vector class, it finds MI described
# below in "OUTPUTS."
 
# 
# INPUTS:
# 
# data: A number of cases (in rows), each with a measurement for a group of
#    variables (in columns). The data should be discretized into integers 1
#    through k. The columns are considered variables V1, V2, ...
# class: an additional measurement of class C. A column vector of length 
#    "cases" with integer values 1,2...
# 
# OUTPUTS (all type double, >0)
# 
# mi_vc: a row vector whose ith value is MI(Vi,C).
# mi_vv: Symmetric matrix with values MI(Vi,Vj).
# mi_vc_v: Non-sym matrix with values MI(ViC|Vj).
# 
# CALLED FUNCTIONS
# 
# mutualinfo and condmutualinfo are from the mutualinfo package (c) 2002 by
# Hanchuan Peng, <penghanchuan@yahoo.com>.
 
# Calculate the value MI(Vi,C)

rows<-nrow(data)
cols<-ncol(data)
mi_vc<-array(0,cols)
 
for (v in 1:cols){
mi_vc[v]<- MutualInfo(data[,v],class) } #using MutualInfo function

# For each variable Vj, calculate MI(Vi,Vj) and MI(ViC|Vj)
 
        mi_vv<-matrix(0,cols,cols)
        mi_vc_v<-matrix(0,cols,cols)
        
for (i in 1:cols){
    for (j in 1:cols){
        mi_vv[i,j]<-MutualInfo (data[,i],data[,j])
        mi_vc_v[i,j]<- CondMutualInfo (data[,i],class,data[,j])}
                     }
 return(list(mi_vc= mi_vc, mi_vv= mi_vv, mi_vc_v= mi_vc_v))

}

