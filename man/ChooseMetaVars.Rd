\name{ChooseMetaVars}
\alias{ChooseMetaVars}
\title{ChooseMetaVars Attempts To Combine Variables into Better Variables}
\description{
  Finds the V-V pairs in the adjacency matrix, and attempts
  to combine them into a metavariable with a higher mutual
  information than either variable alone. If it is possible to do
  this, it returns a new data matrix with the variables combined. 
}
\usage{
ChooseMetaVars(data, class, adj)
}
\arguments{
  \item{data}{double array of discrete integer values, cases in rows and variables in columns.}
  \item{class}{double column vector. Classification of each case}
  \item{adj}{logical adjacency matrix, variables+1 by variables. Last row is class node. Logical meaning "there is an arc from i to j.}
}
\value{
  \item{metamatrix}{logical whose (i,j) means "variable j was combined into variable i (and erased)"}
  \item{finaldata}{double data matrix with the variable combined and rebinned}
  \item{leftbound}{vector,the new left boundary (vector) for binning}
  \item{rightbound}{vector,the new right boundary (vector) for binning}
}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\note{CALLED FUNCTIONS
      opt3bin(rebins combined variables to determine highest MI)}
\examples{
data(traingrp, traingrpclass, adjmat) # load example input data in the package
result <- ChooseMetaVars (traingrp, traingrpclass, adjmat)
         }
