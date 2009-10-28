\name{CondMutualInfo}
\alias{CondMutualInfo}
\title{Calculates The Mutual Information Of Two Variables Conditioned On A Third }
\description{Calculates the mutual information of two  variables conditioned on a third}
\usage{
CondMutualInfo(V1, V2, condV)
}
\arguments{
  \item{V1}{ one of the two vectors }
  \item{V2}{ one of the two vectors }
  \item{condV}{ given condition}
}
\value{mutual information of two variables conditioned on a third}
\author{Qian Si, College of William and Mary, Dept. of Physics, 2009}
\examples{
a<-c(1 ,2, 1, 2, 1)
b<-c(2, 1, 2, 1, 1)
condV<-c(2, 1, 2, 2, 1)
MIxyz <- CondMutualInfo(a, b, condV)
}
