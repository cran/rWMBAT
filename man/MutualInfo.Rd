\name{MutualInfo}
\alias{MutualInfo}
\title{Calculates Mutual Information Of The Two Variables}
\description{Calculates mutual information of the two variables}
\usage{
MutualInfo(v1, v2)
}
\arguments{
  \item{v1}{one of the two vectors }
  \item{v2}{one of the two vectors }
}
\value{mutual information of the two vectors}
\author{Bill Cooke, College of William and Mary, Dept. of Physics, 2009 }
\examples{
a<-c(1 ,2, 1, 2, 1)
b<-c(2, 1, 2, 1, 1)
MI <- MutualInfo(a, b)
}
