\name{Repmat}
\alias{Repmat}
\title{Creates A Large Matrix Consisting Of An m-by-n Tiling Of Copies Of x}
\description{
 Creates a large matrix  consisting of an m-by-n tiling of copies of X
}
\usage{
Repmat(X, m, n)
}
\arguments{
  \item{X}{ matrix used to make copies of }
  \item{m}{ row copy number of the result matrix}
  \item{n}{ colomn copy number of the result matrix}
}
\value{matrix, m-by-n tiling of copies of input matrix }
\author{Qian Si,College of William and Mary, Dept. of Physics, 2009 }
\examples{
Repmat(diag(3),2,3)
# diag(3) creates a 3-by-3 indentity matrix
}

