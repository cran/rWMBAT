\name{looklr}
\alias{looklr}
\title{Finds Boundary}
\description{
 Given a start position, finds another boundary (to create 3 bins) that
 maximizes MI with the class.
}
\usage{
looklr(data, class, startbd, steps)
}
\arguments{
  \item{data}{double array, cases in rows and variables in columns}
  \item{class}{double column vector, values 1:c representing classification of each case  }
  \item{startbd}{double vector, given start position for each case }
  \item{steps}{integer, number of steps to test at while finding maximum MI}
}
\value{
  \item{miout }{double vector, recorded highest MI value}
  \item{nextboundary }{double vector, boundary (to create 3 bins) that maximizes MI with the class}
  \item{binned }{double matrix, recorded the binned value}
}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(traingrp,traingrpclass,boundary)#load the example data input data from package
result<- looklr (traingrp,traingrpclass,boundary, 150)
         }
