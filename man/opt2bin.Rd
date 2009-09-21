\name{opt2bin}
\alias{opt2bin}
\title{Finds The Best Single Boundary For Each Variable To Maximize MI}
\description{
This function takes an array of continuous data, with cases in rows
and variables in columns, along with a vector "class" which holds
the known class of each of the cases, and returns an array
"binneddata" that holds the 2 bin discretized data.
}
\usage{
opt2bin(rawdata, class, steps, typesearch, minint = NA, maxint = NA)
}
\arguments{
  \item{rawdata}{double array of continuous values, cases in rows and variables in columns. Distribution is unknown}
  \item{class}{double column vector, values 1:c representing classification of each case}
  \item{steps}{integer, number of steps to test at while finding maximum MI}
  \item{typesearch}{=0 starting bndry based on data's actual max/min values
                   =1  use the value passed in max as maximum (right) value
                   =-1 use the value passed in min as minimum (left) value
                   =2  used values passed via max, min }
  \item{minint}{vectors whose values limit the range of search for each variables boundaries}
  \item{maxint}{vectors whose values limit the range of search for each variables boundaries}
}
\details{
       The discretization bin boundary is found by maximizing the mutual
       information with the class the resulting MI and boundary are also
       returned. The starting boundaries for the search can be given in
       the vectors min and max, or either one, or neither, in which case
       the data values determine the search boundaries.
}
\value{
  \item{mi}{row vector holding the maximum values of MI(CVi) found}
  \item{boundary}{double vector, the location used to bin the data to get max MI}
  \item{binneddata}{resulting data binned into "1" (low) or "2" (hi)}
}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(traingrp,traingrpclass) #load example input data from package
result <- opt2bin(traingrp,traingrpclass,150,2)
         }
