\name{opt3bin}
\alias{opt3bin}
\title{Find 3 Bin Boundaries Optimizing the MI of Each Variable}
\description{
This function takes an array of continuous sample data of size
cases (rows) by variables (columns), along with a class vector of
integers, each integer specifying the class. The class vector 
has the same number of cases as the data.  The function outputs the
position of the 2 bin boundaries (3 bins) that optimize the mutual
information of each variable's data vector with the class vector.     
}
\usage{
opt3bin(data, class)
}
\arguments{
  \item{data}{double array of continuous values, cases in rows and  variables in columns. Distribution is unknown}
  \item{class}{double column vector, values 1:c representing classification  of each case}
}
\value{
  \item{l}{double row vector of left boundary position for each var}
  \item{r}{double row vector of right boundary position for each var}
  \item{binned}{double data array discretized using boundaries in l and r}
  \item{mi }{double row vector of mutual info between each discr. variable and class }
}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(traingrp,traingrpclass) #load example data in package       
result <- opt3bin(traingrp,traingrpclass)

}
