\name{findmutualinfos}
\alias{findmutualinfos}
\title{Finds Various Mutual Info Combos among Variables}
\description{Given a set of data (many cases, each with values for many variables) and an additional value stored in the vector class, it finds MI described below in "Value."
}
\usage{
findmutualinfos(data, class)
}
\arguments{
  \item{data}{a number of cases (in rows), each with a measurement for a group of variables (in columns). The data should be discredited into integers 1 through k. The columns are considered variables V1, V2, ...}
  \item{class}{a column vector of length "cases" with integer values 1,2..., an additional measurement of class C}
}
\value{
  \item{mi_vc}{a row double vector whose ith value is MI(Vi,C)}
  \item{mi_vv}{double symmetric matrix with values MI(Vi,Vj)}
  \item{mi_vc_v}{double nonsym matrix with values MI(Vi,C;Vj)}
 }
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\note{CALLED FUNCTIONS
     MutualInfo(Calculate Mutual Information Of Two Variables) 
     CondMutualInfo(Calculate Mutual Information of Two Variables Conditioned On a Third)}
\examples{
data(traingrpbin, traingrpclass) # load example input data from package
result <- findmutualinfos(traingrpbin, traingrpclass)
         }
