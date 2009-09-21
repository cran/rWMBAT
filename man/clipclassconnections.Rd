\name{clipclassconnections}
\alias{clipclassconnections}
\title{Delinks Variables From Class}
\description{
 Where two variables are connected to each other and also 
 to the class, attempt to select one as the child of the other and
 disconnect it from the class. Use MI(Vi;C|Vj)<<MI(Vi;C) as a test. 
}
\usage{
clipclassconnections(adj, mivc_vec, mivcv, dropthreshold)
}
\arguments{
  \item{adj}{logical matrix where "true" entries at (i,j) mean "an arc exists from the Bayesian network node Vi to Vj." The class variable C is added at row (number of V's + 1). "0" values mean no arc.}
  \item{mivc_vec}{double row vector containing MI(C;Vi) for each variable}
  \item{mivcv}{double array whose (i,j) entry is MI(Vi,C|Vj)}
  \item{dropthreshold}{double percentage drop from MI(Vj;C) to MI(Vj;C|Vi) before declaring that Vi is between C and Vj}
}
\value{copy of adj with the appropriate arcs removed}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(adjmat2,mi_vc,mi_vc_v,drop) #load input example data from package
adjout <- clipclassconnections(adjmat2,t(mi_vc),mi_vc_v,drop)
         }
