\name{getarcs}
\alias{getarcs}
\title{Builds The Adjacency Matrix For A Set Of Variables}
\description{
By comparing mutual information between two variables to thresholds
determined seperately, this function declares there to be an arc in
a Bayesian network. Arcs are stored in an adjacency matrix,
described below.
}
\usage{
getarcs(mvc, vcthreshold, mvv, vvthreshold)
}
\arguments{
  \item{mvc}{double vector array with mutual information between variables and the class (variables and other variables). The (i,j) entries of mvv are MI(Vi,Vj)}
  \item{vcthreshold}{scalar threshold used to test for existence linkz }
  \item{mvv}{double vector matrix}
  \item{vvthreshold}{scalar threshold used to test for existence linkz}
}
\details{
  The primary tests are
  MI[Vi,Cj]>>vcthreshold (tests for links between Vi and the class)
  MI[Vi,Vj]>>vvthreshold (tests the links between variables)
}
\value{logical matrix whose entries "1" at (i,j) mean "an arc exists from the Bayesian network node Vi to Vj." The class variable C is added at row (number of V's + 1). "0" values mean no arc}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(mi_vc, vcthreshold, mi_vv, vvthreshold) # load the example input data from the package
adjmat<-getarcs(mi_vc, vcthreshold, mi_vv, vvthreshold)
         }
