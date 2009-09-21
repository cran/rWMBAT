\name{clearirrarcs}
\alias{clearirrarcs}
\title{Clears Arcs That Are Not C->V Or C->V<->V }
\description{
  Given an adjacency matrix with V<->V arcs in a square matrix and an
  additional row representing C->V (class to variable), this function
  clears out all V1->V2 arcs where V1 is not a member of the set of
  V's that are class-connected, i.e. have arcs in the final row.
}
\usage{
clearirrarcs(adjin)
}
\arguments{
  \item{adjin}{logical array where a true value at position (i,j) means that there is an arc in a directed acyclic graph between (variable) i and variable j. }
}
\value{copy of adjin with unneeded arcs cleared} 
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke,College of William and Mary, Dept. of Physics, 2009}
\examples{
data(adjmat1) # load input example data from package
adjout <- clearirrarcs( adjmat1)
}
