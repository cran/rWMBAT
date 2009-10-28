\name{BuildBayesNet}
\alias{BuildBayesNet}
\title{Selects Features And Metafeatures Based On Mutual Info}
\description{BuildBayesNet selects features and metafeatures based on mutual info}
\usage{BuildBayesNet(data, class, ffactor, drop)}
\arguments{
  \item{data}{ integer array containing the data used to build the Bayes net,cases in rows, variables in cols}
  \item{class}{ double column vector, the known class variable for each case}
  \item{ffactor}{ multiple of auto MI to use to threshold C->V connections }
  \item{drop}{ double, MI loss percentage threshold for testing independence. Set to.75 and adjust to filter too few/too many variable-to-variable connections}
}
\details{
  This function takes a set of training data and an additional
  variable called "class" and tries to learn a Bayesian Network
  Structure by examining Mutual Information.
}
\value{
  a matrix of zeros and ones, where one in row i, column j
  denotes a directed link in a Bayesian network between 
  variable i and variable j. The class variable is the last
  row/column.
}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke , College of William and Mary, Dept. of Physics, 2009.}
\note{CALLED FUNCTIONS:
       automi(finds an MI threshold based on data)
       findmutualinfos(finds all values MI(VC), MI(VV) and MI(VC|V))
}
\examples{
data(traingrpbin, traingrpclass, ff, drop) #load example data in package 
adjmat <- BuildBayesNet( traingrpbin, traingrpclass, ff, drop )
}
