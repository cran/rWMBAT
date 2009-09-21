\name{TestCases}
\alias{TestCases}
\title{Uses Bayes Rule To Classify A Case}
\description{
Tests each of a set of data vectors by looking up P(data|class) in
a probability table, then finding P(case|class) by multiplying each
of those values in a product.  Then uses Bayes' rule to calculate
P(class|data) for each possible value of class.  Reports this as an
array of class probabilities for each case.
}
\usage{
TestCases(p, prior, data)
}
\arguments{
  \item{p}{3-D double array of probabilities (c,d,v).  The first dimension is the class, the second is the data value, the third is the variable number. The entry is P(var v=value d | class=value c)}
  \item{prior}{double columns vector, prior probabilities for each cases in data}
  \item{data}{double array of discrete integer (1:n) values, cases in rows and variables in columns }
}
\value{2-D double array whose value is P(class=c|data) for each case. Cases are in rows, class in cols}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(ptable, prior, testgrpbin) # load example input data from package
classprobtable <- TestCases (ptable, prior, testgrpbin)
         }
