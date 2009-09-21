\name{FindProbTables}
\alias{FindProbTables}
\title{Estimates The Probabilities P(class=c|data=D)}
\description{
Input a training group of data arranged with cases in rows and 
variables in columns, as well as the class value c for that vector. 
Each case represents a data vector V.  For each possible data value 
vi, and each variable Vi, it calculates P(C=c|Vi=vi) and stores 
that result in a 3-D table.  The table is arranged with the 
dimensions (class value, data value, variable number).
}
\usage{
FindProbTables(data, class)
}
\arguments{
  \item{data}{double array of discrete integer  values, cases in rows and variables in columns}
  \item{class}{double column vector. Classification of each case}
}
\value{3-D array whose (c,d,v) value is P(class=c|data=p) for variable v}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(finaldata, traingrpclass) #load the example input data from package
probtable <- FindProbTables(finaldata, traingrpclass)
}
