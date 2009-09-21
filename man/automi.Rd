\name{automi}
\alias{automi}
\title{Automi Finds A Threshold For Randomized MI(V C)}
\description{
 Finds the threshold of a data set's mutual information with a class
 vector, above which a variable's MI(class, variable) can be
 expected  to be significant. 
}
\usage{
automi(data, class, repeats)
}
\arguments{
  \item{data}{double array of discrete integer (1:n) values, cases in rows and variables in columns.}
  \item{class}{double column vector, also 1:n. Classification of each case.}
  \item{repeats}{integer, the number of times to repeat the randomization}
}
\details{
  The threshold for mi (significance 
 level) is found by taking the data set and randoomizing the class
 vector, then calculating MI(CV) for all the variables. This is
 repeated a number of times. The resulting list of length (repeats
 *variables) is sorted,  and the 99th percentile max MI is taken
 as the threshold.
}
\value{
 a threshold for randomized MI(V C) 
       }
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009.}
\examples{
       data(traingrpbin, traingrpclass) # load the example input data from package
       threshold <- automi( traingrpbin, traingrpclass, 10 )
          }


