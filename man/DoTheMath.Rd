\name{DoTheMath}
\alias{DoTheMath}
\title{Takes A Data Set And Performs Feature Selection}
\description{
 DoTheMath takes a data array, class vector, and other information
 and builds and assesses a Bayesian network after selecting features
 from within the data array.  
}
\usage{
DoTheMath(InputStructure)
}
\arguments{
           \item{InputStructure}{a list containing following inputs:
           Class:Vector of length "cases", with discrete values identifying class of each case (may be integer).
           ID:double patient ID array of length cases, with one or more cols.
           MZ:vector of length "variables" holding labels for variables.
           Options:logical 6x1 array. Options are:
              1. Normalize on population total ion count (sum across rows).
              2. Remove negative data values by setting them to zero.
              3. After normalizing, before binning, average cases with same ID.
              4. Find the MI threshold by randomization.
              5. Take log(data) prior to binning.  Negative values set to 1.
              6. Remove Low Signal cases.
                 NOT DONE 3 Bin (2 Bin if False).
           n:integer,the "n" in n-fold cross validation.
           repeats:integer,times to repeat the whole process (e.g. re-crossvalidate).
           threshold:double, factor by which the maximum "random" MI us multiplied to find the minimum "significant" MI (double, 1.0-5.0).}
           }
\details{
  This is the umbrella script that loops a specified number of times
 (see "repeats" above), each time doing a full n-fold cross
 validation and recording the results.  All input and output data
 are stored in a single data structure, described below.
}
\value{
       all the fields of InputStructure, plus, ErrorRate:Vector containing misclassification rate for each repeat.KeyFeatures:Index to vector MZ that identifies features selected
       }
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\note{CALLED FUNCTIONS
       InitialProcessing(Applies the options listed above)
       BuildBayesNet(Learns a Bayesian Network from the training data)
       ChooseMetaVars(Combines variables that may not be physically separate molecules)
       TestCases(Given the BayesNet, tests the "test group" to determine the probability of being in each class)
       opt3bin(Discretizes continuous data into 3 bins, optimizing MI)
       FindProbTables(Learns the values P(C,V) for each variable)
}
\examples{
data(In) #load input example data from the package
OutputDataStructure <- DoTheMath (In) ## Running this example may take some time, about 10 minutes on a Intel(R) Core(TM)2 CPU, T55000 @ 1.66GHz,1.49 BG of RAM laptop.
         }
