\name{WMBAT}
\alias{WMBAT}
\title{The William and Mary Bayesian Analysis Tool}
\description{
       WMBAT takes an array of mass spec peak intensities, a vector
       describing which of two classes each sample belongs to, 
       and other information and builds and assesses a Bayesian network
       after selecting features (peaks) from within the data array that
       are diagnostic of the class. The primary output is an adjacency
       matrix describing the resulting Bayesian network.
}
\usage{
WMBAT(tofListMetaData, alignedPeakList, Options, nfold, repeats, threshold)
}
\arguments{
  \item{alignedPeakList}{contain information about mass spec peak intensity values}
  \item{tofListMetaData}{a list containing information about every spectra, the information related to the spectra we needed in this package is class and ID.
                         Class,Integer vector, values 1 or 2 identifying the class of each case, such as "disease, non-disease"
                         ID,Double one or two column array contains the sample IDfor each case. Second column is optional and would identifyreplicates of the same sample. }
  \item{Options}{Logical 6x1 array. Options are
          1. Normalize on population total ion count (sum across rows)
          2. Remove negative data values by setting them to zero
          3. After normalizing, before binning, average cases with same ID
          4. NOT USED - SET TO FALSE
          5. Take log(data) prior to binning.  Negative values set to 1.
          6. NOT USED - SET TO FALSE }
  \item{nfold}{integer, the "n" in n-fold cross validation (integer 4-10). 10 is recommended}
  \item{repeats}{integer, times to repeat the whole process (e.g.re-crossvalidate). 100 is recommended }
  \item{threshold}{double,factor by which the maximum "random" MI is multiplied to
           find the minimum "significant" MI (double, 1.0-5.0). We
           recommend starting with 1 and increasing until a "reasonable"
           number of diagnostic peaks is reached and error rates are
           minimized. This setting is dependant on the data and the
          correlations between variables }
}
\value{
  \item{IntOut }{double Intensities input array, after processing by the various options selected by the logical Options above}
  \item{IDOut }{double vector,the ID number of each row in the IntOut array. With no replicate averaging, each ID will be preserved (but reformatted) from the input.  With replicate averaging, only the primary ID number remains.}
  \item{PredClass }{double matrix,the predicted class of each case, during each of the trials (from input "repeats") }
  \item{Class2Vars }{vector whose ith value is the fraction of times peak i (from the vector MZ) was selected as being connected to the class. The maximum times it could have been selected was nfold*repeats}
  \item{Var2Vars }{integer array whose (i,j) entry is the fraction of times a second level link was found from peak i to peak j, when peak i was connected to the class, as found in SumLvl1}
  \item{MetaVars }{integer array whose (i,j) entry is the fraction of times a metavariable was created using peak i and peak j and stored in the level}
  \item{TrialErr }{double vector, the error rate for each of the "repeats" possible trials. Records the percentage of cases where PredClass was not equal to the input Class.}
}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\note{CALLED FUNCTIONS
DoTheMath(Learns a Bayesian Network from the data)
}
\examples{
data(tofListMetaData, alignedPeakList, Options, nfold, repeats, threshold)#load input example data from the package 
result <- WMBAT (tofListMetaData, alignedPeakList, Options, nfold, repeats, threshold)# Running this example may take some time, about one hour on a Intel(R) Core(TM)2 CPU, T55000 @ 1.66GHz,1.49 BG of RAM laptop.
         }
