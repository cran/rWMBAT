\name{InitialProcessing}
\alias{InitialProcessing}
\title{INITIALPROCESSING Inital Prep of Data from Signal Processing}
\description{Takes peak list which has been imported into R and prepares the peak list For Bayesian Analysis}
\usage{
InitialProcessing(StructIn)
}
\arguments{
  \item{StructIn}{List with the following double-typed arrays:
  Intensities:double n x m real-valued array with variables (peaks) in columns, cases (samples) in rows.
  MZ:double list of the labels (m/z value) for each of the variables.Must be the same size as the number of variables in Intensities.
  Class:intefer vector, classification of each sample (disease state)-- 1 or 2--must be the same size as the number of cases in Intensities.
  ID:double column array, case or patient ID number, same size as class.  May have second column, so each row is [ID1 ID2] where ID2 is replicate number.
  Options:logical Array of processing options with elements
           1. Normalize.
           2. Clip Data (remove negatives).
           3. Replicate Average.
           4. Auto threshold MI.
           5. Use Log of Data.
           6. Remove Low Signal cases.
           NOT DONE 3 Bin (2 Bin if False).}
                 
}
\value{
  \item{RawData}{Intensities as input}
  \item{ClipData}{RawData where all values less than 1 are set to 1}
  \item{NormData}{ClipData normalized by total ion count, i.e. divided by the sum of all variables for each case}
  \item{LogData}{Natural logarithm of NormData}
  \item{Class}{Same as input}
  \item{MZ}{Same as input}
  \item{ID}{SIngle column. If replicates are not averaged, the entries are now ID1.ID2. If replicates averaged, then just ID1}
  \item{DeltaMZ}{difference in peak m/z values to look for adducts}
  \item{RatioMZ}{ratios of m/z values ot look for satellites}
}
\references{http://kwkusc.people.wm.edu/dissertation/dissertation.htm }
\author{Karl Kuschner, Qian Si and William Cooke, College of William and Mary, Dept. of Physics, 2009}
\examples{
data(In) #load input example data from the package
StructOut <- InitialProcessing(In)
}
