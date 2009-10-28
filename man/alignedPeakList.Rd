\name{alignedPeakList}
\alias{alignedPeakList}
\docType{data}
\title{Example Data}
\description{
One of the input example data for main function "WMBAT", a list containing information of time labels for peaks and peak intensity values for every spectrum.

}
\usage{data(alignedPeakList)}
\examples{
data(alignedPeakList)
names(alignedPeakList)

plot(alignedPeakList$peaks,alignedPeakList$data[[1]]$Intensities,"l")
#"alignedPeakList$peaks" is an array containing time labels for peaks.
#"alignedPeakList$data[[1]]$Intensities" is an array containing peak intensity values of the 1st spectrum.

}
\keyword{datasets}
