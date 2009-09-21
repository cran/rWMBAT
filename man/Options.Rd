\name{Options}
\alias{Options}
\docType{data}
\title{Example Data}
\description{Logical 6x1 array}
\usage{data(Options)}
\details{Options are
          1. Normalize on population total ion count (sum across rows)
          2. Remove negative data values by setting them to zero
          3. After normalizing, before binning, average cases with same ID
          4. NOT USED - SET TO FALSE
          5. Take log(data) prior to binning.  Negative values set to 1.
          6. NOT USED - SET TO FALSE
}
\examples{
data(Options)
}
\keyword{datasets}
