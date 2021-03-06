% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step4_outlier_filter.R
\name{outlierRejection}
\alias{outlierRejection}
\title{Remove outliers and fill them in with NA}
\usage{
outlierRejection(input)
}
\arguments{
\item{input}{image matrix}
}
\value{
A matrix with non-breechface areas set to 0, and outliers set to NA.
}
\description{
For each pixel, we consider a local patch consisting all pixels in a 21 x 21
pixel matrix, centered at the pixel being considered. Outliers are defined as
pixels that have values larger than 3 standard deviations from the mean of
the pixel values in its local patch. This is translated from code by Joseph
Roth and Zach Richardson.
}
\examples{
\dontrun{
# first run this:
croppedExample <- cropBorders(removedExample, centeredExample$centeredPrimer)
# then
outlierNAexample <- outlierRejection(croppedExample)
}

}
