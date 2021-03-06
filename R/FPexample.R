#' Matrix of pixel values after selecting only the breechface impression
#'
#' "NBIDE R BF 118.png" after selecting only the breechface impression region,
#' as returned by \code{findFP("NBIDE R BF 118.png", primer = primerExample)}.
#' This is part of the first pre-processing step. Areas that are not part of the
#' breechface impression are set to 0, and the remaining pixel values are left
#' as is.
#'
#' @format A 1919 x 1919 matrix, in the original orientation, i.e. a cropped
#'   version of the original 1944 x 2592 image.

#' @source \url{https://tsapps.nist.gov/NRBTD}
"FPexample"
