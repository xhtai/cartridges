#' Matrix of pixel values of processed image from the NBIDE study.
#'
#' A matrix of pixel values after processing, from an example image taken from
#' the NBIDE study, one of the data sets available in the NIST Ballistics
#' Toolmark Research Database. This image is a 2D breechface ring light image,
#' and the original filename in the NIST download is "NBIDE R BF 129.png". The
#' cartridge being imaged is the second test fire from a Ruger gun (gun 3 in the
#' study), using PMC ammunition (NBIDE R BF 118.png is the first test fire).
#' This image was processed using the function \code{\link{allPreprocess}}.
#' Theoretically, possible values are -510 to 510, but in actual fact the values
#' are much smaller in magnitude. Pixels that are not part of the breechface
#' marks are set to 0.
#'
#' @format A matrix with dimensions 1347 x 1341 Each entry contains pixel
#'   values after processing.
#' @source \url{https://tsapps.nist.gov/NRBTD}
"processedExample2"
