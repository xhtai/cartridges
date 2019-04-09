#' Matrix of pixel values after the first three pre-processing steps
#'
#' "NBIDE R BF 118.png" after the first three pre-processing steps: selecting
#' the breechface impression, leveling, and removing circular symmetry. The
#' following is the full set of code to obtain this matrix. \preformatted{
#' # step 1
#' primerExample <- findPrimer(system.file("extdata", "NBIDE R BF 118.png", package = "cartridges"))
#' FPexample <- findFP(system.file("extdata", "NBIDE R BF 118.png", package = "cartridges"), primer = primerExample)
#' # step 2
#' centeredExample <- centerBFprimer(FPexample, primerExample)
#' leveledExample <- levelBF(centeredExample$centeredBF)
#' # step 3
#' removedExample <- removeCircular(leveledExample)
#' }
#' Areas that are not part of the breechface impression are set to NA. If this is
#' run after \code{levelBF}, removing residuals again means that possible values
#' are -510 to 510.
#'
#' @format A 1769 x 1769 image matrix.

#' @source \url{https://tsapps.nist.gov/NRBTD}
"removedExample"
