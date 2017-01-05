#' Perform all pre-processing steps
#'
#' Images must be of the standard format as provided in the NIST Ballistics and
#' Research Toolmark Database. The dimensions are 1944 x 2592 pixels, and images
#' are on a 255-grayscale. Here we use ring light images. This function performs
#' all the pre-processing steps, consisting of \itemize{ \item{Automatically
#' selecting the breechface marks} \item{Leveling the image} \item{Removing
#' circular symmetry} \item{Removing outliers and filtering.} }
#'
#' @param fileName location of image file. This function has been tested on
#'   \code{png} files.
#'
#' @return An image matrix, cropped to the smallest rectangle containing the
#'   valid breechface area. Non-breechface pixels are set to 0. After all
#'   pre-processing, possible pixel values are -510 to 510.
#'
#' @examples
#' \dontrun{
#' processedExample <- allPreprocess(system.file("extdata", "NBIDE R BF 118.png", package = "cartridges"))
#' }
#'
#' @export

allPreprocess <- function(fileName) {
    # step 1
    primer <- findPrimer(fileName) # primer is an EBImage image
    tmp <- findFP(fileName, primer = primer)

    # step 2
    out <- centerBFprimer(tmp, primer)
    tmp <- levelBF(out$centeredBF)

    # step 3
    tmp <- removeCircular(tmp)

    # step 4
    cropped <- cropBorders(tmp, out$centeredPrimer)
    tmp <- outlierRejection(cropped)
    tmp <- inpaint_nans(tmp)

    nonBF <- is.na(cropped)
    out <- gaussianFilter(tmp, nonBF)
    return(out)
}
