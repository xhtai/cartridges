#' Center image
#'
#' Given a 1919 x 1919 image matrix, and its binary primer image (in
#' \code{EBImage} format), computes the center of the image as the center of the
#' primer region, and returns a 1769 x 1769 image matrix where the center of the
#' primer is located at (884, 884).
#'
#' @param inputImage 1919 x 1919 matrix
#' @param primer binary image in \code{EBImage} format, as returned by
#'   \code{findPrimer}
#'
#' @return A list with two items: 1) \code{centeredBF}: A 1769 x 1769 matrix,
#'   with pixel values ranging from 0-255. Areas that are not part of the
#'   breechface impression are set to NA. 2) \code{centeredPrimer}: A 1769 x
#'   1769 binary matrix denoting the primer region, after centering. This is to
#'   be used in a later step.
#' @examples
#' \dontrun{
#' centeredExample <- centerBFprimer(FPexample, primerExample)
#' }
#'
#' @export

centerBFprimer <- function(inputImage, primer) {
    primer <- t(primer@.Data)
    indices <- which(primer == 1, arr.ind = TRUE)
    # centers are mean(indices[, "row"]), mean(indices[, "col"])
    centeri <- round(mean(indices[, "row"]))
    centerj <- round(mean(indices[, "col"]))

    centeredBF <- matrix(0, nrow = 1769, ncol = 1769)
    centeredBF[51:1719, 51:1719] <- inputImage[(centeri - 834):(centeri + 834), (centerj - 834):(centerj + 834)] # here we are assuming that biggest BF is < 1669px diameter
    centeredBF[centeredBF == 0] <- NA

    centeredPrimer <- matrix(0, nrow = 1769, ncol = 1769)
    centeredPrimer[51:1719, 51:1719] <- primer[(centeri - 834):(centeri + 834), (centerj - 834):(centerj + 834)]
    centeredPrimer[centeredPrimer == 0] <- NA

    return(list(centeredBF = centeredBF, centeredPrimer = centeredPrimer))
}

#' Level image
#'
#' Given a 1769 x 1769 image matrix, fits a plane through the image and returns
#' the residuals. Areas that are not part of the breechface impression are set
#' to NA. Removing residuals means that possible values are -255 to 255.
#'
#' @param inputImage 1769 x 1769 matrix
#'
#' @return A leveled 1769 x 1769 image matrix.
#'
#' @examples
#' \dontrun{
#' # first run this:
#' centeredExample <- centerBFprimer(FPexample, primerExample)
#' # then:
#' leveledExample <- levelBF(centeredExample$centeredBF)
#' }
#'
#' @export

levelBF <- function(inputImage) {
    dataframe <- expand.grid(1:1769, 1:1769, KEEP.OUT.ATTRS = FALSE)
    center <- floor(1769/2) # floor in case dim is odd
    names(dataframe) <- c("row", "col")
    dataframe$pixelValue <- inputImage[cbind(dataframe$row, dataframe$col)]
    dataframe$x <- dataframe$col - center
    dataframe$y <- center - dataframe$row

    fitPlane <- lm(pixelValue ~ x + y, data = dataframe)
    dataframe$leveled <- NA
    dataframe$leveled[!is.na(dataframe$pixelValue)] <- dataframe$pixelValue[!is.na(dataframe$pixelValue)] - fitPlane$fitted.values
    leveled <- matrix(dataframe$leveled, nrow = 1769)

    return(leveled)
}


