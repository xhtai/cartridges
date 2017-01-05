#' Remove circular symmetry of an image
#'
#' Given a 1769 x 1769 image matrix, fits a circularly symmetric model on the
#' image and returns the residuals. Areas that are not part of the breechface
#' impression are set to NA. If this is run after \code{levelBF}, removing
#' residuals again means that possible values are -510 to 510.
#'
#' @param inputImage 1769 x 1769 matrix
#'
#' @return A 1769 x 1769 image matrix.
#'
#' @examples
#' \dontrun{
#' # first run this:
#' centeredExample <- center(FPexample, primerExample)
#' leveledExample <- levelBF(centeredExample$centeredBF)
#' # then
#' removedExample <- removeCircular(leveledExample)
#' }
#'
#' @export

removeCircular <- function(inputImage) {
    basis <- fitBasis(inputImage, basis1769)
    basis$basisFnNumber <- 1:236463
    out <- loess(basis$coef ~ basis$basisFnNumber)
    fittedOut <- fitted(out)
    basis$fitted <- NA
    basis$fitted[!is.na(basis$coef)] <- fittedOut

    ####### get fitted and residuals
    fitted <- getFitted(basis$fitted, basis1769, 1769)
    residuals <- inputImage - fitted

    return(residuals)
}

#' Get ij coordinates for circularly symmetric basis functions
#'
#' A circularly symmetric image has the same pixel values for pixels that are
#' the same distance from the center. A circularly symmetric basis spans the
#' space of circularly symmetric matrices. Each matrix in the basis takes the
#' value 1 for pixels that are the same distance from the center, and zero
#' otherwise. Basis are enumerated from center outwards. Functions that take an
#' ij coordinate as an input and return the value 0 or 1 are termed circularly
#' symmetric basis functions.
#'
#' @param dimension size of image (\code{dimension} x \code{dimension}), where
#'   \code{dimension} is odd
#'
#' @return A list with length = number of basis functions for an image of size
#'   \code{dimension} x \code{dimension}. Each list item is a matrix whose rows
#'   are the ij coordinates of pixels in each basis function, and i2j2 is the
#'   square of the distance from the center of the pixels in each basis.
#' @examples
#' basis3 <- getBasisFunctions(3)
#' @export

getBasisFunctions <- function(dimension) {
    m <- dimension
    n <- (m - 1)/2
    values <- expand.grid(i = -n:n, j = -n:n, KEEP.OUT.ATTRS = FALSE)
    values$i2j2 <- values$i^2 + values$j^2
    values$i <- values$i + (m + 1)/2
    values$j <- values$j + (m + 1)/2
    if (requireNamespace("plyr", quietly = TRUE)) {
        basis <- plyr::dlply(values, 'i2j2', identity)
        attr(basis, 'split_type') <- NULL
        attr(basis, 'split_labels') <- NULL
        for (i in 1:length(basis)) {
            attr(basis[[i]], 'vars') <- NULL
        }
    } else {
        basis <- split(values, f = values$i2j2)
        for (i in 1:length(basis)) {
            rownames(basis[[i]]) <- NULL
        }
    }

    return(basis)
}


#' Get statistics for each basis function
#'
#' Get statistics on pixel values in each basis function: mean of the pixel
#' values (basis function coefficient), residual sum of squares using fitted
#' coefficient, sum of pixel values.
#'
#' @param ij ij coordinates of the basis function. Input is typically one item
#'   of the list of basis functions
#' @param image image that statistics are to be calculated for
#' @param requestedStats vector with each element being a requested statistic,
#'   e.g. \code{c("coef","RSS")}. Possible options are \code{"numPixels"} for
#'   the number of pixels in the input basis function, \code{"coef"} for the
#'   coefficient as described above, \code{"RSS"} for the residual sum of
#'   squares, \code{"sum"} for the sum of pixel values, and \code{"max"} and
#'   \code{"min"} for the maximum and minimum pixel values.
#'
#' @return A list with the requested statistics in the order provided.
#' @examples
#' basis3 <- getBasisFunctions(3)
#' stats <- statisticsByBasisFunction(basis3[[2]],
#'          matrix(1:9, nrow = 3), c("numPixels", "coef"))
#' @export

statisticsByBasisFunction <- function(ij, image, requestedStats){
    ret <- vector(mode = "list", length = length(requestedStats))
    pixelValues <- image[cbind(ij$i, ij$j)]
    pixelValues <- pixelValues[!is.na(pixelValues)]

    for (i in 1:length(requestedStats)) {
        ret[[i]] <- switch(requestedStats[i],
                         numPixels = nrow(ij),
                         coef = mean(pixelValues),
                         RSS = sum((pixelValues - mean(pixelValues))^2),
                         sum = sum(pixelValues),
                         max = max(pixelValues),
                         min = min(pixelValues)
        )
    }
    return(ret)
}

#' Fit basis functions
#'
#' Fit basis functions to a centered image.
#'
#' @param image centered image. Image has to be square.
#' @param basis list of basis functions of the appropriate dimension, such as
#'   those produced by \code{getBasisFunctions()}.
#'
#' @return A data frame with the number of pixels in each basis function, basis
#'   function coefficient and residual sum of squares for each basis function.
#' @examples
#' basis3 <- getBasisFunctions(3)
#' sampleBasis <- fitBasis(matrix(1:9,nrow=3),basis3)
#' @export

fitBasis <- function(image, basis){
    imageBasis <- unlist(lapply(basis, function(x) statisticsByBasisFunction(x, image, c("numPixels", "coef", "RSS"))), use.names = FALSE)
    imageBasis <- data.frame(matrix(imageBasis, byrow = TRUE, ncol = 3))
    names(imageBasis) <- c("numPixels", "coef", "RSS")
    return(imageBasis)
}


#' Get fitted image from basis function coefficients
#'
#' From basis function coefficients, produce a matrix of pixel values.
#'
#' @param basisCoefficients vector of basis function coefficients, such as
#'   the \code{coef} column produced by \code{fitBasis()}
#' @param basis list of basis functions of the same length as the number of rows
#'   in \code{basisCoefficients}. Can be produced by \code{getBasisFunctions()}.
#' @param dimension dimension of output image. This has to match the number of
#'   basis functions in \code{basis}.
#'
#' @return A matrix of pixel values
#' @examples
#' \dontrun{
#' # first run this:
#' centeredExample <- center(FPexample, primerExample)
#' # then
#' basis <- fitBasis(centeredExample, basis1769)
#' fitted <- getFitted(basis$fitted, basis1769, 1769)
#' }
#' @export

getFitted <- function(basisCoefficients, basis, dimension) {
    fittedImage <- matrix(NA, nrow = dimension, ncol = dimension)
    for (k in 1:length(basis)) {
        ij <- basis[[k]]
        fittedImage[cbind(ij$i, ij$j)] <- basisCoefficients[k]
    }
    return(fittedImage)
}


