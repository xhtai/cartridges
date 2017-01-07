#' Read cartridge case image and find primer region
#'
#' Images must be of the standard format as provided in the NIST Ballistics and
#' Research Toolmark Database. The dimensions are 1944 x 2592 pixels, and images
#' are on a 255-grayscale. Here we use ring light images. This function finds
#' the primer region, which is the circular region containing the breechface
#' marks and firing pin impression. By default images are cropped to 1919 x 1919
#' pixels. For more information on the steps involved in finding this region,
#' refer to the GitHub page.
#'
#' This function uses the \code{EBImage} package, which stores the image as a
#' 2592 x 1944 matrix. The primer region that is returned is similarly an
#' \code{EBImage} image, where the primer region is stored in a matrix with
#' binary values. If using this matrix directly on the image matrix in the
#' original orientation (1944 x 2592), note that the should first be transposed.
#'
#'
#' @param fileName location of image file. This function has been tested on
#'   \code{png} files.
#'
#' @return A 1919 x 1919 binary image in \code{EBImage} format. Values 1
#'   indicate the primer region.
#' @examples
#' \dontrun{
#' primerExample <- findPrimer(system.file("extdata", "NBIDE R BF 118.png", package = "cartridges"))
#' }
#'
#' @export

findPrimer <- function(fileName) {
    ebimage <- EBImage::readImage(fileName)
    I <- ebimage[337:2255, 13:1931]

    sigma <- 25
    blurred <- EBImage::gblur(I, sigma, radius = 2 * ceiling(2 * sigma) + 1)

    equalized <- ifelse(blurred <= median(blurred), 0, 1)

    # imfill
    tmp <- EBImage::fillHull(equalized)

    # bwlabel -- choose the center region
    L = EBImage::bwlabel(tmp)
    ij <- matrix(c(959, 600, 1200, 959, 959, 960, 960, 960, 600, 1200), ncol = 2)

    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    tmp[L != Mode(L[cbind(ij)])] <- 0

    ##### dilate, fill and erode
    se <- EBImage::makeBrush(size = 139, shape = "disc") # should be odd
    tmp <- EBImage::dilate(tmp, se)
    filled <- EBImage::fillHull(tmp)
    se <- EBImage::makeBrush(size = 269, shape = "disc") # should be odd
    outer <- EBImage::erode(filled,se)
    return(outer)
}


#' Find and remove firing pin region
#'
#' After finding the primer region using \code{findPrimer()}, we now remove the
#' firing pin impression, returning a 1919 x 1919 matrix, where areas that are
#' not part of the breechface impression are set to 0. Pixel values range from 0
#' to 255. For more information on the steps involved in finding the firing pin
#' region, refer to the GitHub page.
#'
#' @param fileName location of image file. This function has been tested on
#'   \code{png} files.
#' @param primer binary image in \code{EBImage} format, as returned by
#'   \code{findPrimer}
#'
#' @return A 1919 x 1919 matrix, in the original orientation, i.e. a cropped
#'   version of the original 1944 x 2592 image.
#' @examples
#' \dontrun{
#' FPexample <- findFP(system.file("extdata", "NBIDE R BF 118.png", package = "cartridges"), primer = primerExample)
#' }
#'
#' @export
#' @import methods

findFP <- function(fileName, primer) {
    if (!requireNamespace("imager", quietly = TRUE) || !requireNamespace("purrr", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) {
        stop("This function uses a Canny edge detector, and this implementation requires the packages imager, purrr, and dplyr. Please install them. For an example, see help(FPexample).", call. = FALSE)
    }
    ebimage <- EBImage::readImage(fileName)
    J <- ifelse(ebimage <= median(ebimage), 0, 1)
    J <- J[337:2255, 13:1931]

    # canny
    sigma <- 22 # changed from 32
    blurred <- EBImage::gblur(J, sigma)

    tmp <- blurred@.Data
    imagerImage <- imager::as.cimg(tmp)
    out <- cannyEdges(imagerImage*255, q1 = .89*.4, q2 = .89)

    out[primer == 0] <- 0

    if (sum(out) < 2000) {
        thres <- .89
        while (sum(out) < 2000) {
            thres <- thres - .02
            out <- cannyEdges(imagerImage*255, q1 = thres*.4, q2 = thres)
            out[primer == 0] <- 0
        }
    }

    outtmp <- as.matrix(out)
    outebimage <- EBImage::Image(outtmp, colormode = 'Grayscale') # back to ebimage

    se <- EBImage::makeBrush(size = 399, shape = "disc") # should be odd
    tmp <- EBImage::dilate(outebimage, se)
    tmp <- EBImage::fillHull(tmp)
    se <- EBImage::makeBrush(size = 379, shape = "disc") # should be odd
    FP <- EBImage::erode(tmp,se)

    while (sum(FP == 1) < 675000) { #3000000 -- fix this -- percentage of BF to remove
        se <- EBImage::makeBrush(9, 'disc')
        FP <- EBImage::dilate(FP, se)
    }

    # second pass for FP
    newJ <- J
    newJ[primer == 0] <- 1 # 1 instead of 255
    newJ[FP == 1] <- 1 # 1 instead of 0 -- make the FP white not black

    sigma <- 2.5 # changed from 2 -- blur more
    blurred <- EBImage::gblur(newJ, sigma)
    tmp <- blurred@.Data
    imagerImage <- imager::as.cimg(tmp)
    tmp <- cannyEdges(imagerImage*255, q1 = .1, q2 = .98) # change q1 to .1 -- let more weakly connected components be found

    outtmp <- as.matrix(tmp) + FP@.Data # combine with earlier FP
    outtmp[outtmp > 1] <- 1
    outebimage <- EBImage::Image(outtmp, colormode = 'Grayscale') # back to ebimage

    se <- EBImage::makeBrush(139, 'disc')
    tmp <- EBImage::dilate(outebimage,se)
    tmp <- EBImage::fillHull(tmp)
    se <- EBImage::makeBrush(139, 'disc')
    FP <- EBImage::erode(tmp, se)

    newI <- ebimage[337:2255, 13:1931]
    newI[FP == 1] <- 0 #FP==1 (white) -- remove
    newI[primer == 0] <- 0 #BW3==0 (black) -- remove (outer regions)

    # convert to matrix, 0-255
    newI <- t(newI@.Data * 255)
    return(newI)
}

# ACKNOWLEDGEMENTS: the rest of the code in this file was adapted from Simon
# Barthelmé, "A loop-free canny edge detector"
# (http://dahtah.github.io/imager/canny.html)
#' @importFrom magrittr "%>%"

fillInit <- function(strong) {
    lab <- imager::label(strong, TRUE)*strong
    as.data.frame(lab) %>% dplyr::filter(value > 0) %>% dplyr::group_by(value) %>% dplyr::summarize(x = x[1], y = y[1])
}

#Starts a fill at each successive location, and accumulates the results
rescueFill <- function(strong, weak) {
    v <- strong
    v[weak == 1] <- .9
    loc <- fillInit(strong)
    #Transform the data.frame into a list of locations
    tmp <- dplyr::select(loc, -value)
    loc <- purrr::transpose(tmp)
    #Fold
    out <- purrr::reduce(loc, function(v,l) imager::bucketfill(v, l$x, l$y, color = 1, sigma = .1, high = TRUE), .init = v)
    imager::as.cimg(out == 1)
}

nonmax <- function(gr) {
    mag <- with(gr, sqrt(x^2 + y^2))
    ang <- with(gr, atan2(y, x))
    grs <- list(x = gr$x/mag, y = gr$y/mag)
    X <- imager::Xc(gr$x)
    Y <- imager::Yc(gr$y)
    val.bwd <- imager::interp(mag, data.frame(x = as.vector(X - grs$x),
                                     y = as.vector(Y - grs$y)))
    val.fwd <- imager::interp(mag,data.frame(x = as.vector(X + grs$x),
                                     y = as.vector(Y + grs$y)))

    throw <- (mag < val.bwd) | (mag < val.fwd)
    mag[throw] <- 0
    mag
}


cannyEdges <- function(im, q1, q2){
    mag <- imager::imgradient(im, "xy") %>% nonmax
    t1 <- q1*max(mag)
    t2 <- q2*max(mag)
    strong <- imager::as.cimg(mag > t2)
    weak <- imager::as.cimg(mag >= t1 & mag <= t2)

    rescueFill(strong, weak)

}

