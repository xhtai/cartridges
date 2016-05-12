#' Read in a 2592 x 1944 TIFF cartridge image
#'
#' Images from the NIST Ballistics and Research Toolmark
#' Database are in this standard format, with dimensions
#' 2592 x 1944, and on a 255-grayscale. Images are read
#' in, returning a matrix with pixel values taking values
#' 0-255. By default images are cropped to square (1919 x
#' 1919 pixels) for further processing.
#'
#' @param TIFFfilename location of TIFF file to be read in
#' @param crop a logical value indicating whether the image should
#' be cropped to 1919 x 1919 pixels
#' @return A 1919 x 1919 (cropped) or 2592 x 1944 matrix
#' containing pixel values from 0 to 255.
#' @examples
#' LL1_3 <- readCropTIFF(system.file("extdata", "LL1_3.tif",
#'  package="cartridges"))
#' @export
#'
readCropTIFF<-function(TIFFfilename,crop=TRUE){
    if (!requireNamespace("rgdal", quietly = TRUE)) {
        stop("rgdal needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (crop==TRUE) {
        ret<-raster::as.matrix(raster::raster(TIFFfilename))[13:1931,337:2255]
    } else {
        ret<-raster::as.matrix(raster::raster(TIFFfilename))
    }
    return(ret)
}


#' Shift a 1919 x 1919 image
#'
#' @param image image that is to be shifted
#' @param centeri i-coordinate of desired center
#' @param centerj j-coordinate of desired center
#'
#' @return A shifted 1919 x 1919 image. Pixel values
#' that are not available are set to 0.
#'
#' @examples
#' shifted <- shiftedImage(LL1_3, 910, 910)
#' @export
#'
shiftedImage<-function(image,centeri,centerj){
    # default center is at (959,959)
    # convert desired ij center into dx and dy
    dy<-centeri-959
    dx<-959-centerj

    newGrid<-matrix(0,nrow=nrow(image),ncol=ncol(image))
    if (dx>=0 & dy>=0){
        newGrid[1:(1919-dy),(1+dx):1919]<-image[(1+dy):1919,1:(1919-dx)]
    } else if (dx>=0 & dy<0){
        newGrid[(1-dy):1919,(1+dx):1919]<-image[1:(1919+dy),1:(1919-dx)]
    } else if (dx<0 & dy>=0) {
        newGrid[1:(1919-dy),1:(1919+dx)]<-image[(1+dy):1919,(1-dx):1919]
    } else if (dx<0 & dy<0){
        newGrid[(1-dy):1919,1:(1919+dx)]<-image[1:(1919+dy),(1-dx):1919]
    }
    return(newGrid)
}


#' Get fitted image from basis function coefficients
#'
#' From basis function coefficients, produce a matrix of
#' pixel values.
#'
#' @param basisCoefficients data frame of basis function
#' coefficients, such as those produced by
#' \code{fitBasis()}
#' @param basis list of basis functions of the same length
#' as the number of rows in \code{basisCoefficients}. Can
#' be produced by \code{getBasisFunctions()}.
#' @param dimension dimension of output image. This has to
#' match the number of basis functions in \code{basis}.
#'
#' @return A matrix of pixel values
#' @examples
#' LL1_3_basis <- fitBasis(LL1_3,basis1919)
#' fittedImage <- getFittedImage(LL1_3_basis,basis1919,1919)
#' @export

getFittedImage<-function(basisCoefficients,basis,dimension){
    fittedImage<-matrix(NA,nrow=dimension,ncol=dimension)
    for (k in 1:length(basis)){
        ij<-basis[[k]]
        fittedImage[cbind(ij$i,ij$j)]<-basisCoefficients$coef[k]
    }
    return(fittedImage)
}


#' Plot a square image
#'
#' @param image image to be plotted. Image has to be square.
#' @param type either \code{"original"} for images on the
#' original scale (pixel values 0-255), or \code{"residuals"}
#' for residual pixel values (-255 to 255).
#' @param grayscale logical value indicating whether or not
#' the grayscale is to be plotted.
#' @param main title for plots with grayscale
#'
#' @examples
#' plotImage(LL1_3,"original",grayscale=FALSE)
#' @export


plotImage<-function(image,type,grayscale,main){
    dimension<-nrow(image)
    if (grayscale==FALSE){
        par(mar=c(0,0,0,0))
        plot(c(1,dimension), c(1,dimension),type="n",xaxt="n",yaxt="n",xlab = "", ylab = "",yaxs="i",xaxs="i")
        if (type=="original"){
            rasterImage(image/255,1,1,dimension,dimension,interpolate=FALSE)
        } else if (type=="residuals"){
            rasterImage((image+255)/510,1,1,dimension,dimension,interpolate=FALSE)
        }
    } else if (grayscale==TRUE){
        par(mar=c( 2.1, 2.1, 6.1, 6.1))
        plot(c(1,dimension), c(1,dimension),type="n",xaxt="n",yaxt="n",xlab = "", ylab = "",yaxs="i",xaxs="i",main=main)
        if (type=="original"){
            rasterImage(as.raster(image/255),1,1,dimension,dimension,interpolate=FALSE)
            fields::image.plot(legend.only=TRUE, zlim=c(0,255),col=gray(0:255/255))
        } else if (type=="residuals"){
            rasterImage(as.raster((image+255)/510),1,1,dimension,dimension,interpolate=FALSE)
            fields::image.plot(legend.only=TRUE, zlim=c(-255,255),col=gray(0:255/255),                  axis.args=list(at=seq(-250, 250, 50),
                                       labels=seq(-250, 250, 50),
                                       cex.axis=0.8))
        }
     }
}


