#' Get ij coordinates for circularly symmetric basis functions
#'
#' A circularly symmetric image has the same pixel values
#' for pixels that are the same distance from the center.
#' A circularly symmetric basis spans the space of circularly
#' symmetric matrices. Each matrix in the basis takes the
#' value 1 for pixels that are the same distance from the
#' center, and zero otherwise. Basis are enumerated from
#' center outwards.
#'
#' For a 3 x 3 matrix, there are a total of 3 basis. The
#' first has 1 for the (2, 2) pixel and 0 otherwise. The
#' second basis has 1's for coordinates (1, 2), (2, 1),
#' (2, 3) and (3, 2), and 0 otherwise. The third has
#' 1's for coordinates (1, 1), (1, 3), (3, 1), (3, 3),
#' and 0's otherwise. Functions that take an ij coordinate
#' as an input and return the value 0 or 1 are termed
#' circularly symmetric basis functions.
#'
#' @param dimension size of image (\code{dimension} x
#' \code{dimension}), where \code{dimension} is odd
#'
#' @return A list with length = number of basis functions
#' for an image of size \code{dimension} x \code{dimension}.
#' Each list item is a matrix whose rows are the ij
#' coordinates of pixels in each basis function, and i2j2
#' is the distance from the center of the pixels in each
#' basis.
#' @examples
#' basis3 <- getBasisFunctions(3)
#' @export

getBasisFunctions<-function(dimension){
    m<-dimension
    n<-(m-1)/2
    values<-expand.grid(i=-n:n,j=-n:n,KEEP.OUT.ATTRS=FALSE)
    values$i2j2<-values$i^2+values$j^2
    values$i<-values$i+(m+1)/2
    values$j<-values$j+(m+1)/2
    basis<-split(values,f=values$i2j2)
    for (i in 1:length(basis)){
        rownames(basis[[i]])<-NULL
    }
    return(basis)
}


#' Get subset of basis functions
#'
#' Get a subset of basis functions produced by
#' \code{getBasisFunctions()}, subsetting based
#' on distance from center and the number of pixels
#' in each basis function.
#'
#'
#' @param originalBasis list object produced by
#' \code{getBasisFunctions()}
#' @param distFromCenter vector of length 2 indicating
#' the minimum and maximum distances from the center
#' @param numPixels integer value indicating the required
#' number of pixels in each basis function
#' @return A subset of the original basis with basis
#' functions of the required distances from the center
#' and number of pixels.
#' @examples
#' basis3 <- getBasisFunctions(3)
#' basis3_rad1_4pixels <- subsetBasis(basis3,c(1,1),4)
#'
#' @export

subsetBasis<-function(originalBasis,distFromCenter,numPixels){
    if (!missing(distFromCenter)){
        radii<-unlist(lapply(originalBasis,FUN=function(x) {sqrt(x[1,"i2j2"])}),use.names=FALSE)
        indicesToKeep<-which(radii>=distFromCenter[1] & radii<=distFromCenter[2])
        basis<-originalBasis[indicesToKeep]
    }
    if (!missing(numPixels)){
        pixels<-unlist(lapply(basis,FUN=nrow),use.names=FALSE)
        basis<-basis[which(pixels==numPixels)]
    }
    return(basis)
}



#' Get statistics for each basis function
#'
#' Get statistics on pixel values in each basis function:
#' mean of the pixel values (basis function coefficient),
#' residual sum of squares using fitted coefficient,
#' sum of pixel values.
#'
#'
#' @param ij ij coordinates of the basis function.
#' Input is typically one item of the list of basis
#' functions
#' @param image image that statistics are to be calculated
#' for
#' @param requestedStats vector with each element
#' being a requested statistic, e.g. \code{c("coef","RSS")}.
#' Possible options are \code{"numPixels"} for the number
#' of pixels in the input basis function, \code{"coef"} for
#' the coefficient as described above, \code{"RSS"} for
#' the residual sum of squares, and \code{"sum"} for
#' the sum of pixel values.
#'
#' @return A list with the requested statistics in the
#' order provided.
#' @examples
#' basis3 <- getBasisFunctions(3)
#' stats <- statisticsByBasisFunction(basis3[[2]],
#'          matrix(1:9,nrow=3),c("numPixels","coef"))
#' @export
#'
statisticsByBasisFunction<-function(ij,image,requestedStats){
    ret<-vector(mode="list",length=length(requestedStats))
    pixelValues<-image[cbind(ij$i,ij$j)]

    for (i in 1:length(requestedStats)){
        ret[[i]]<-switch(requestedStats[i],
                         numPixels=nrow(ij),
                         coef=mean(pixelValues),
                         RSS=sum((pixelValues-mean(pixelValues))^2),
                         sum=sum(pixelValues)
                         )
    }
    return(ret)
}

#' Fit basis functions
#'
#' Fit basis functions to a centered image.
#'
#' @param image centered image
#' @param basis list of basis functions of the appropriate
#' dimension, such as those
#' produced by \code{getBasisFunctions()}.
#'
#' @return A data frame with the number of pixels in
#' each basis function, basis function coefficient and
#' residual sum of squares for each basis function.
#' @examples
#' basis3 <- getBasisFunctions(3)
#' sampleBasis <- fitBasis(matrix(1:9,nrow=3),basis3)
#' @export

fitBasis<-function(image,basis){
    imageBasis<-unlist(lapply(basis,function(x) statisticsByBasisFunction(x,image,c("numPixels","coef","RSS"))),use.names=FALSE)
    imageBasis<-data.frame(matrix(imageBasis,byrow=TRUE,ncol=3))
    names(imageBasis)<-c("numPixels","coef","RSS")
    return(imageBasis)
}

