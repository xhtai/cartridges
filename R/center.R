
#' Rough search for the center of a 1919 x 1919 image
#'
#' This method is intended only as a rough search for a starting pixel location
#' to be used by \code{gridSearch()}. An appropriate center would fit the ring
#' in the image well, and the goal is to find a center producing the longest
#' continuous string of fitted small coefficients where the ring is expected to
#' be. An approximate location of the ring is 450-1000 pixels from the center,
#' so only the subset of basis functions corresponding to this distance is
#' recommended.
#'
#' Along an input row (recommended rows are 940 and 980), we consider 201
#' possible center locations from columns 860 to 1060. Then after selecting an
#' appropriate column, we consider 201 possible rows (100 above and 100 below
#' the input row). We recommend use of a further subset of basis functions
#' (those with 8 pixels in each basis function, since these are the most common)
#' for speed. Also, instead of fitting each basis function (computing the mean
#' of the pixel values) and searching for the longest string of small values, we
#' convert the images from 256-grayscale to binary values, with values from 0 to
#' 127 taking the value 0 and those from 128 to 255 taking the value 1. Now at
#' each possible center, we compute the sum of pixel values of each basis
#' function, and find the center resulting in the longest continuous string of
#' zeros.
#'
#'
#' @param image image that rough center is to be determined for
#' @param basis subset of basis functions of the appropriate dimension (1919 x
#'   1919). Recommended restrictions are basis functions with distances from
#'   center of 450-1000 pixels, and with 8 pixels in each basis function.
#' @param startingRow row to be searched
#'
#' @return A data frame with 402 rows, each row with the coordinates of the
#'   center considered, and the maximum number of continuous zeros in the sum of
#'   binary pixels for each basis function.
#'
#' @examples
#' \dontrun{
#' basis1919_8pixels_rad450_1000 <- subsetBasis(basis1919,c(450,1000),8)
#' center <- roughCenter(LL1_3, basis1919_8pixels_rad450_1000, 940) # warning:
#' this is a long computation (may take half an hour or more, depending on your
#' system)
#' }
#' @export
#'

roughCenter<-function(image,basis,startingRow){
    image<-floor(image/(2^7))

    columnNumber<-860
    rowNumber<-startingRow

    runLengths<-data.frame(matrix(NA,nrow=402,ncol=3))
    names(runLengths)<-c("i","j","maxRun")

    for (j in 1:201){
        shifted<-shiftedImage(image,rowNumber,columnNumber)
        allSums<-unlist(lapply(basis,function(x) statisticsByBasisFunction(x,shifted,c("sum"))),use.names=FALSE)

        tmp<-rle(allSums)
        runLengths[j,]<-c(rowNumber,columnNumber,max(tmp$lengths))
        # do a check for values: if not 0, put max run as 0
        if (tmp$values[which.max(tmp$lengths)]!=0) {runLengths[j,3]<-0}

        columnNumber<-columnNumber+1
    }

    columnNumber<-runLengths[which.max(runLengths$maxRun),"j"]
    # note: if tied, go with first one
    rowNumber<-startingRow-100

    for (i in 1:201){
        shifted<-shiftedImage(image,rowNumber,columnNumber)
        allSums<-unlist(lapply(basis,function(x) statisticsByBasisFunction(x,shifted,c("sum"))),use.names=FALSE)

        tmp<-rle(allSums)
        runLengths[201+i,]<-c(rowNumber,columnNumber,max(tmp$lengths))
        if (tmp$values[which.max(tmp$lengths)]!=0) {runLengths[201+i,3]<-0}

        rowNumber<-rowNumber+1
    }

    return(runLengths)
}


#' Grid search for the center of a 1919 x 1919 image
#'
#' At a grid of potential centers, fit basis functions corresponding to the ring
#' (approximately 450-1000 pixels from the center) to the image, and calculate
#' the total residual sum of squares of the fit. Small values indicate a better
#' fit and a more appropriate center.
#'
#' @param image image that center is to be determined for
#' @param basis subset of basis functions of the appropriate dimension (1919 x
#'   1919), corresponding approximately to the outer ring. Recommended distances
#'   from center are 450-1000 pixels.
#' @param startingRow starting row to be searched
#' @param startingColumn starting column to be searched
#' @param plot logical value indicating whether a surface plot of total RSS
#'   values is desired
#'
#' @return A data frame with the ij coordinates of the center considered, and
#'   the total residual sum of squares of the fit.
#'
#' @examples
#' \dontrun{
#' basis1919_rad450_1000<-subsetBasis(basis1919,c(450,1000))
#' RSSdtf<-gridSearch(LL1_3,basis1919_rad450_1000,970,962,plot=TRUE) # warning:
#' this is a long computation and may take an hour or more, depending on your
#' system
#' }
#' @export

gridSearch<-function(image,basis,startingRow,startingColumn,plot){
    RSSdtf<-expand.grid(centeri=(startingRow-5):(startingRow+5),centerj=(startingColumn-5):(startingColumn+5))
    RSSdtf$totalRSS<-NA

    for (i in 1:nrow(RSSdtf)){
        shifted<-shiftedImage(image,RSSdtf$centeri[i],RSSdtf$centerj[i])
        allRSS<-unlist(lapply(basis,function(x) statisticsByBasisFunction(x,shifted,c("RSS"))),use.names=FALSE)
        RSSdtf$totalRSS[i]<-sum(allRSS)
    }
    tmp<-RSSdtf[which.min(RSSdtf$totalRSS),]

    multiplier<-2
    while ((tmp$centeri %in% c(min(RSSdtf$centeri),max(RSSdtf$centeri)) | tmp$centerj %in% c(min(RSSdtf$centerj),max(RSSdtf$centerj)))){

        addRSSdtf<-expand.grid(centeri=(startingRow-(5*multiplier)):(startingRow+(5*multiplier)),centerj=(startingColumn-(5*multiplier)):(startingColumn+(5*multiplier)))
        addRSSdtf<-addRSSdtf[abs(addRSSdtf$centeri-startingRow)>(5*(multiplier-1)) | abs(addRSSdtf$centerj-startingColumn)>(5*(multiplier-1)),]
        addRSSdtf$totalRSS<-NA
        for (i in 1:nrow(addRSSdtf)){
            shifted<-shiftedImage(image,addRSSdtf$centeri[i],addRSSdtf$centerj[i])
            allRSS<-unlist(lapply(basis,function(x) statisticsByBasisFunction(x,shifted,c("RSS"))),use.names=FALSE)
            addRSSdtf$totalRSS[i]<-sum(allRSS)
        }
        RSSdtf<-rbind(RSSdtf,addRSSdtf)
        tmp<-RSSdtf[which.min(RSSdtf$totalRSS),]

        multiplier<-multiplier+1
        if (multiplier>10) stop("Center not found after long search")
    }
    if (plot==TRUE){
        surfacePlot(RSSdtf,main=paste0(deparse(substitute(image)),": Total RSS of Basis Functions\nwith Radius of 450-1000 Pixels"))
    }
    return(RSSdtf)
}

#' Surface plot of total residual sums of squares
#'
#' Produces a surface plot of the total residual sum of squares of all the
#' centers tested
#'
#' @param RSS data frame of RSS values, in the format produced by
#'   \code{gridSearch()}
#' @param main title of plot
#'
#' @examples
#' \dontrun{
#' basis1919_rad450_1000<-subsetBasis(basis1919,c(450,1000))
#' RSSdtf<-gridSearch(LL1_3,basis1919_rad450_1000,970,962) # warning: this is
#' a long computation and may take an hour or more, depending on your system
#' surfacePlot(RSSdtf,main="Total RSS of Basis Functions \n
#'              with Radius of 450-1000 Pixels")
#' }
#' @export

surfacePlot<-function(RSS,main){
    par(mar=c(2.1,2.1,5.1,2.1))
    makeGrid<-matrix(NA,nrow=sqrt(nrow(RSS)),ncol=sqrt(nrow(RSS)))
    for (i in 1:nrow(RSS)){
        makeGrid[RSS$centeri[i]-min(RSS$centeri)+1,RSS$centerj[i]-min(RSS$centerj)+1]<-RSS$totalRSS[i]
    }
    minI<-min(RSS$centeri)
    minJ<-min(RSS$centerj)
    maxI<-max(RSS$centeri)
    maxJ<-max(RSS$centerj)
    tmp<-rep(NA,4)
    tmp[1]<-RSS$totalRSS[RSS$centeri==minI & RSS$centerj==minJ]
    tmp[2]<-RSS$totalRSS[RSS$centeri==maxI & RSS$centerj==minJ]
    tmp[3]<-RSS$totalRSS[RSS$centeri==maxI & RSS$centerj==maxJ]
    tmp[4]<-RSS$totalRSS[RSS$centeri==minI & RSS$centerj==maxJ]
    thetaToUse<-switch(EXPR=which.min(tmp),30,120,210,300)

    graphics::persp(z=makeGrid,xlab=paste0("Row: ",min(RSS$centeri)," to ",max(RSS$centeri)),ylab=paste0("Column: ",min(RSS$centerj)," to ",max(RSS$centerj)),zlab="Total RSS",col="lightblue",main=main,theta=thetaToUse)
    mtext(paste0("Smallest RSS at (",RSS$centeri[which.min(RSS$totalRSS)],", ",RSS$centerj[which.min(RSS$totalRSS)],")"))
}


