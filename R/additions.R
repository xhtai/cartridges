#' Compute maximum correlation for two images, taking into account a sequence of
#' rotation angles and all integer-valued translations (updated)
#'
#' This is an updated version of calculateCCFmax(). It normalizes the images
#' after rotation (mean 0, variance 1), and uses a slightly different grid
#' search to reduce runtime.
#'
#' @param image1 first image matrix
#' @param image2 second image matrix. The two images do not need to be the same
#'   size.
#' @param pad logical value for whether the smaller image should be zero-padded.
#'   The default is TRUE. Only set this to FALSE if the images are of the same
#'   size and are square.
#'
#' @return A list with four items: 1) maximum correlation taking into account
#'   many possible rotations and translations of the second image, 2) the
#'   corresponding dx, the horizontal translation, 3) dy, the vertical
#'   translation, and 4) theta, the corresponding rotation angle.
#' @examples
#' \dontrun{
#' calculateCCFmaxSearch(processedExample, processedExample2)
#' }
#'
#' @export
calculateCCFmaxSearch <- function(image1, image2, pad = TRUE) {
    image1_small <- EBImage::resize(image1, w = floor(dim(image1)[1]/4), h = floor(dim(image1)[2]/4))
    image2_small <- EBImage::resize(image2, w = floor(dim(image2)[1]/4), h = floor(dim(image2)[2]/4))

    # pad if different sizes
    if (pad == TRUE) {
        dim1 <- dim(image1_small)
        dim2 <- dim(image2_small)
        max_size = max(c(dim1, dim2))
        padded1 <- matrix(0, nrow = max_size, ncol = max_size)
        padded2 <- matrix(0, nrow = max_size, ncol = max_size)

        row_start <- floor((max_size - dim1[1]) / 2)
        col_start <- floor((max_size - dim1[2]) / 2)
        padded1[(row_start + 1):(row_start + dim1[1]), (col_start + 1):(col_start + dim1[2])] <- image1_small

        row_start <- floor((max_size - dim2[1]) / 2)
        col_start <- floor((max_size - dim2[2]) / 2)
        padded2[(row_start + 1):(row_start + dim2[1]), (col_start + 1):(col_start + dim2[2])] <- image2_small

        image1_small <- padded1
        image2_small <- padded2
    }

    thetas <- c(seq(from = -175, to = 180, by = 5), -7.5, -2.5, 2.5, 7.5)
    allResults <- data.frame(thetas = thetas, dx = NA, dy = NA,
        corr = NA)

    for (i in 1:length(thetas)) {
        rotated <- cartridges:::bilinearInterpolation(image2_small,
            thetas[i])
        mean1 <- mean(image1_small[image1_small != 0])
        image1_small[image1_small != 0] <- image1_small[image1_small !=
            0] - mean1
        mean2 <- mean(rotated[rotated != 0])
        rotated[rotated != 0] <- rotated[rotated != 0] - mean2
        fact <- sqrt(sum(image1_small^2))
        image1_small <- image1_small/fact
        fact <- sqrt(sum(rotated^2))
        rotated <- rotated/fact
        out <- cartridges:::comparison(image1_small, rotated)
        allResults[i, ] <- c(thetas[i], out$dx, out$dy, out$corr)
    }
    tmp <- thetas[which.max(allResults$corr)]
    if (tmp %in% seq(from = -10, to = 10, by = 2.5)) {
        fineThetas <- seq(from = tmp - 2, to = tmp + 2, by = 0.5)
    } else {
        fineThetas <- seq(from = tmp - 4, to = tmp + 4, by = 1)
    }
    fineResults <- data.frame(fineThetas = fineThetas, dx = NA,
        dy = NA, corr = NA)
    for (i in 1:length(fineThetas)) {
        rotated <- cartridges:::bilinearInterpolation(image2_small,
            fineThetas[i])
        mean1 <- mean(image1_small[image1_small != 0])
        image1_small[image1_small != 0] <- image1_small[image1_small !=
            0] - mean1
        mean2 <- mean(rotated[rotated != 0])
        rotated[rotated != 0] <- rotated[rotated != 0] - mean2
        out <- cartridges:::comparison(image1_small, rotated)
        fineResults[i, ] <- c(fineThetas[i], out$dx, out$dy,
            out$corr)
    }
    index <- which.max(fineResults$corr)
    ret <- list(corr = fineResults$corr[index], dx = fineResults$dx[index],
        dy = fineResults$dy[index], theta = fineResults$fineThetas[index])
    return(ret)
}

#' Given a data set with all pairwise comparisons including A-B and B-A
#' comparisons, take the larger of the similarity score
#'
#' @param allResults needs to have columns `compare`, `newImage`, and a column
#'   with similarity scores (see next parameter)
#' @param similarityCol name of column with similarity scores, e.g. `"corr"`
#'
#' @return a data frame with one row per comparison, adding a column for the B-A
#'   match and a column `corrMax`
#'
#' @export
#' @importFrom dplyr inner_join

removeDups <- function(allResults, similarityCol) {
    out <- dplyr::inner_join(allResults[, c("compare", "newImage", similarityCol, "match")], allResults[, c("compare", "newImage", similarityCol)], by = c("compare" = "newImage", "newImage" = "compare")) # corr.x and corr.y

    out$corrMax <- NA
    for (i in 1:nrow(out)) {
        tmp <- out[i, c("compare", "newImage")]
        tmp <- sort(tmp)
        out[i, c("compare", "newImage")] <- tmp
        out$corrMax[i] <- max(out[i, c(paste0(similarityCol, ".x"), paste0(similarityCol, ".y"))])
    }
    out <- out[duplicated(out[, c("compare", "newImage", "corrMax")]) == FALSE, ]
    return(out)
}


# from GraphAT package
clust2Mat<-function(memb){
    N<-length(memb)
    return(as.numeric(outer(memb, memb, FUN="=="))-outer(1:N,1:N,"=="))}


#' Hierarchical clustering to generate final clusters
#'
#' Given similarities/predictions for each pair, generate transitive closures
#' using hierarchical clustering, returning pairwise match prediction
#'
#' @param pairs pairs with `predCol` column, which are similarity
#'   scores/predictions (higher = more similar)
#' @param predCol name of column with similarities/predictions, input as
#'   character, e.g. "myPreds"
#' @param myCutoff similarity cutoff (>= myCutoff corresponds to a "match")
#' @param myMethod type of linkage, can be "single", "complete", "average",
#'   "minimax"
#' @param hash1 name of column with first item in comparison, e.g. "image1"
#' @param hash2 name of column with second item in comparison, e.g. "image2"
#' @return vector (logical) indicating if pair is matched or not, length is the
#'   same as `pairs`
#' @export
#' @importFrom protoclust protocut protoclust

linksAnalysis <- function(pairs, predCol, myCutoff, myMethod, hash1, hash2) {
    tmp <- which(pairs[, predCol] >= myCutoff)
    if (length(tmp) == 0) {
        return(rep(0, nrow(pairs)))
    } else {
        forHierarchical <- cbind(pairs[tmp, c(hash1, hash2)], preds = pairs[tmp, predCol])
        # cat("Original number of links: ", nrow(forHierarchical), "\n")
        hashes <- unique(c(forHierarchical[, hash1], forHierarchical[, hash2]))
        hashes <- sort(hashes)

        distMat <- matrix(NA, nrow = length(hashes), ncol = length(hashes))
        distMat[lower.tri(distMat)] <- 1

        for (ii in 1:nrow(forHierarchical)) {
            tmpi <- which(hashes == forHierarchical[ii, hash1])
            tmpj <- which(hashes == forHierarchical[ii, hash2])
            i <- max(tmpi, tmpj)
            j <- min(tmpi, tmpj)
            distMat[i, j] <- 1 - forHierarchical$preds[ii]
        }

        distObj <- as.dist(distMat, diag = FALSE, upper = FALSE)

        if (myMethod != "minimax") {
            hcluster <- hclust(distObj, method = myMethod)
            clustersAll <- cutree(hcluster, h = 1 - myCutoff)
        } else if (myMethod == "minimax") {
            hcluster <- protoclust::protoclust(distObj)
            clustersAll <- protoclust::protocut(hcluster, h = 1 - myCutoff)$cl
        }

        adjacencyMatrix <- cartridges:::clust2Mat(clustersAll)

        ################ need to get it to return the link col
        adjacencyMatrix[lower.tri(adjacencyMatrix)] <- 0
        getHashes <- data.frame(which(adjacencyMatrix == 1, arr.ind = TRUE))

        getHashes[, hash1] <- hashes[getHashes[, "row"]]
        getHashes[, hash2] <- hashes[getHashes[, "col"]]

        tmp <- dplyr::left_join(pairs[, c(hash1, hash2)], getHashes[, c("row", hash1, hash2)], by = c(hash1, hash2))
        tmp$link <- 0
        tmp$link[!is.na(tmp$row)] <- 1
        return(tmp$link)
    }
}


### modified from plot.PRROC from the PRROC package
#' @importFrom DescTools AUC
#' @import PRROC
plotPRROC <- function(x, xlim = c(0, 1), ylim = c(0, 1), auc.main = TRUE,
    # auc.type = c("integral", "davis.goadrich"),
    legend = ifelse(is.logical(color) &
        color == TRUE, 4, NA), xlab = NULL, ylab = NULL, main = NULL,
    color = TRUE, lwd = 3, add = FALSE, scale.color = hsv(h = seq(0,
        1, length = 100) * 0.8, s = 1, v = 1), max.plot = FALSE,
    min.plot = FALSE, rand.plot = FALSE, fill.area = (max.plot &
        min.plot), maxminrand.col = grey(0.5), fill.color = grey(0.95),
    ...) {
    # auc.type <- match.arg(auc.type)
    tmpAUC <- DescTools::AUC(c(1, x[, 1], 0), c(0, x[, 2], 1), na.rm = FALSE)
    min <- 0
    max <- 1
    if (is.null(xlab)) {
        my.xlab <- "Recall"
    }
    else {
        my.xlab <- xlab
    }
    if (is.null(ylab)) {
        my.ylab <- "Precision"
    }
    else {
        my.ylab <- ylab
    }
    if (is.null(main)) {
        my.main <- "PR curve"
    }
    else {
        my.main <- main
    }
    if (auc.main) {
        my.main <- paste(my.main, "\nAUC = ", format(tmpAUC), sep = "", collapse = "")
    }

            cols <- PRROC:::getColor(scale.color, x[, 3], min, max)
            plotscale.color = T
            segment = T

    if (!add & !is.na(legend) & (is.numeric(legend) | suppressWarnings(legend ==
        TRUE)) & plotscale.color) {
        if (is.logical(legend)) {
            legend <- 4
        }
        m <- NULL
        widths <- rep(1, 2)
        heights <- rep(1, 2)
        if (legend == 1) {
            m <- matrix(c(1, 2), nrow = 2)
            heights <- c(4, lcm(2))
        }
        else if (legend == 2) {
            m <- matrix(c(2, 1), nrow = 1)
            widths = c(lcm(2.5), 4)
        }
        else if (legend == 3) {
            m <- matrix(c(2, 1), nrow = 2)
            heights = c(lcm(2), 4)
        }
        else {
            m <- matrix(c(1, 2), nrow = 1)
            widths = c(4, lcm(2.5))
        }
        layout(mat = m, widths = widths, heights = heights)
    }
    if (!add) {
        plot(0, xlim = xlim, ylim = ylim, col = 0, xlab = my.xlab,
            ylab = my.ylab, main = my.main, ...)
    }
    d = nrow(x)
    if (segment) {
        segments(x[1:(d - 1), 1], x[1:(d - 1), 2], x[2:d, 1],
            x[2:d, 2], col = cols, lwd = lwd, ...)
    }
    else {
        lines(x[, 1], x[, 2], col = cols, lwd = lwd, ...)
    }
    if (!add & legend & !is.numeric(color) & color == TRUE) {
        scale <- seq(min, max, length = 100)
        cols <- PRROC:::getColor(scale.color, scale, min, max)
        bak <- par("mar")
        on.exit(par(mar = bak))
        if (legend == 2 | legend == 4) {
            if (legend == 4) {
                par(mar = c(5, 1, 4, 2) + 0.1)
            }
            else {
                par(mar = c(5, 2, 4, 1) + 0.1)
            }
            image(c(1), scale, matrix(scale, nrow = 1), col = cols,
                xlab = "", ylab = "", axes = F)
        }
        else {
            if (legend == 1) {
                par(mar = c(2, 4, 0, 2) + 0.1)
            }
            else {
                par(mar = c(0, 4, 2, 2) + 0.1)
            }
            image(scale, c(1), matrix(scale, ncol = 1), col = cols,
                xlab = "", ylab = "", axes = F)
        }
        axis(legend)
        layout(1)
    }
}


#' Going from links to clusters
#'
#' Given pairwise links as generated by `linksAnalysis()`, produce a list of
#' individual items and their cluster membership
#'
#' @param pairs pairs with `preds` column, a binary indicator for whether the
#'   pair is linked after hierarchical clustering
#' @param preds name of column with link, input as character, e.g. "minimax0.5"
#' @param hash1 name of column with first item in comparison, e.g. "image1"
#' @param hash2 name of column with second item in comparison, e.g. "image2"
#' @param sizes whether to return cluster sizes or not, default TRUE
#' @return data frame with three columns: `image`, the name of the image,
#'   `cluster`, the cluster number the image is a member of, and `clusterSize`,
#'   the size of that cluster
#' @export
#' @importFrom dplyr group_by summarize left_join
#' @importFrom magrittr "%>%"
getClust <- function(pairs, preds, hash1, hash2, sizes = TRUE) {
    # tmp <- which(preds == 1)
    forHierarchical <- pairs[, c(hash1, hash2, preds)]
    hashes <- unique(c(forHierarchical[, hash1], forHierarchical[, hash2]))
    hashes <- sort(hashes)

    distMat <- matrix(NA, nrow = length(hashes), ncol = length(hashes))
    distMat[lower.tri(distMat)] <- 1

    for (ii in 1:nrow(forHierarchical)) {
        tmpi <- which(hashes == forHierarchical[ii, hash1])
        tmpj <- which(hashes == forHierarchical[ii, hash2])
        i <- max(tmpi, tmpj)
        j <- min(tmpi, tmpj)
        distMat[i, j] <- 1 - forHierarchical[ii, preds]
    }

    distObj <- as.dist(distMat, diag = FALSE, upper = FALSE)

    hcluster <- hclust(distObj, method = "single") # doesn't matter because everything is already linked properly
    clustersAll <- cutree(hcluster, h = .5)
    names(clustersAll) <- hashes

    # out <- data.frame(clustersAll) %>% dplyr::group_by(clustersAll) %>% dplyr::summarize(clusterSize = length(clustersAll))
    # tmp <- table(out$clusterSize)
    # largest <- as.numeric(names(tmp)[length(tmp)])

    outClusters <- data.frame(image = names(clustersAll), cluster = clustersAll, stringsAsFactors = FALSE)
    if (sizes == TRUE) {
        rownames(outClusters) <- NULL
        `%>%` <- magrittr::`%>%`
        clusterSizes <- outClusters %>% dplyr::group_by(cluster) %>% dplyr::summarize(clusterSize = length(cluster))
        myClusters <- dplyr::left_join(outClusters, clusterSizes)
        myClusters <- myClusters[order(-myClusters$clusterSize, myClusters$cluster), ]
        return(myClusters)
    } else {
        return(outClusters)
    }
}
