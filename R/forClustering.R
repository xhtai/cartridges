#' Hierarchical clustering given pairwise distances
#'
#' @param allPairwise name of data frame containing all pairwise comparisons.
#'   This needs to have at least three columns, one representing the first item
#'   in the comparison, one representing the second item, and the last
#'   representing a similarity or distance metric. These are enumerated in the
#'   next two parameters.
#' @param pairColNums vector of length 2 indicating the column numbers in
#'   `allPairwise` of 1. item 1 in comparison, 2. item 2 in comparison
#' @param distSimCol name of column in `allPairwise` indicating distances or
#'   similarities, input as character, e.g. "l2dist". If this is a similarity
#'   and not a difference, input `myDist` parameter to be FALSE
#' @param linkage one of "single", "complete", "average", "centroid", "minimax"
#' @param myDist is `distSimCol` a distance or similarity measure? Default TRUE,
#'   i.e. distance measure
#'
#' @return tree (dendrogram) after hierarchical clustering
#' @export
#' @importFrom protoclust protoclust

getHcluster <- function(allPairwise, pairColNums, distSimCol, linkage, myDist = TRUE) {
    distMat <- longToSquare(allPairwise, pairColNums, distSimCol, myDist)
    distObj <- as.dist(distMat)

    if (linkage != "minimax") {
      hcluster <- hclust(distObj, method = linkage)
    } else if (linkage == "minimax") {
      hcluster <- protoclust::protoclust(distObj)
    }
    return(hcluster)
}

#' Generate distance matrix from pairwise list and distance (or similarity)
#' column
#'
#' @param allPairwise name of data frame containing all pairwise comparisons.
#'   This needs to have at least three columns, one representing the first item
#'   in the comparison, one representing the second item, and the last
#'   representing a similarity or distance metric. These are enumerated in the
#'   next two parameters.
#' @param pairColNums vector of length 2 indicating the column numbers in
#'   `allPairwise` of 1. item 1 in comparison, 2. item 2 in comparison
#' @param distSimCol name of column in `allPairwise` indicating distances or
#'   similarities, input as character, e.g. "l2dist". If this is a similarity
#'   and not a difference, input `myDist` parameter to be FALSE
#' @param myDist is `distSimCol` a distance or similarity measure? Default TRUE,
#'   i.e. distance measure
#'
#' @return symmetric distance matrix
#' @export

longToSquare <- function(allPairwise, pairColNums, distSimCol, myDist = TRUE) {
  hashes <- unique(c(allPairwise[, pairColNums[1]], allPairwise[, pairColNums[2]]))
  hashes <- sort(hashes)

  distMat <- matrix(0, nrow = length(hashes), ncol = length(hashes))
  distMat[lower.tri(distMat)] <- 1 # don't actually need this

  if (myDist == FALSE) {
    allPairwise[, distSimCol] <- 1 - allPairwise[, distSimCol]
  }

  for (ii in 1:nrow(allPairwise)) { # filing in lower.tri
    tmpi <- which(hashes == allPairwise[ii, pairColNums[1]])
    tmpj <- which(hashes == allPairwise[ii, pairColNums[2]])
    i <- max(tmpi, tmpj)
    j <- min(tmpi, tmpj)
    distMat[i, j] <- allPairwise[ii, distSimCol]
  }

  # make symmetric (fill in upper.tri)
  distMat[upper.tri(distMat)] <- t(distMat)[upper.tri(distMat)]
  return(distMat)
}

#' Generate binary column representing links between pairs, given dendrogram and
#' cut specification
#'
#' @param allPairwise name of data frame containing all pairwise comparisons.
#'   This needs to have at least two columns, one representing the first item in
#'   the comparison, and one representing the second item.
#' @param pairColNums vector of length 2 indicating the column numbers of 1.
#'   item 1 in comparison, 2. item 2 in comparison
#' @param htree tree (dendrogram) after hierarchical clustering, returned by
#'   `getHcluster`
#' @param numClust number of desired clusters in clustering, can be a number
#'   from 1 to n, where n is the number of items to be clustered
#' @param linkage one of "single", "complete", "average", "centroid", "minimax"
#'
#' @return binary vector representing whether links exist between pairs in
#'   `allPairwise`
#' @export
#' @importFrom protoclust protocut
#' @importFrom dplyr left_join

makeLinkCol <- function(allPairwise, pairColNums, htree, numClust, linkage) {
    if (linkage != "minimax") {
      set.seed(0)
      clustersAll <- cutree(htree, k = numClust)
    } else if (linkage == "minimax") {
      set.seed(0)
      clustersAll <- protoclust::protocut(htree, k = numClust)$cl
    }

    adjacencyMatrix <- clust2Mat(clustersAll)

    ################ need to get it to return the link col
    adjacencyMatrix[lower.tri(adjacencyMatrix)] <- 0
    getHashes <- data.frame(which(adjacencyMatrix == 1, arr.ind = TRUE))

    hashes <- unique(c(allPairwise[, pairColNums[1]], allPairwise[, pairColNums[2]]))
    hashes <- sort(hashes)

    colName1 <- names(allPairwise)[pairColNums[1]] # say "hash1"
    colName2 <- names(allPairwise)[pairColNums[2]]
    getHashes[, colName1] <- hashes[getHashes[, "row"]]
    getHashes[, colName2] <- hashes[getHashes[, "col"]]

    tmp <- dplyr::left_join(allPairwise[, pairColNums], getHashes[, c("row", colName1, colName2)], by = c(colName1, colName2))
    tmp$link <- 0
    tmp$link[!is.na(tmp$row)] <- 1
    return(tmp$link)
}

# from GraphAT package
clust2Mat<-function(memb){
  N<-length(memb)
  return(as.numeric(outer(memb, memb, FUN="=="))-outer(1:N,1:N,"=="))
}


