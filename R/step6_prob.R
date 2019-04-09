#' Compute the probability of obtaining a higher correlation by chance
#'
#' Given a reference distribution of correlation values of non-matches, compute
#' the probability of obtaining a higher correlation, under the assumption that
#' the two images are not a match, i.e. not from the same gun.
#'
#' @param corr maximum correlation between two images, taking into account all
#'   translations and many rotation angles, e.g. the output of
#'   \code{calculateCCFmax()}.
#' @param empiricalNull the non-match distribution to compare the value of
#'   \code{corr} against. We suggest using the empirical null distribution based
#'   on all pairwise comparisons of non-matches in a known database.
#'
#' @return the probability of obtaining a higher similarity score by chance.
#' @examples
#' \dontrun{
#' computeProb(.1, rnorm(50, .07, 1))
#' }
#'
#' @export

computeProb <- function(corr, empiricalNull) {
    return(sum(empiricalNull > corr) / length(empiricalNull))
}
