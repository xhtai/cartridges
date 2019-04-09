#' List of basis functions for a 1769 x 1769 image.
#'
#' A dataset containing the pixels in each circularly symmetric basis function.
#' For more details see \url{https://github.com/xhtai/cartridges}. These data
#' are produced by \code{\link{getBasisFunctions}}.
#'
#' @format A list of length 236463. Each list item corresponds to one basis
#'   function, and contains a data frame with the number of rows being the
#'   number of pixels in the basis function, and 3 columns:
#' \describe{
#'   \item{i}{i coordinate of pixel (ranges from 1 to 1769)}
#'   \item{j}{j coordinate of pixel (ranges from 1 to 1769)}
#'   \item{i2j2}{square of distance of pixel from the center}
#' }
#' Basis functions are enumerated from center outwards, i.e. the first basis
#' function is distance 0 from the center, the second is distance 1 from the
#' center, etc.
"basis1769"
