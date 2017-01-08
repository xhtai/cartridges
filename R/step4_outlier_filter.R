# ACKNOWLEDGEMENTS: all the code in this file was adapted from MATLAB code
# provided by Joseph Roth of Michigan State University. The MATLAB code was
# worked on by Joseph Roth and Zach Richardson, as part of their Breech Face
# Ballistics work.

#' Crop borders of an image
#'
#' Given an image with NA values in the borders, crops out these borders to
#' leave the smallest rectangular region containing valid data.
#'
#' @param image 1769 x 1769 matrix
#' @param primer A 1769 x 1769 binary matrix denoting the primer region
#'
#' @return A matrix, smaller than 1769 x 1769, containing only the valid primer
#'   region. Not necessarily square.
#'
#' @examples
#' \dontrun{
#' croppedExample <- cropBorders(removedExample, centeredExample$centeredPrimer)
#' }
#'
#' @export

cropBorders <- function(image, primer) {
    primerRows <- which(Matrix::rowSums(primer != 0) > 0)
    primerCols <- which(Matrix::colSums(primer != 0) > 0)

    cropped <- image[primerRows, primerCols]
    return(cropped)
}

#' Remove outliers and fill them in with NA
#'
#' For each pixel, we consider a local patch consisting all pixels in a 21 x 21
#' pixel matrix, centered at the pixel being considered. Outliers are defined as
#' pixels that have values larger than 3 standard deviations from the mean of
#' the pixel values in its local patch. This is translated from code by Joseph
#' Roth and Zach Richardson.
#'
#' @param input image matrix
#'
#' @return A matrix with non-breechface areas set to 0, and outliers set to NA.
#'
#' @examples
#' \dontrun{
#' # first run this:
#' croppedExample <- cropBorders(removedExample, centeredExample$centeredPrimer)
#' # then
#' outlierNAexample <- outlierRejection(croppedExample)
#' }
#'
#' @export

outlierRejection <- function(input) {

    valid <- !is.na(input) # will give a matrix
    input[is.na(input)] <- 0
    H <- matrix(1, nrow = 21, ncol = 21)
    N <- filterViaFFT(valid, H) # just counts the number of points in neighborhood
    N[N < 1] <- 0 # numerical error
    mu <- filterViaFFT(input, H)/N
    mu[N == 0] <- 0
    #suppressWarnings(sigma <- sqrt(filterViaFFT( (input - mu)^2, H ) / N)) # this is wrong

    firstTerm <- filterViaFFT(input^2, H)/N
    suppressWarnings(sigma <- sqrt(firstTerm - mu^2)) # this line gives a warning:
    sigma[N == 0] <- 0

    input[((input - mu)/sigma) > 3 & valid != 0] <- NA

    return(input)
}

# ACKNOWLEDGEMENTS: inpaint_nans was translated from John D'Errico's MATLAB code
# (https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans,
# Copyright (c) 2009, John D'Errico). This code was optimized for MATLAB and
# this current R implementation is much slower (e.g. on an example image it runs
# almost instantly on MATLAB but takes 6 minutes in R).

#' Fill in NA values
#'
#' This was translated from John D'Errico's MATLAB code:
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans},
#' Copyright (c) 2009, John D'Errico.
#'
#' @param A image matrix
#'
#' @return A matrix with NA values filled in.
#'
#' @examples
#' \dontrun{
#' # first run this:
#' croppedExample <- cropBorders(removedExample, centeredExample$centeredPrimer)
#' outlierNAexample <- outlierRejection(croppedExample)
#' # then
#' inpaintedExample <- inpaint_nans(outlierNAexample)
#' }
#'
#' @export

inpaint_nans <- function(A) {
    originalA <- A
    n <- dim(A)[1]
    m <- dim(A)[2]
    A <- c(A)
    nm <- n*m
    k <- is.na(A)

    # % list the nodes which are known, and which will
    # % be interpolated
    nan_list <- which(k == TRUE)
    known_list <- which(k == FALSE)

    # % how many nans overall
    nan_count <- length(nan_list)

    # % convert NaN indices to (r,c) form
    # % nan_list==find(k) are the unrolled (linear) indices
    # % (row,column) form
    rowcol <- which(is.na(originalA), arr.ind = TRUE)

    # % both forms of index in one array:
    #     % column 1 == unrolled index
    # % column 2 == row index
    # % column 3 == column index
    nan_list <- cbind(nan_list, rowcol)

    # % horizontal and vertical neighbors only
    talks_to <- matrix(c(-1, 0, 0, -1, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
    neighbors_list <- identify_neighbors(n, m, nan_list, talks_to)

    # % list of all nodes we have identified
    all_list <- rbind(nan_list, neighbors_list)

    # % generate sparse array with second partials on row
    # % variable for each element in either list, but only
    # % for those nodes which have a row index > 1 or < n
    L <- which((all_list[, 2] > 1) & (all_list[, 2] < n))
    nl <- length(L)
    if (nl > 0) {
        fda <- Matrix::sparseMatrix(i = rep(all_list[L, 1], 3),
                                    j = rep(all_list[L, 1], 3) + rep(c(-1, 0, 1), each = nl),
                                    x = rep(c(1, -2, 1), each = nl), dims = c(nm, nm))
    } else {
        fda <- Matrix::Matrix(0, nrow = nm, ncol = nm, sparse = TRUE)
    }

    # % 2nd partials on column index
    L <- which((all_list[, 3] > 1) & (all_list[, 3] < m))
    nl <- length(L)
    if (nl > 0) {
        fda <- fda + Matrix::sparseMatrix(i = rep(all_list[L,1], 3),
                                          j = rep(all_list[L, 1], 3) + rep(c(-n, 0, n), each = nl),
                                          x = rep(c(1, -2, 1), each = nl), dims = c(nm, nm))
    }

    # % eliminate knowns
    rhs <- -fda[, known_list] %*% A[known_list]
    k <- which(Matrix::rowSums(fda[, nan_list[, 1], drop = FALSE] != 0) > 0)

    # % and solve...
    B <- A
    B[nan_list[, 1]] <- as(qr.solve(fda[k, nan_list[, 1], drop = FALSE], matrix(rhs[k], ncol = 1)), "vector")

    # % all done, make sure that B is the same shape as
    # % A was when we came in.
    B <- matrix(B, nrow = n, ncol = m)

    return(B)
}

identify_neighbors <- function(n, m, nan_list, talks_to) {

    if (length(nan_list) > 0) {
        # % use the definition of a neighbor in talks_to
        nan_count <- nrow(nan_list)
        talk_count <- nrow(talks_to)

        nn <- matrix(0, nrow = nan_count*talk_count, ncol = 2)
        j <- c(1, nan_count)
        for (i in 1:talk_count) {
            nn[j[1]:j[2], ] <- nan_list[, 2:3] + matrix(rep(talks_to[i, ], nan_count), ncol = 2, byrow = TRUE)
            j <- j + nan_count
        }

        # % drop those nodes which fall outside the bounds of the
        # % original array
        nn <- nn[nn[, 1] %in% 1:n & nn[, 2] %in% 1:m, ]

        # % form the same format 3 column array as nan_list
        neighbors_list <- matrix(c((nn[, 2] - 1)*n + nn[, 1], nn), ncol = 3)
        colnames(neighbors_list) <- c("nan_list", "row", "col")

        # % delete replicates in the neighbors list
        neighbors_list <- unique(neighbors_list)

        # % and delete those which are also in the list of NaNs.
        neighbors_list <- neighbors_list[neighbors_list[, "nan_list"] %in% setdiff(neighbors_list[, 1], nan_list[, 1]), ]
    }

    else {
        neighbors_list <- c()
    }
    return(neighbors_list)
}

#' Filter image.
#'
#' This function applies a low and high-pass Gaussian filter. The low-pass
#' filter uses an 11 x 11 matrix with \code{sigma = 6}, and the high-pass filter
#' uses a 149 x 149 matrix with \code{sigma = 44}. This is translated from code
#' by Joseph Roth and Zach Richardson.
#'
#' @param inputImage image matrix
#' @param nonBF binary matrix of the same dimensions as \code{inputImage}, where
#'   1s indicate the non-breechface region.
#'
#' @return The image matrix after filtering.
#'
#' @examples
#' \dontrun{
#' # first run this:
#' croppedExample <- cropBorders(removedExample, centeredExample$centeredPrimer)
#' # then
#' nonBF <- is.na(croppedExample)
#' processedExample <- gaussianFilter(inpaintedExample, nonBF)
#' }
#'
#' @export

gaussianFilter <- function(inputImage, nonBF) {
    short_F <- EBImage::makeBrush(size = 11, shape = "Gaussian", sigma = 6)
    long_F <- -EBImage::makeBrush(size = 149, shape = "Gaussian", sigma = 44)
    long_F[ceiling(nrow(long_F)/2), ceiling(ncol(long_F)/2)] <- long_F[ceiling(nrow(long_F)/2), ceiling(ncol(long_F)/2)] + 1

    # % Amplifies information at a certain wavelength
    outputImage <- filterViaFFT(inputImage, long_F)
    # % handle boundary of circle case
    tmp_ <- filterViaFFT(nonBF, long_F)
    outputImage <- outputImage + tmp_ * inputImage
    outputImage[nonBF] <- 0

    # % Smooths image to remove noise
    # % TODO: Can try smoothing before amplifying.
    outputImage <- filterViaFFT(outputImage, short_F)
    outputImage[nonBF] <- 0

    return(outputImage)
}



filterViaFFT <- function(A, B) {
    # size of full filter
    m <- dim(A)
    n <- dim(B)
    x <- m + n - 1

    # pad images with 0 so that we do not have circular issues with FFT
    padA <- matrix(0, nrow = x[1], ncol = x[2])
    padB <- matrix(0, nrow = x[1], ncol = x[2])
    padA[1:m[1], 1:m[2]] <- A
    padB[1:n[1], 1:n[2]] <- B

    # Filter in frequency domain
    C <- circshift(fftshift( fft(fft(padA)*Conj(fft(padB)), inverse = TRUE)/(prod(x)) ), round2((n - m)/2, 0))

    half_n <- round2(n/2, 0)
    C <- C[half_n[1]:(half_n[1] + m[1] - 1), half_n[2]:(half_n[2] + m[2] - 1)]
    if (all.equal(c(Im(C)), rep(0, prod(dim(C)))) == FALSE) {
        stop("Non-zero imaginary part")
    }
    return(Re(C))
}

# to copy behavior of round() in matlab
# http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
}

# http://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r
fftshift <- function(input_matrix) {

    rows <- dim(input_matrix)[1]
    cols <- dim(input_matrix)[2]

    swap_up_down <- function(input_matrix) {
        rows_half <- ceiling(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- ceiling(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    input_matrix <- swap_up_down(input_matrix)
    return(swap_left_right(input_matrix))

}

# http://stackoverflow.com/questions/18791212/equivalent-to-numpy-roll-in-r
circshift <- function(x, vec) {
    dimx <- dim(x)
    # row first
    if (vec[1] != 0) {
        #out <- rbind(x[(dimx[1] - vec[1] + 1):dimx[1], ], x[1:(dimx[1] - vec[1]), ])
        tmp <- c(t(x))
        x <- matrix(c( tail(tmp, vec[1]*dimx[2]) , head(tmp, -vec[1]*dimx[2]) ), byrow = TRUE, nrow = dimx[1])
    }
    # col
    if (vec[2] != 0) {
        tmp <- c(x)
        x <- matrix(c( tail(tmp, vec[2]*dimx[1]) , head(tmp, -vec[2]*dimx[1]) ), nrow = dimx[1])
    }
    return(x)
}
