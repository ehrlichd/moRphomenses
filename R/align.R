#' Array Data
#'
#'
#' Construct a ragged array (containing missing data) of a specified length (up/down sampling individuals to fit).
#'
#' @name mm_ArrayData
#' @param IDs A vector that contains indiviaul IDs repeated for muliple days of collection.
#' @param Days A vector that contains information on time, IE Day 1, Day 2, Day 3. Note: this vector should include integers, continuous data might produce unintended results.
#' @param Value A vector containing the variable sampled.
#' @param Mid A vector containng the midpoint day for each individual. Note: Mid must have the same number of ervations as unique Individuals.
#' @param StartDay Default is starting at Day 1, can specify other values to subsample.
#' @param EndDay If NULL (default), the highest Day is used for each individual.
#' @param avgLength Integer. Number of days to up/down sample ervations to using \code{\link{mm_GetInterval}}.
#' @param avgMid If NULL (default) 0 will be centered using mm_interval.
#'
#' @return Returns a 3D array of data to be analyzed with individuals in the 3rd dimension.
#' @export
#'
mm_ArrayData <-
  function(IDs,
           Days,
           Value,
           Mid,
           StartDay = 1,
           EndDay = NULL,
           avgLength,
           avgMid = NULL) {

    IDs <- as.factor(IDs)
    IDlevs <- levels(IDs)

    if (length(Mid) != length(IDlevs)) {
      stop("length of 'Mid' not equal to number of individuals")
    }
    if (!all(
      length(Days) == length(Value) |
      length(Days) == length(IDs) |
      length(Value == length(IDs))
    )) {
      stop("length of observations not equal")
    }

    aDat <- array(dim = c(avgLength, 2, length(IDlevs)))
    dimnames(aDat)[[3]] <- IDlevs
    mshpx <- mm_GetInterval(days = avgLength, day0 = avgMid)

    dat1 <- cbind(IDs, Days, Value)

    for (i in 1:length(IDlevs)) {
      ## ID, Day, Val
      ss <- dat1[IDs == IDlevs[i], ]
      if (is.null(EndDay)) {
        maxi <- max(ss[, 2], na.rm = T) ## max day
      } else {
        maxi <- EndDay
      }
      seq1 <- StartDay:maxi
      anchi <- Mid[i]

      mati <- cbind(seq1, rep(NA, length(seq1)))

      ## fill in matrix for days with values present
      mati[seq1 %in% ss[, 2], 2] <- ss[, 3]

      ## center values based on anchor day
      ## IE align by day of ovluation (with ovulation being the midpoint)
      seq2 <- seq1 - anchi
      mati[, 1] <- seq2

      lh <- mati[seq2 < 0, ]
      cent <- mati[seq2 == 0, ]
      uh <- mati[seq2 > 0, ]

      ## scaled  days
      lh[, 1] <- lh[, 1] / length(lh[, 1])
      uh[, 1] <- uh[, 1] / length(uh[, 1])

      mat2 <- rbind(lh, cent, uh)

      fill <- matrix(nrow = avgLength, ncol = 2)
      for (j in 1:avgLength) {
        val <- which.min(abs(mat2[, 1] - mshpx[j]))
        fill[j, 1] <- mshpx[j]
        fill[j, 2] <- mat2[val, 2]
      }
      aDat[, , i] <- fill
    }
    return(aDat)
  }


#' Min-Max Scaling
#'
#' Scale a vector from 0,1 based on its minimum and maximum values.
#'
#' @param x A Numeric vector to be scaled. Missing values are allowed and ignored.
#'
#' @return Returns a scaled vector
#'
#' @examples
#' mm_MinMaxScale(1:10)
#' @export
#'

mm_MinMaxScale <- function(x){
  return((x- min(x, na.rm = T)) /(max(x, na.rm = T)-min(x, na.rm = T)))
}

#' Geometric Scaling
#'
#' Calculate the geometric mean of a vector and scale all values by it.
#'
#' @param x A numeric vector to be scaled. Missing values will produce NA, conduct knn imputation using mm_FillMissing first.
#'
#' @examples
#' mm_GeomScale(1:10)
#' @export
#'
#
mm_GeomScale <- function(x){
  return(x/(prod(x)^(1/length(x))))
}




#' Impute Missing Data
#'
#'
#' Fill in a ragged away by nearest neighbor imputation
#' @name mm_FillMissing
#' @param A A ragged array (IE, contains missing cells), presumably constructed with \code{\link{mm_ArrayData}}.
#' @param knn Number of nearest neighbors to draw on for imputation (default = 3).
#' @param scale Type of scaling to implement (or not). Must be one of "none", "MinMax", "Geom", "log10", "logE", "zscore".
#'
#' @export
#'

mm_FillMissing <- function(A,
                           knn = 3,
                           scale = c("none: impute missing based on raw values",
                                     "MinMax: impute missing after applying Min-Max Scaling",
                                     "Geom: impute missing after scaling by geometric mean",
                                     "log10: impute missing after applying log-base 10 transformation",
                                     "logE: impute missing after applying natural-log tranformation",
                                     "zscore: impute missing after scaling by individual z-score")){

  n <- dim(A)[[3]]
  if (is.null(dimnames(A)[[3]])){
    dimnames(A)[[3]] <- paste("Spec",1:n, sep = "")
  }

  missing <- apply(A, 3, anyNA)
  intA <- A

  for (i in 1:n){
    if(anyNA(A[,,i])){
      ps <- is.na(A[,2,i]) ## identify points for interpolation
      tar <- A[!ps,,i] ## target individual

      all <- A[!ps,,!missing] ##all individuals

      rank <- numeric(dim(all)[[3]])
      names(rank) <- dimnames(all)[[3]]
      for (j in 1:dim(all)[[3]]){
        rank[j] <- sum((tar-all[,,j])^2)
      }
      ch <- names(rank[order(rank)])[1:knn] ## chose nunmber of neighbors
      rep <- apply(A[,,ch], c(1,2), mean) ## calculate the replacement values
      intA[ps,,i] <- rep[ps,] ## replace only the missing values
    }
  }


  ## It should be possible to automate this. perhaps with do.call??
  ## logic for this section should actually just be a bunch of standalone if calls to modift intA[,2,i]
  if(scale == "MinMax"){
    for (i in 1:n){
      intA[,2,i] <- mm_MinMaxScale(intA[,2,i])
    }
  }

  if(scale == "Geom"){
    for (i in 1:n){
      intA[,2,i] <- mm_GeomScale(intA[,2,i])
    }
  }

  if(scale == "log10"){
    for (i in 1:n){
      ## figure out log10 tranformation
      # intA[,2,i] <- mm_GeomScale(intA[,2,i])
    }
  }


  if(scale == "logE"){
    for (i in 1:n){
      ## figure out natural log tranformation
      # intA[,2,i] <- mm_GeomScale(intA[,2,i])
    }
  }


  if(scale == "zscore"){
    for (i in 1:n){
      ## figure out z scores
      # intA[,2,i] <- mm_GeomScale(intA[,2,i])
    }
  }



  mshp <- apply(intA, c(1,2), mean)
  outliers <- data.frame("ID" = dimnames(intA)[[3]], "nmis" = numeric(n), "error" = numeric(n))

  for (i in 1:n){
    outliers$nmis[[i]] <- sum(is.na(A[,2,i]))
    outliers$error[[i]] <- sum((intA[,2,i] - mshp[,2])^2)
  }
  out <- list("dat" = intA, "info" = outliers, "grandM" = mshp,  "scaleType" = scale)

  return(out)
}


#' Create equallly spaced intervals.
#'
#'
#' Create a sequence from -1:1 of specified length. Midpoint (day0) can be
#' @name mm_GetInterval
#' @param days The number of days(divisions) fit between -1 and 1 (inclusive)
#' @param day0 If NULL (default), the median integer will be calculated. This produces symmetrical ranges when days = odd number. Can be specified for asymmetric ranges.
#'
#' @examples
#' mm_GetInterval(15) ## Symmetrical sequence from -1 to 1 with 0 in the middle.
#' mm_GetInterval(15, day0 = 8) ## The same sequence, explicitly specifying the midpoint
#'
#' mm_GetInterval(15, day0 = 3) ## 15 divisions with an asymmetric distribution.
#'
#' @export
#'
mm_GetInterval <- function(days, day0 = NULL){
  seq1 <- 1:days
  if (is.null(day0)){
    mid <- as.integer(median(seq1))
  } else {
    mid <- day0
  }
  ## center
  seq2 <- seq1 - mid

  lh <- seq2[seq2 < 0]
  uh <- seq2[seq2 > 0]

  ## scale
  lh2 <- lh/length(lh)
  uh2 <- uh/length(uh)

  seq3 <- c(lh2, 0, uh2)
  return(seq3)
}




