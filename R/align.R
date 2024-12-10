#' Array Data
#'
#'
#' Construct a ragged array (containing missing data) of a specified length (up/down sampling individuals to fit).
#'
#' @name mm_ArrayData
#' @param IDs A vector that contains individual IDs repeated for multiple days of collection.
#' @param DAYS A vector that contains information on time, IE Day 1, Day 2, Day 3. Note: this vector should include integers, continuous data might produce unintended results.
#' @param VALUE A vector containing the variable sampled.
#' @param MID Am optional vector of midpoints to center each individuals profile. These should be unique to each individual and repeated for each observation of DAYS, VALUE, and IDs. If NULL (defualt), data will not be centered on any day.
#' @param targetLENGTH Integer. Number of days to up/down sample observations to using \code{\link{mm_get_interval}}.
#' @param targetMID If NULL (default) data will not be centered and will range from 0 to 1. If specified, data will be centered on 0 ranging from -1 to 1.
#' @param transformation Which (if any) data transformation to apply. Our reccomendation is minmax, but Geometric mean, Zscore, natural log and log10 transformations are available, if desired.
#' @param impute_missing Integer. If not null, number of nearest-neighbors to use to impute missing data (Default = 3).
#'
#' @return Returns a 3D array of data to be analyzed with individuals in the 3rd dimension.
#' @export
#'
mm_ArrayData <-
  function(IDs,
           DAYS,
           VALUE,
           MID=NULL,
           targetLENGTH,
           targetMID = NULL,
           transformation = c("minmax", "geom", "zscore", "log", "log10"),
           impute_missing = 3){

    IDs <- as.factor(IDs)
    IDlevs <- levels(IDs)

    if (!all(
      length(DAYS) == length(VALUE) |
      length(DAYS) == length(IDs) |
      length(VALUE == length(IDs))

    )) {
      stop("length of observations not equal")
    }



    aDat <- array(dim = c(targetLENGTH, 2, length(IDlevs)))
    unsDat <- array(dim = c(targetLENGTH, 2, length(IDlevs)))
    dimnames(aDat)[[3]] <- IDlevs
    dimnames(unsDat)[[3]] <- IDlevs

    sizeDat <- data.frame(
      "size_x" = numeric(length=length(IDlevs)),
      "size_y" = numeric(length=length(IDlevs))
    )



    if(is.null(MID)){
      ## range will be 0 to 1
      mshpx <- seq(from = 0, to = 1, length.out = targetLENGTH)
      dat1 <- data.frame(IDs, DAYS, VALUE)
    } else {
      ## range will be -1 to 1 with 0 at targetMID
      mshpx <- mm_get_interval(days = targetLENGTH, day0 = targetMID)
      dat1 <- data.frame(IDs, DAYS, VALUE, MID)
    }




    ## For each INDIVIDUAL #####

    for (i in 1:length(IDlevs)) {
      ## ID, Day, Val
      ss <- dat1[IDs == IDlevs[i], ]

      full_days <- data.frame(
        "DAYS" = seq(from=1, to = max(ss$DAYS, na.rm = T))
        )

      full_days2 <- merge(full_days, ss, by = "DAYS", all.x = TRUE)

      ss <- full_days2


      if(!is.null(MID)){
      ss_mid <- median(ss$MID,na.rm = T)

      x_centered <- ss$DAYS - ss_mid

      neg_x <- x_centered[x_centered <= 0] ## we WANT to include 0 at this step
      pos_x <- x_centered[x_centered >= 0]

      }


      ## if no scaling

      mms_y <- ss$VALUE
      uns_y <- ss$VALUE
      mms_x <- ss$DAYS
      size_x <- length(ss$DAYS)
      size_y <- max(ss$VALUE, na.rm = T)

      ## apply scaling

      if(transformation=="minmax"){
        mms_y <- mm_transf_minmax(ss$VALUE)
        size_y <- max(ss$VALUE, na.rm = T)

        if(is.null(MID)){
          mms_x <- mm_transf_minmax(ss$DAYS)
        } else {
          scl_neg_x <- mm_transf_minmax(abs(neg_x))*-1
          scl_pos_x <- mm_transf_minmax(pos_x)
          scl_pos_x <- scl_pos_x[-1] ## drop the duplicate 0

          mms_x <- as.numeric(c(scl_neg_x, scl_pos_x))
        }
      }

      if(transformation=="geom"){
        mms_y <- mm_transf_geom(ss$VALUE)
        size_y <- (prod(ss$VALUE)^(1/length(ss$VALUE)))
        if(is.null(MID)){
          mms_x <- mm_transf_geom(ss$DAYS)
        } else {
          scl_neg_x <- mm_transf_geom(abs(neg_x))*-1
          scl_pos_x <- mm_transf_geom(pos_x)
          scl_pos_x <- scl_pos_x[-1] ## drop the duplicate 0


          mms_x <- as.numeric(c(scl_neg_x, scl_pos_x))
        }
      }

      if(transformation=="zscore"){
        mms_y <- mm_transf_zscore(ss$VALUE)
        size_y <- mean(mms_y) ## avg Z score??
        if(is.null(MID)){
          mms_x <- mm_transf_zscore(ss$DAYS)
        } else {
          scl_neg_x <- mm_transf_zscore(abs(neg_x))*-1
          scl_pos_x <- mm_transf_zscore(pos_x)
          scl_pos_x <- scl_pos_x[-1] ## drop the duplicate 0


          mms_x <- as.numeric(c(scl_neg_x, scl_pos_x))
        }
      }

      if(transformation=="log"){
        mms_y <- mm_transf_log(ss$VALUE)
        size_y <- max(log(ss$VALUE), na.rm=T) ## max log value??
        if(is.null(MID)){
          mms_x <- mm_transf_log(ss$DAYS)
        } else {
          scl_neg_x <- mm_transf_log(abs(neg_x))*-1
          scl_pos_x <- mm_transf_log(pos_x)
          scl_pos_x <- scl_pos_x[-1] ## drop the duplicate 0


          mms_x <- as.numeric(c(scl_neg_x, scl_pos_x))
        }
      }

      if(transformation=="log10"){
        mms_y <- mm_transf_log10(ss$VALUE)
        size_y <- max(log10(ss$VALUE), na.rm=T)
        if(is.null(MID)){
          mms_x <- mm_transf_log10(ss$DAYS)
        } else {
          scl_neg_x <- mm_transf_log10(abs(neg_x))*-1
          scl_pos_x <- mm_transf_log10(pos_x)
          scl_pos_x <- scl_pos_x[-1] ## drop the duplicate 0


          mms_x <- as.numeric(c(scl_neg_x, scl_pos_x))
        }
      }



      ## slide x

      shape_mat <- cbind(mms_x, mms_y)
      uns_mat <- cbind(mms_x, uns_y)

      shape_fill <- matrix(nrow = targetLENGTH, ncol = 2)
      uns_fill <- matrix(nrow = targetLENGTH, ncol = 2)

      for (j in 1:targetLENGTH) {

        # threshold <- round((1/targetLENGTH),3) ## this seems to be causing problems
        threshold <- .05 ## previously .036
        sel_min <- mshpx[j]-threshold
        sel_max <- mshpx[j]+threshold

        sel_vals <- shape_mat[shape_mat[,1] > sel_min & shape_mat[,1] < sel_max,1]


        if(!length(sel_vals)==0){
          ## handle multiple values
          sel_val <- sel_vals[which.min(abs(sel_vals - mshpx[j]))]
          shape_fill[j, 2] <- shape_mat[shape_mat[,1]==sel_val,2]
          uns_fill[j,2] <- uns_mat[shape_mat[,1]==sel_val,2]

        } else {
          shape_fill[j, 2] <- NA_real_
          uns_fill[j,2] <- NA_real_
        }

        shape_fill[j, 1] <- mshpx[j]
        uns_fill[j,1] <- mshpx[j]

      }

      aDat[, , i] <- shape_fill
      unsDat[,,i] <- uns_fill
      sizeDat[i,] <- c(size_x, size_y)

      # if(!is.null(MID)){
      #   left_scale <- ss_MID[1]
      #   right_scale <- size_x-ss_MID[1]
      #   to_fill_left <- shape_fill[shape_fill[,1] < 0,1] * left_scale
      #   unsDat[1:length(to_fill_left),1,i] <- to_fill_left
      #   unsDat[to_fill_left:dim(unsDat)[[1]],1,i] <- shape_fill[shape_fill[,1] < 0,1] * left_scale
      # } else {
      #   unsDat[,1,i] <- shape_fill[,1]*size_x
      # }
      # unsDat[,2,i] <- shape_fill[,2]*size_y
      }

      out <- list(
        "Shape_data" = aDat,
        "Size_data" = sizeDat,
        "Unscaled_y" = unsDat,
        "scaleType" = transformation
      )


    if(!is.null(impute_missing)){
      out$shape_data_wNA <- aDat
      knn_dat <- mm_FillMissing(aDat, knn=impute_missing)
      out$Shape_data <- knn_dat$dat


      out$knn_info <- knn_dat$info
    }




    return(out)
  }





#' Impute Missing Data
#'
#'
#' Fill in a ragged away by nearest neighbor imputation
#' @name mm_FillMissing
#' @param A A ragged array (IE, contains missing cells), presumably constructed with \code{\link{mm_ArrayData}}.
#' @param knn Number of nearest neighbors to draw on for imputation (default = 3).
#' @return Returns an array of the same dimensions with all missing data filled.
#'
#' @export
#'

mm_FillMissing <- function(A,
                           knn = 3){

  n <- dim(A)[[3]]
  if (is.null(dimnames(A)[[3]])){
    dimnames(A)[[3]] <- paste("Obs",1:n, sep = "")
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



  mshp <- apply(intA, c(1,2), mean)
  outliers <- data.frame("ID" = dimnames(intA)[[3]], "nmis" = numeric(n), "error" = numeric(n))

  for (i in 1:n){
    outliers$nmis[[i]] <- sum(is.na(A[,2,i]))
    outliers$error[[i]] <- sum((intA[,2,i] - mshp[,2])^2)
  }
  out <- list("dat" = intA, "info" = outliers, "grandM" = mshp)

  return(out)
}







#' Create equallly spaced intervals.
#'
#'
#' Create a sequence from -1:1 of specified length. MIDpoint (day0) can be
#' @name mm_get_interval
#' @param days The length of the sequence to return, inclusive of the endpoints (-1,1)
#' @param day0 If NULL (default), the median integer will be calculated, centering the range on 0. Specifying a value will set 0 to that value, creating asymmetric ranges.
#' @return Returns a numeric vector of specified length, ranging from -1 to 1
#'
#' @examples
#' mm_get_interval(15) ## Symmetrical sequence from -1 to 1 with 0 in the middle.
#' mm_get_interval(15, day0 = 8) ## The same sequence, explicitly specifying the midpoint
#' mm_get_interval(15, day0 = 3) ## 15 divisions with an asymmetric distribution.
#'
#' @export
#'
mm_get_interval <- function(days, day0 = NULL){
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




## scaling functions #####




#' Min-Max Scaling
#'
#' Scale a vector from 0,1 based on its minimum and maximum values.
#'
#' @param x A Numeric vector to be scaled. Missing values are allowed and ignored.
#'
#' @return Returns a scaled vector
#'
#' @examples
#' mm_transf_minmax(1:10)
#' @export
#'

mm_transf_minmax <- function(x){
  return((x- min(x, na.rm = T)) /(max(x, na.rm = T)-min(x, na.rm = T)))
}

#' Geometric Scaling
#'
#' Calculate the geometric mean of a vector and scale all values by it.
#'
#' @param x A numeric vector to be scaled. Missing values will produce NA, conduct knn imputation using mm_FillMissing first.
#'
#' @return Returns a scaled vector
#'
#' @examples
#' mm_transf_geom(1:10)
#' @export
#'
#
mm_transf_geom <- function(x){
  return(x/(prod(x)^(1/length(x))))
}


#' Z scores
#'
#' Calculate and return z-scores given a numeric vector.
#'
#' @param x A numeric vector to be scaled. Missing values will produce NA, conduct knn imputation using mm_FillMissing first.
#'
#' @return Returns a scaled vector
#' @examples
#' mm_transf_zscore(1:10)
#' @export
#'
#
mm_transf_zscore <- function(x){
  return(mean(x, na.rm=T)/sd(x, na.rm=T))
}


#' natural log transform
#'
#' Transform a vector by the natural log.
#'
#' @param x A numeric vector to be scaled. Missing values will produce NA, conduct knn imputation using mm_FillMissing first.
#'
#' @return Returns a scaled vector
#' @examples
#' mm_transf_log(1:10)
#' @export
#'
#
mm_transf_log <- function(x){
  return(log(x))
}


#' Common log transform
#'
#' Transform a vector by the common log (base 10).
#'
#' @param x A numeric vector to be scaled. Missing values will produce NA, conduct knn imputation using mm_FillMissing first.
#'
#' @return Returns a scaled vector
#' @examples
#' mm_transf_log10(1:10)
#' @export
#'
#
mm_transf_log10 <- function(x){
  return(log10(x))
}




