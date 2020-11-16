#' Construct a ragged array (containing missing data) of a specified length
#'
#'
#' @param ObsIDs A vector that contains indiviaul IDs repeated for muliple days of collection
#' @param ObsDays A vector that contains information on time, IE Day 1, Day 2, Day 3. Note: this vector should include integers, continuous data might produce unintended results.
#' @param ObsValue A vector containing the variable sampled.
#' @param ObsMid A vector containng the midpoint day for each individual. Note: ObsMid must have the same number of observations as unique Individuals.
#' @param startDay Default is starting at ObsDay 1, can specify other values to subsample.
#' @param endDay If NULL (default), the highest ObsDay is used for each individual.
#' @param scaleTo Integer. Number of days to up/down sample observations to using mm_interval.
#' @param scaleToMid If NULL (default) 0 will be centered using mm_interval.
#'
#' @return Returns a 3D array of data to be analyzed. Data array is in the form [StudyDay , Value , Indivual]
#'  @export
#'
#'
mm_arrayDat <- function(ObsIDs, ObsDays, ObsValue, ObsMid, startDay = 1, endDay = NULL, scaleTo, scaleToMid = NULL){
  ObsIDs <- as.factor(ObsIDs)
  IDlevs <- levels(ObsIDs)
  if(length(ObsMid) != length(IDlevs)){
    stop("length of 'ObsMid' not equal to number of individuals")
  }
  if(!all(length(ObsDays) == length(ObsValue) |
          length(ObsDays)== length(ObsIDs) |
          length(ObsValue == length(ObsIDs)))){
    stop("length of observations not equal")
  }

  aDat <- array(dim = c(scaleTo, 2, length(IDlevs)))
  mshpx <- mm_intervals(days = scaleTo, day0 = scaleToMid)

  dat1 <- cbind(ObsIDs, ObsDays, ObsValue)

  for (i in 1:length(IDlevs)){

    ## ID, ObsDay, ObsVal
    ss <- dat1[ObsIDs==IDlevs[i],]
    if (is.null(endDay)){
      maxi <- max(ss[,2], na.rm = T) ## max day
    } else {
      maxi <- endDay
    }
    seq1 <- startDay:maxi
    anchi <- ObsMid[i]

    mati <- cbind(seq1, rep(NA, length(seq1)))

    ## fill in matrix for days with values present
    mati[seq1 %in% ss[,2],2] <- ss[,3]

    ## center values based on anchor day
    seq2 <- seq1 - anchi
    mati[,1] <- seq2

    lh <- mati[seq2 < 0,]
    cent <- mati[seq2 == 0,]
    uh <- mati[seq2 > 0,]

    ## scale fractional days
    lh[,1] <- lh[,1]/length(lh[,1])
    uh[,1] <- uh[,1]/length(uh[,1])

    mat2 <- rbind(lh, cent, uh)

    fill <- matrix(nrow = scaleTo, ncol = 2)
    for (j in 1:scaleTo){
      val <- which.min(abs(mat2[,1] - mshpx[j]))
      fill[j,1] <- mshpx[j]
      fill[j,2] <- mat2[val,2]
    }
    aDat[,,i] <- fill
  }
  return(aDat)
}



#' Create a sequence from -1:1 of specified length
#'
#' @param days The number of days(divisions) fit between -1 and 1 (inclusive)
#' @param day0 If NULL (default), the median integer will be calculated. This produces (nearly) symmetrical ranges. Can be specified for asymmetric ranges.
#'
#'@export
#'
mm_intervals <- function(days, day0 = NULL){
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
