

#' new alignment function
#'
#' hopefully more flexible
#'
#' @name new_array
#' @param ObsIDs A vector that contains indiviaul IDs repeated for muliple days of collection.
#' @param ObsDays A vector that contains information on time, IE Day 1, Day 2, Day 3. Note: this vector should include integers, continuous data might produce unintended results.
#' @param ObsValue A vector containing the variable sampled.
#' @param ObsFixed A vector containing the days to be treated as fixed. Note: ObsMid must have the same number of observations as unique Individuals.
#' @param tar_array A standard target to align individuals. **EXPLAIN**
#' @param tar_fixed A sequence identifying fixed landmarks in the target
#' @param tar_scl A sequence identifying the scale
#'
#' @return Returns a 3D array of data to be analyzed with individuals in the 3rd dimension.



new_array <-
  function (ObsIDs,
            ObsDays,
            ObsValue,
            ObsFixed,
            tar_array,
            tar_fixed,
            tar_scl
            ) {


  ObsIDs <- as.factor(ObsIDs)
  IDlevs <- levels(ObsIDs)

  if (!all(length(ObsDays) == length(ObsValue) | length(ObsDays) ==
           length(ObsIDs) | length(ObsValue == length(ObsIDs)))) {
    stop("length of observations not equal")
  }


  # if(all(tar_scl!=0i)){
  #   ## something to calculate centroid
  # }
  #


  aDat <- array(dim = c(dim(tar_array)[[1]], 2,length(IDlevs)))

  dimnames(aDat)[[3]] <- IDlevs


  mshpx <- mm_GetInterval(days = dim(aDat)[[1]], day0 = tar_fixed[centered])

  dat1 <- cbind(ObsIDs, ObsDays, ObsValue)

  for (i in 1:length(IDlevs)) {
    ss <- dat1[ObsIDs == IDlevs[i], ]


    rr <- range(ObsFixed[[i]])

    seq1 <- rr[1]:rr[2]

    ## would need another if(is.null(centered)){ anchi <- mean(seq1)}

    cent <- median(ObsFixed)
    ## first step is to get it to 0


    mati <- cbind(seq1, rep(NA, length(seq1)))
    mati[seq1 %in% ss[, 2], 2] <- ss[, 3]
    seq2 <- seq1 - anchi
    mati[, 1] <- seq2

    ## scale the fixed landmarks
    flmi <- ObsFixed[[i]]

    iscl <- tar_scl/flmi



    seqs <- list()
    for(j in 1:(length(flmi)-1)){
      nn <- nrow(mati[flmi[j]:flmi[j+1],])
      seqs[[j]] <- seq(from = tar_scl[j], to = tar_scl[j+1], length.out = nn)

    }

    seqs <- unique(unlist(seqs))



    mat2 <- mati
    mat2[,1] <- seqs

    fill <- matrix(nrow = nrow(mat2), ncol = 2)
    thresh <- abs(seqs[1]) - abs(seqs[2])
    thresh <- thresh/2
    for (k in 1:nrow(fill)) {
      ## calulate explicit threshold

      chk <- abs(mat2[,1] - mshpx[k])
      val <- chk[chk < thresh]
      if(length(val)==0){

        fill[k, 1] <- mshpx[k]
        fill[k, 2] <- NA
      } else {

        fill[k, 1] <- mshpx[k]
        fill[k, 2] <- mat2[mat2[,1]==val, 2]
      }

    }
    aDat[, , i] <- fill
  }
  return(aDat)
  }

