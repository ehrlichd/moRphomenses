### trying to work out new array function

sub <- DataToScale[DataToScale$id %in% head(unique(DataToScale$id)),]
ind_ovday <- ind_ovday[1:6]

rm(DataToScale, data_filepath, i, sel, subset, unique_id)



end_day <- NULL
for (i in IDlevs){
  end_day[i] <- max(sub$cycleDay[sub$id==i])

}

fixed <- as.data.frame(rbind(1,ind_ovday, end_day))

fixed2 <- list()

for(i in IDlevs){
  fixed2[[i]] <- fixed[,i]
}



ObsIDs <- sub$id
ObsDays <- sub$cycleDay
ObsValue <- sub$E1G..ng.ml._adjSG
ObsFixed <- fixed2

tar_array <- array(dim = c(28,2,1))
tar_fixed <- c(1,16,28)
tar_scl <- c(-1,0,1)



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

    anchi <- ObsFixed[[i]][centered]
    ## first step is to get it to 0


    mati <- cbind(seq1, rep(NA, length(seq1)))
    mati[seq1 %in% ss[, 2], 2] <- ss[, 3]
    seq2 <- seq1 - anchi
    mati[, 1] <- seq2

    ## scale the fixed landmarks
    flmi <- ObsFixed[[i]]

    iscl <- tar_scl/flmi
    iscl


    seqs <- list()
    for(j in 1:(length(flmi)-1)){
      nn <- nrow(mati[flmi[j]:flmi[j+1],])
      seqs[[j]] <- seq(from = tar_scl[j], to = tar_scl[j+1], length.out = nn)

    }

    seqs <- unique(unlist(seqs))
    seqs


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

test_array <- new_array(
  ObsIDs = sub$id,
  ObsDays = sub$cycleDay,
  ObsValue = sub$E1G..ng.ml._adjSG,
  ObsFixed = fixed2,

  tar_array = array(dim = c(28,2,1)),
  tar_fixed = c(1,16,28),
  tar_scl = c(-1,0,1)
  )

  ### test out limiting the range of the cycle using fixed landmarks
