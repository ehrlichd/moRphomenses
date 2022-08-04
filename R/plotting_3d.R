## separete file for 3D plots


#' Calculate elipses by groups
#'
#' Calculate coordinates of 3D elipsoids using \code{rgl}
#'
#' @param dat A 2D matrix, most likely PC scores
#' @param grps A vector of group IDs
#' @param confidence CI to draw ellipses at. Default is 95%
#' @return If successful, returns a mesh3d object (list) to be plotted with rgl functions. See \code{mm_add_3D_ellipses}.
#'
#'
#'
#'
mm_calc_Ellipses <- function(dat, grps, confidence = .95){
  fgrps <- as.factor(grps)
  k <- length(levels(fgrps))
  out <- vector(mode = "list", length = k)
  for ( i in seq_along(out)){
    out[[i]] <- NULL

    out[[i]] <- try(rgl::ellipse3d(x= cov(dat[fgrps==i,]),
                                   centre = colMeans(dat[fgrps==i,]),
                                   level = confidence), silent = TRUE)
  }
  names(out) <- paste0("g",1:k)
  return(out)
}

#' Add ellipses to plot
#'
#' Once mesh3d objects are calculated, add them to the current rgl window
#'
#' mesh3d ellipses must be clalculated first using \code{mm_calc_ellipses()} and the output saved to an object. That object is used for \code{elist}.
#'
#' @param elist Output of \code{mm_calc_ellipses}, or an otherwise defined mesh3d object.
#' @param cols a vector of colors to plot each ellipse
#' @param alpha optional value to specify transparency 0 = invisible, 1 = opqaue.
#'
#'
#'
mm_add_3D_ellipses <- function(elist, cols, alpha){
  l <- length(elist)
  if(length(cols)==1){
    allCols <- rep(cols, l)
  } else {
    allCols <- cols
  }

  if(length(alpha)==1){
    allAlpha <- rep(alpha,l)
  } else {
    allAlpha <- alpha
  }

  for(i in seq_along(elist)){
    if(class(elist[[i]]) != "mesh3d"){
      next
    } else {
      shade3d(elist[[i]],col = allCols[i], alpha = allAlpha[i])
    }
  }
}
