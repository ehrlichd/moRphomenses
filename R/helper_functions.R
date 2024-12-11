


#' Flatten Array
#'
#' Convert a 3D array to 2D matrix suitable for PCA, etc. Note, this function is identical to geomorph::two.d.array, reproduced here for convenience.
#'
#' @param A, an array to be flattened
#' @param sep Separator to be used for column names
#'
#' @return Returns a flattened array
#' @export
#'
mm_FlattenArray <- function(A, sep = "."){
  pxk <- dim(A)[1]*dim(A)[2]
  n <- dim(A)[3]
  tmp <- aperm(A, c(3,2,1))
  dim(tmp) <- c(n,pxk)
  rownames(tmp)<-dimnames(A)[[3]]
  colnames(tmp)<-as.vector(t(outer(dimnames(A)[[1]],
                                   dimnames(A)[[2]],
                                   FUN = paste, sep=sep)))
  return(tmp)
}


#' Scree Plot
#'
#' Plot total within group sum of squares to evalaute clusters
#'
#' @param x Input data for cluster analysis (IE, PCA)
#' @param maxC Maximum clusters to evaluate
#' @param ... Additional arguments to be passed to plot
#' @return No value, produces diagnostic plot.
#' @export
#'

mm_ScreePlot <- function(x, maxC = 15, ...) {
  wss <- 0
  max_i <- maxC
  for (i in 1:max_i) {
    km.model <- kmeans(x, centers = i, nstart = 20)
    wss[i] <- km.model$tot.withinss
  }
  plot(1:max_i, wss, type = "b",
       xlab = "Number of Clusters",
       ylab = "Mean w/i group SSE",
       las = 2, ...)
}


#' Silhouette Width Plot
#'
#' Plot average silhouete widths to evaluate clusters
#' @param x Input data for cluster analysis (IE PCA)
#' @param maxC Maximum clusters to evaluate
#' @param ... additional arguments passed to plot
#' @return No value, produces diagnostic plot.
#' @export
#'
mm_SilPlot <- function(x, maxC=15, ...) {
  sw <- 0
  max_i <- maxC
  for (i in 2:max_i) {
    km.model <- cluster::pam(x, k = i)
    sw[i] <- km.model$silinfo$avg.width
  }
  sw <- sw[-1]
  plot(2:max_i, sw, type = "b",
       xlab = "Number of Clusters",
       ylab = "Avg. sil. wid.", ...)
}




## NOTE: mm_pheno should be broken up into mm_pheno and mm_diagnostics.



#' Distance from Centroid
#'
#' Calculate and plot group distance from centroid (grand mean)
#'
#' @param dat a 2d matrix of data. Presumably PC scores
#' @param grps a vector defining group IDs
#' @param plots Logical. Should distances be plotted as boxplots? If FALSE, distance calculations are still performed
#' @return A list containing individual distances from the sample mean shape. If
#'   `plots=TRUE`, will also visualize results
#' @export
#'

mm_grp_dists <- function(dat, grps, plots  = TRUE){
  fgrps <- as.factor(grps)
  k <- length(levels(fgrps))
  n <- nrow(dat)

  out <- list()
  evals <- vector(mode = "list", length = k)

  for(i in seq_along(evals)){
    ss <- dat[fgrps==i,]
    cent <- colMeans(ss)
    l <- dim(ss)[1]
    evals[[i]] <- numeric(l)

    for(j in 1:l){
      evals[[i]][[j]] <- sqrt(sum((ss[j,]-cent)^2))
    }
  }
  names(evals) <- paste("g",1:k,sep="")

  evals$Grand <- numeric(n)

  ss <- dat[,]
  cent <- colMeans(ss)
  for(i in 1:n){
    evals$Grand[i] <- sqrt(sum((ss[i,]-cent)^2))
  }

  ## reformat for plotting
  bplots <- data.frame("error" = numeric(n*2), "grps" = factor(n*2))
  bb1 <- NULL
  bb2 <- NULL

  ntab <- as.numeric(table(fgrps))
  ntab <- c(ntab, n)

  for (i in seq_along(evals)){
    bb1 <- c(bb1, evals[[i]])
    bb2 <- c(bb2, rep(names(evals)[i], each = ntab[i]))
  }
  bplots[,1] <- bb1
  bplots[,2] <- bb2

  if(plots == TRUE){
    boxplot(bplots$error ~ bplots$grps, col = c(rainbow(k, s = .4), "grey"), xlab = "grps", ylab = "Distance from Centroid",notch = TRUE)
  }
  out$evals <- evals
  out$plotting <- bplots
  return(out)

}




#' Color leves of a dendrogram
#'
#' Specify color order approriately for a dendrogram
#'
#' Leaves of a dendrogram will be re-ordered compared to most input classifiers. This function takes the study-ordered colors and correctly applies them to the dendrogram using \code{dendextend}
#'
#' @param dendro A dendrogram or hclust class object
#' @param cols a vector of colors
#' @return A dendrogram class object with leaves colored as specified.
#' @export

mm_ColorLeaves <- function(dendro, cols){
  if(!class(dendro) %in% c("dendrogram")){
    dendro <- as.dendrogram(dendro)
  }
  xcols <- cols[order.dendrogram(dendro)]
  dendextend::labels_colors(dendro) <- xcols
  return(dendro)
}




#' Take a color and modify it
#'
#' Modify color/transparency using hsv syntax
#'
#' @param cols a vector of colors, eg: "#0066FF"
#' @param s Either a single value or a vector of same length as cols specifying
#'   a new saturation (range 0-1). colors darken to black (0).
#' @param v Either a single value or a vector of same length as cols specifying
#'   a new value (range 0-1). colors lighten to white (0)
#' @param alpha Either a single value or a vector of same length as cols
#'   specifying a transparency value (range 0-1). colors translucent at 0.
#' @return A vector of colors that have been modified in saturation, value, or
#'   alpha
#' @export

mm_mute_cols <- function(cols, s=NULL,v=NULL,alpha=.4){
  tmp_col <- data.frame(t(rgb2hsv(col2rgb(cols))))
  tmp_col$alpha <- alpha
  if(!is.null(s)){
    tmp_col$s <- s
  }
  if(!is.null(v)){
    tmp_col$v <- v
  }

  out_col <- hsv(tmp_col$h,
                 tmp_col$s,
                 tmp_col$v,
                 tmp_col$alpha)

  return(out_col)
}


#'  Add confidence ellipses to an active scatterplot.
#'
#'
#' @param dat A matrix of data to draw an ellipses around.
#' @param ci Percentage of data to capture. Must be one of c(67.5, 90, 95, 99).
#' @param linesCol Border color of the shape.
#' @param fillCol Fill color of the shape.
#' @param smoothness Lower values will look jagged, higher value will make
#'   smoother lines, but may take a long time to plot. Default value is 20.
#' @return No value. Will add an ellipses of a given size to the current plot.
#' @export

mm_ellipse <- function (dat, ci = c(67.5, 90, 95, 99), linesCol = "black",
                        fillCol = "grey", smoothness = 20)
{
  sm <- smoothness
  if (ci == 90) {
    chi.v <- 4.605
  }
  else if (ci == 95) {
    chi.v <- 5.991
  }
  else if (ci == 99) {
    chi.v <- 9.21
  }
  else if (ci == 67.5) {
    chi.v <- 2.25
  }
  else {
    stop("Invalid CI, please choose either 90,95, or 99")
  }
  cov.dat <- cov(dat)
  cent <- t(colMeans(dat))
  tr <- sum(cov.dat[1, 1], cov.dat[2, 2])
  det <- ((cov.dat[1, 1] * cov.dat[2, 2]) - (cov.dat[1, 2] *
                                               cov.dat[2, 1]))
  ei1 <- (tr + ((tr^2) - 4 * det)^0.5)/2
  ei2 <- tr - ei1
  ei.a <- (ei1^0.5) * (chi.v^0.5)
  ei.b <- (ei2^0.5) * (chi.v^0.5)
  th <- atan2((ei1 - cov.dat[1, 1]), cov.dat[1, 2])
  q <- matrix(nrow = 2, ncol = 2)
  q[1, 1] <- cos(th)
  q[2, 2] <- cos(th)
  q[2, 1] <- sin(th)
  q[1, 2] <- sin(th) * -1
  circ <- matrix(nrow = 2 * sm + 1, ncol = 3)
  for (i in 1:dim(circ)[1]) {
    circ[i, 1] <- (i - 1) * (pi/sm)
    circ[i, 2] <- (q[1, 1] * ei.a * cos(circ[i, 1]) + q[1,
                                                        2] * ei.b * sin(circ[i, 1])) + cent[1, 1]
    circ[i, 3] <- (q[2, 1] * ei.a * cos(circ[i, 1]) + q[2,
                                                        2] * ei.b * sin(circ[i, 1])) + cent[1, 2]
  }
  if (is.null(fillCol)) {
    lines(circ[, c(2, 3)], col = linesCol)
  }
  else {
    polygon(circ[, c(2, 3)], col = fillCol)
    lines(circ[, c(2, 3)], col = linesCol, lwd = 1.5)
  }
}



#'  Visualize shape of target coordinates
#'
#' @param A A landmark array used for the pca
#' @param PCA output of prcomp. Should contain $transormation
#' @param target_coords A single set of X,Y coordinates.
#' @param target_PCs Integer identifying which pc to use on the X and Y axis.
#'   Default is c(1,2) for PC1 on x and PC2 on y
#' @return A landmark array representing the hypothetical shape of a given set
#'   of coordinates.
#' @export
#'

mm_coords_to_shape <- function (A, PCA, target_coords, target_PCs = c(1,2)){
  if(!is.numeric(target_coords)){
    target_coords <- as.numeric(target_coords)
  }

  mshp <- apply(A, c(1,2), mean)
  nr <- nrow(mshp)

  if(length(target_coords)==2){
    coords_mat <- matrix(target_coords, ncol = 2, nrow = 1)
  } else {
    coords_mat <- matrix(target_coords)
  }

  new_shape_long <- coords_mat %*% t(PCA$rotation[,target_PCs])
  new_shape <- matrix(new_shape_long, nrow = nr, ncol = 2, byrow = TRUE)
  new_shape <- new_shape + mshp
  new_shape

}



#' Print basic summary
#'
#' @param aln An object created with mm_ArrayData
#' @param grps (Optional) A numeric vector that defines groupings
#' @return A character vector with basic descriptive information, to be used with
#'   [print()]. If `grps=TRUE`, will return a list of character vectors.
#' @export
#'

print_summary <- function(aln, grps = NULL){

    if(is.null(grps)){
      out <- paste(
        "Summary:", "\n\n",
      "n-Individuals:", "\n",
      dim(aln$Shape_data)[[3]], "\n\n",

      "Mean (sd) X", "\n",
      paste0(
        round(mean(aln$Size_data$size_x),2), " +/- (",
        round(sd(aln$Size_data$size_x),2), ")"
      ),
      "\n\n",

      "Mean (sd) Y",
      "\n",
      paste0(
        round(mean(aln$Size_data$size_y),2), " +/- (",
        round(sd(aln$Size_data$size_y),2), ")"
      ),
      "\n\n",

      "Mean (sd) Error (dist to mean)",
      "\n",
      paste(
        round(mean(aln$knn_info$error),2), "+/- (",
        round(sd(aln$knn_info$error),2), ")"
      )
    )

  } else {
    uu_grps <- unique(grps)

    out <- vector(mode = "list", length = length(uu_grps))

    for(ii in seq_along(uu_grps)){

      out[[ii]] <- paste(
        paste0("Group",ii), "\n",
        "n Individuals:", "\n","  ",
        dim(aln$Shape_data[,,grps == uu_grps[ii]])[[3]], "\n",

        "Mean (sd) X", "\n", "  ",
        paste0(
          round(mean(aln$Size_data$size_x[grps == uu_grps[ii]]),2), " +/- (",
          round(sd(aln$Size_data$size_x[grps == uu_grps[ii]]),2), ")"
        ),
        "\n",

        "Mean (sd) Y", "\n", "  ",
        paste0(
          round(mean(
            aln$Size_data$size_y[grps == uu_grps[ii]]),2)," +/- (",
          round(sd(
            aln$Size_data$size_y[grps == uu_grps[ii]]),2), ")"
        ), "\n",

        "Mean (sd) Error (dist to mean)", "\n", "  ",
        paste(
          round(mean(aln$knn_info$error[grps == uu_grps[ii]]),2), "+/- (",
          round(sd(aln$knn_info$error[grps == uu_grps[ii]]),2), ")"
        ), "\n\n\n"
      )

    }

  }

  out

}

