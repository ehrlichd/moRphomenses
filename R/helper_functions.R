


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
#'
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
       ylab = "Within groups sum of squares",
       las = 2, ...)
}


#' Silouhete Width Plot
#'
#' Plot average silhouete widths to evaluate clusters
#' @param x Input data for cluster analysis (IE PCA)
#' @param maxC Maximum clusters to evaluate
#' @param ... additional arguments passed to plot
#'
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
       ylab = "Average silhouette width", ...)
}




## NOTE: mm_pheno should be broken up into mm_pheno and mm_diagnostics.



#' Distance from Centroid
#'
#' Calculate and plot group distance from centroid (grand mean)
#'
#' @param dat a 2d matrix of data. Presumably PC scores
#' @param grps a vector defining group IDs
#' @param plots Logical. Should distances be plotted as vusing boxplots? If FALSE, distance calculations are still performed
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
#' @param s Either a single value or a vector of same length as cols specifying a new saturation (range 0-1). colors darken to black (0)
#' @param v Either a single value or a vector of same length as cols specifying a new value (range 0-1). colors lighten to white (0)
#' @param alpha Either a single value or a vector of same length as cols specifying a transparency value (range 0-1). colors translucent at 0.
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





#' Generate Diagnostic Plots WIP
#'
#' Generate various diagnostic plots
#'
#'
#'
#'


## NOTE: need to figure out what the actual input/output of this should be.

## data array is one main object; how to simplyfy the results of phenotyping??
#
# layout(matrix(1))
# barplot(summary(PCA)$importance[2,1:maxPC], main = "PC Loadings")
# abline(h = .05, col = "red")
# abline(h = .01, col = "dark red")
# pairs(PCA$x[,1:maxPC])
#
# layout(matrix(1))
# plot(PCA$x[,1:2])
#
# layout(matrix(1))
# plot(hcl)
#
# layout(matrix(1:2, ncol = 2))
# mm_ScreePlot(PCA$x[,1:maxPC])
# mm_SilPlot(PCA$x[,1:maxPC])
# layout(matrix(1))





