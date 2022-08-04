#' Check Imputation
#'
#' Plot Raw (alinged) data along side by side with imputed data.
#'
#'
#' @name mm_CheckInt
#' @param A1 An aligned array, containng missing data (presumably made with \code{\link{mm_ArrayData}})
#' @param A2 An aligned and imputed array (presumbably made with \code{\link{mm_FillMissing}})
#' @param ObO One-by-One. If TRUE (default), individuals will be plotted one at a time, requiring the user to advance/exit the operation. If FALSE, all plots will be generated at once to be browsed/exported from  the Plot History
#'
#'
#'
#' @export

mm_CheckInt <- function(A1, A2, ObO = TRUE){
  if (!identical(dim(A1), dim(A2))){
    stop("Arrays must match")
  }

  n <- dim(A1)[[3]]

  lbl <- character(n)
  if(is.null(dimnames(A1)[[3]])){
    lbl <- paste("Spec",1:n, sep="")
  } else {
    lbl <- dimnames(A1)[[3]]
  }

  y1 <- range(A1[,2,], na.rm = T)
  y2 <- range(A2[,2,], na.rm = T)

  for (i in 1:n){
    layout(matrix(1:2))
    ps <- is.na(A1[,2,i])
    plot(A1[,1:2,i], main = lbl[i], pch = 16, col = "grey", ylab = "Raw Values", ylim = y1, xlab = "")

    plot(A2[,1:2,i],  pch = 16, col = "grey", ylab = "Imputed Val.", ylim = y2, xlab = "")
    points(A2[ps,,i], col = "red")
    if(ObO){
      readline("Press [ENTER] to Continue, [ESC] to exit")
    }
  }
}


#' Plot Array
#' Plot individuals and optionally mean form
#' @name mm_PlotArray
#'
#' @param A An array to be plotted
#' @param MeanShape Logical. Should the Mean Shape be calculated and plotted
#' @param AllCols Either a single color for all individuals, or a vector specifying colors for each individual. If NULL (default) individuals will be plotted in grey
#' @param MeanCol A single color for the mean shape. If Null (default) mean shape will be plotted in black
#' @param type Should the data be plotted as points or lines.
#' @param lbl A title (main =) for the plot. If NULL (default) the name of the array will be used.
#' @param yr Y-range, in the form c(0,100)
#' @export
#'
#'
#'

mm_PlotArray <- function(A, MeanShape = TRUE, AllCols = NULL, MeanCol = NULL, type = c("points", "lines"), lbl = NULL, yr = NULL){
  ## added lbl and yr arguments
  ## added flexibility for missing data (i think)
  ## can handle groups of n=1

  if(length(dim(A)) != 3){
    n <- 1
  } else {
    n <- dim(A)[[3]]
  }

  if(is.null(AllCols)){
    AllCols <- "grey"
  }

  if(length(AllCols)==1){
    AllCols <- rep(AllCols, n)
  }

  if(is.null(MeanCol)){
    MeanCol <- "black"
  }

  if(is.null(lbl)){
    lbl <- deparse(substitute(A))
  }

  if(is.null(yr)){
    y <- range(A[,2,], na.rm = T)
  } else {
    y <- yr
  }

  if(n == 1){
    mshp <- A

    if(type == "points"){
      plot( mar = c(1,2,1,1),
            mshp, col = MeanCol, cex = 1.5, main = lbl, ylim = y, xlab = "", ylab = "")
    } else {
      plot( mar = c(1,2,1,1),
            mshp, col = MeanCol, type = "l", lwd = 3, main = lbl, ylim = y, xlab = "", ylab = "")
    }

    return()

  } else {
    mshp <- apply(A, c(1,2),
                  function(x){
                    mean(x, na.rm= T)
                  })
  }

  plot(
    mar = c(1,2,1,1),
    mshp, type = "n", main = lbl, ylim = y, xlab = "", ylab = "")

  if(type == "points"){
    for(i in 1:n){
      points(A[,1:2,i], col = AllCols[i])
    }
    if(MeanShape){
      points(mshp, col = MeanCol, cex = 1.5)
    }

  } else {
    for(i in 1:n){
      points(A[,1:2,i], col = AllCols[i], type = "l")
    }
    if(MeanShape){
      points(mshp, col = MeanCol, type = "l", lwd = 3)
    }
  }
}


#' Plot Arrays of groups
#'
#' Attempts to optimally format a grid of arrays by group
#'
#' 4 Groups will plot as a 2x2 grid, while 9 groups plot in a 3x3.
#' Function is experimental
#'
#' @param A an array to be plotted
#' @param grps a vector defining group IDs to subset along the 3rd dimension of the array
#' @export

mm_grps_PlotArray <- function(A, grps){

  k <- length(levels(as.factor(grps)))

  cols <- rainbow(k, s= .4, v = 1)
  mcols <- rainbow(k, s = 1, v = .4)

  ## Set layout for plotting - this needs work; plot margins can be a big issue
  if (k %in% c(1,2,3)){
    layout(matrix(1:k, ncol = k))
  } else if (k %% 3 == 0){
    layout(matrix(1:k, nrow = 3))
  } else if (k %% 2 == 0){
    layout(matrix(1:k, nrow = 2))
  }  else {
    layout(matrix(1:k, ncol = k))
  }


  for(q in 1:k){
    mm_PlotArray(A=A[,,grps==q],AllCols = cols[q], MeanCol = mcols[q], type = "lines", lbl = paste("g",q,sep=""), yr = range(A[,2,],na.rm = T))
  }
  layout(matrix(1))
}




#' Color leves of a dendrogram
#'
#' Specify color order approriately for a dendrogram
#'
#' Leaves of a dendrogram will be re-ordered compared to most input classifiers. This function takes the study-ordered colors and correctly applies them to the dendrogram using \code{dendextend}
#'
#' @param dendro A dendrogram class object
#' @param cols a vector of colors
#' @export

mm_ColorLeaves <- function(dendro, cols){
  xcols <- cols[order.dendrogram(dendro)]
  dendextend::labels_colors(dendro) <- xcols
  return(dendro)
}



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



#' Get Phenotypes
#'
#' Apply PCA and Hierarchical Clustering to infer phenotypes
#'
#' @param A an array to be analyzed. No missing data allowed (see \code{\link{mm_FillMissing}})
#' @param maxPC Maximum number of PCs to use in analysis. The default (10) is entirely arbitrary. 10 may be too many shape variables to consider, or it may be far from enough.
#' @param k The number of groups for which phenotypes will be drawn. If NULL (default) Diagnostic plots will be drawn to help gauge appropriate k

#'
#' @export

mm_Phenotype <- function(A, maxPC = 10, k = NULL){

  n <- dim(A)[[3]]
  lbl <- character(n)
  if(is.null(dimnames(A)[[3]])){
    lbl <- paste("Spec",1:n, sep ="")
  } else {
    lbl <- dimnames(A)[[3]]
  }

  A1 <- mm_FlattenArray(A)

  PCA <- prcomp(A1)

  PCA$shapes <- list()

  for(i in 1:maxPC){

    PCA$shapes[[i]] <- geomorph::shape.predictor(
        A = A,
        x = PCA$x[,i],
        min = min(PCA$x[,i]),
        max = max(PCA$x[,i]))

  }
  names(PCA$shapes) <- paste0("PC", 1:maxPC)

  for(i in seq_along(PCA$shapes)){
    PCA$shapes[[i]]$min <- matrix(PCA$shapes[[i]]$min, nrow = dim(A)[[1]])
    PCA$shapes[[i]]$max <- matrix(PCA$shapes[[i]]$max, nrow = dim(A)[[1]])
  }

  PCA$eigs <- summary(PCA)$importance[2:3,]



  hcl <- hclust(dist(PCA$x[,1:maxPC]),method = "ward.D2")
  hcl <- as.dendrogram(hcl)


  if(is.null(k)){

    layout(matrix(1))
    barplot(summary(PCA)$importance[2,1:maxPC], main = "PC Loadings")
    abline(h = .05, col = "red")
    abline(h = .01, col = "dark red")
    pairs(PCA$x[,1:maxPC])

    layout(matrix(1))
    plot(PCA$x[,1:2])

    layout(matrix(1))
    plot(hcl)

    layout(matrix(1:2, ncol = 2))
    mm_ScreePlot(PCA$x[,1:maxPC])
    mm_SilPlot(PCA$x[,1:maxPC])
    layout(matrix(1))

    out <- list(
      "PCA" = PCA,
      "Dendro" = hcl
    )
    layout(matrix(1))
    return(out)

  }

  grpID <- list()
  grpShapes <- list()

  for(i in 1:length(k)){
    grpID[[i]] <- data.frame(
      "grpID" = dendextend::cutree(hcl,k = k[i]),
      "grpCol" = character(n)
    )
    grpShapes[[i]] <- list()

    grpID[[i]]$grpCol <- as.character(grpID[[i]]$grpCol)

    cols <- rainbow(k[i], s = .4, v = 1) ## light tint
    mcols <- rainbow(k[i], s = 1, v = .4) ## dark shade

    for(j in 1:k[i]){
      grpID[[i]][grpID[[i]]$grpID==j,2] <- cols[j]
    }

    layout(matrix(1))
    barplot(summary(PCA)$importance[2,1:maxPC], main = "PC Loadings")
    abline(h = .05, col = "red")
    abline(h = .01, col = "dark red")
    pairs(PCA$x[,1:3], col = grpID[[i]][,2])
    layout(matrix(1))
    plot(PCA$x[,1:2], col = grpID[[i]][,2])

    layout(matrix(1))
    plot(mm_ColorLeaves(hcl, grpID[[i]][,2]))

    for(q in 1:k[i]){
      grpShapes[[i]][[q]] <- apply(A[,,grpID[[i]]$grpID==q],c(1,2), mean)
    }

    names(grpShapes)[[i]] <- paste("g",k[i], sep="")
    names(grpID)[[i]] <- paste("g",k[i], sep="")

    layout(matrix(1:2, ncol = 2))
    mm_ScreePlot(PCA$x[,1:maxPC])
    mm_SilPlot(PCA$x[,1:maxPC])
    layout(matrix(1))

  }

  out <- list(
    "PCA" = PCA,
    "Dendro" = hcl,
    "Groups" = grpID,
    "Shapes" = grpShapes
  )
  layout(matrix(1))
  return(out)
}


#' Visualize PC axes
#'
#' Plot a scatterplot and vizualize shape change across the X axis.
#'
#' Meant to be a quick diagnostic plot with minimal customization.
#'
#' @param pheno Output of \code{mm_phenotype}, containing a PCA object with PC shapes
#' @param xPC The PC to be plotted on the x axis
#' @param yPC The PC to be plotted on the y axis
#' @param yr The y-xis range, in the format c(0,1)
#' @param title To be used for the plot


mm_pheno_plot <- function(pheno, xPC = 1, yPC = 2, yr = c(0,1), title = ""){

  eig <- summary(pheno$PCA)$importance
  layout(matrix(1:3, ncol =3))

  x_r <- range(pheno$PCA$x[,xPC])
  y_r <- range(pheno$PCA$x[,yPC])

  yax <- seq(from = y_r[1], to = y_r[2], length.out = 5)
  xax <- seq(from = x_r[1], to = x_r[2], length.out = 5)

  x_lab <- paste("PC", xPC, ": ", eig[2,xPC]*100, "% Var")
  y_lab <- paste("PC", yPC, ": ", eig[2,yPC]*100, "% Var")

  plot(pheno$PCA$shapes[[xPC]]$min, col = "cyan", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "min"))
  points(pheno$PCA$shapes$GrandM, col = "grey", type = "l")


  plot(pheno$PCA$x[,c(xPC,yPC)], main = title, xlab = x_lab, ylab = y_lab, as= T, bty = "n", xaxt = "n", yaxt = "n", type = "n")
  abline(h = 0, col = "grey", lty = 2)
  abline(v = 0, col = "grey", lty = 2)

  #text(rep(0, 5), yax, labels = round(yax),cex = 1, col = "red",pos = 4)
  #text(xax, rep(0, 5), labels = round(xax),cex = 1, col = "blue", pos = 1)

  points(pheno$PCA$x[,c(xPC,yPC)])

  plot(pheno$PCA$shapes[[xPC]]$max, col = "blue", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "max"))
  points(pheno$PCA$shapes$GrandM, col = "grey", type = "l")

}


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

