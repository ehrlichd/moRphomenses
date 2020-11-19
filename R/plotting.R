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
#'
#' @export
#'
mm_PlotArray <- function(A, MeanShape = TRUE, AllCols = NULL, MeanCol = NULL, type = c("pts", "lines")){
  n <- dim(A)[[3]]
  if(is.null(AllCols)){
    AllCols <- "grey"
  }

  if(length(AllCols)==1){
    AllCols <- rep(AllCols, n)
  }

  if(is.null(MeanCol)){
    MeanCol <- "black"
  }

  y <- range(A[,2,], na.rm = T)

  mshp <- apply(A, c(1,2), mean)
  plot(mshp, type = "n", main = deparse(substitute(A)), ylim = y, xlab = "", ylab = "")

  if(type == "pts"){
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
#' @param k The number of groups for which phenotypes will be drawn. If NULL (default) Diagnostic plots will be drawn to help gauge approprite k
#' @param maxPC Maximum number of PCs to use in analysis.
#'
#' @export

mm_Phenotype <- function(A, k = NULL, maxPC = 10){

  n <- dim(A)[[3]]
  lbl <- character(n)
  if(is.null(dimnames(A)[[3]])){
    lbl <- paste("Spec",1:n, sep ="")
  } else {
    lbl <- dimnames(A)[[3]]
  }

  A1 <- mm_FlattenArray(A)

  PCA <- prcomp(A1)

  hcl <- hclust(dist(PCA$x),method = "ward.D2")
  hcl <- as.dendrogram(hcl)


  if(is.null(k)){

    barplot(summary(PCA)$importance[2,1:maxPC], main = "PC Loadings")
    abline(h = .05, col = "red")
    abline(h = .01, col = "dark red")
    pairs(PCA$x[,1:maxPC])
    plot(PCA$x[,1:2])

    plot(hcl)

    layout(matrix(1:2, ncol = 2))
    mm_ScreePlot(PCA$x[,1:maxPC])
    mm_SilPlot(PCA$x[,1:maxPC])
    layout(matrix(1))

    out <- list(
      "PCA" = PCA,
      "Dendro" = hcl
    )
    return(out)

  }

  grpID <- list()
  grpShapes <- list()
  for(i in 1:length(k)){
    grpID[[i]] <- data.frame(
      "grpID" = dendextend::cutree(hcl,k = k[i]),
      "grpCol" = character()
      )

    names(grpShapes[[i]]) <- names(grpID[[i]]) <- paste("g",k[i], sep="")


    cols <- rainbow(k[i], s = .4, v = 1) ## light tint
    mcols <- rainbow(k[i], s = 1, v = .4) ## dark shade

    for(j in 1:k[i]){
      grpID[[i]][grpID[[i]][,1]==j,2] <- cols[j]
    }

    barplot(summary(PCA$x)$importance[2,1:maxPC], main = "PC Loadings")
    abline(h = .05, col = "red")
    abline(h = .01, col = "dark red")
    pairs(PCA$x[,1:3], col = grpID[[i]][,2])
    plot(PCA$x[,1:2], col = grpID[[i]][,2])


    cols <- grpID[[i]][,2][order.dendrogram(hcl)]
    dendextend::labels_colors(hcl) <- cols
    plot(hcl)

    layout(matrix(1:k[i]))
    for(q in 1:k[i]){
      mm_PlotArray(A=A[,,grpID[[i]][,1]==q],AllCols = cols[i], MeanCol = mcols[i], type = "lines")

      grpShapes[[i]][[q]] <- apply(A[,,grpID[[i]][,1]==q],c(1,2), mean)

    }

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
  return(out)
}
