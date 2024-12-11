
# Explore #####

#' Define Shapespace of aligned dataset.
#'
#' Conduct PCA of shape data and visualize major shape trends.
#'
#' @param dat A 3D array of shape data to be analyzed.
#' @param max_Shapes The maximum amount of PCs to visualize. Default 10.
#' @return A list containing the results of shape-pca, including vizualizations of shape extrema for each Principal Component.
#' @export
#'

mm_CalcShapespace <- function(dat, max_Shapes = 10){

  out <- list()

  ## check if its an array or PCA

  if(class(dat) %in% c("prcomp")){
    warning("input must be a 3D array of shape data.")
  } else if (class(dat) %in% c("matrix", "data.frame", "tibble")){
    warning("Convert input data to landmark array first. See `mm_ArrayData`.")
  } else if (class(dat) %in% c("array")) {

    n <- dim(dat)[[3]]
    lbl <- character(n)
    if(is.null(dimnames(dat)[[3]])){
      lbl <- paste("Obs",1:n, sep ="")
    } else {
      lbl <- dimnames(dat)[[3]]
    }



    A1 <- mm_FlattenArray(dat)

    out$PCA <- prcomp(A1)
    out$PCA$eigs <- summary(out$PCA)$importance[2:3,]


  }


  ## shape
  out$Shapes <- list()


  for(i in 1:max_Shapes){


    out$Shapes[[i]] <- list(
      "min" = mm_coords_to_shape(A=dat, PCA = out$PCA, target_coords = -1, target_PCs = i),
      "max" = mm_coords_to_shape(A=dat, PCA = out$PCA, target_coords = 1, target_PCs = i)
    )

  }
  names(out$Shapes) <- paste0("PC", 1:max_Shapes)

  out$Shapes <- lapply(out$Shapes, function(x){
    lapply(x, function(y){
      class(y) <- "matrix"
      y
    })
  })

  ## add mean shape

  out$Shapes$GrandM <- apply(dat, c(1,2), mean)
  out$ALN <- dat


  class(out) <- "mmPCA"
  return(out)

  ## make the plots




}

#' Run a suite of diagnostic analyses.
#'
#' Conduct a set of analyses to make shape-PCA results easier to interpret.
#' Specifically, this will provide a table of eigen values (optional barplot), provide 5-number summary across each PC, conduct a naive Ward's clustering of PC scores (optional dendrogram, along with silhouette plot and scree plot of individual
#' distance to the sample mean
#'
#'
#'
#' @param dat A 3D array or a mmPCA object (output of mm_CalcShapespace).
#' @param max_PC_viz Maximum number of PCs to include in visualizations (EG Eigenplots, or shape trends.
#' @param max_PC_calc By default (NULL), all PCs will be included in calculations. However, if fewer PCs are required users may specify an integer, n, to get the first n PCS.
#' @param hide_plots By default (FALSE), helpful visuals are plotted.
#' @return Returns a list containing the results of:
#' \itemize{
#'   \item eigs - A table containing individual and cumulutive loadings for each PC
#'   \item PC_5_num - A data.frame containing the fivenum summary for each PC
#'   \item TREE - A dendrogram representing the results of a naive-Ward's clustering
#'   }
#' @export
#'
mm_Diagnostics <- function(dat, max_PC_viz=10, max_PC_calc=NULL, hide_plots = FALSE){

  on.exit(layout(matrix(1)))

  out <- list()
  ## check if array or already PCA

  if(any(class(dat) %in% "mmPCA")){
    out$ALN <- dat$ALN
    out$PCA <- dat$PCA
  } else {
    out$PCA <- prcomp(mm_FlattenArray(dat))
    out$ALN <- dat
  }

  nr <- nrow(out$PCA$x)

  if(is.null(max_PC_calc)){
    max_PC_calc <- ncol(out$PCA$x)
  } else if (!is.integer(max_PC_calc)){
    warning("max_PC_calc must be an integer.")
  }

  ## SHAPESPACE

  ## add eigen values
  out$eigs <- summary(out$PCA)$importance[2:3,1:max_PC_calc]

  ## plot eigenvalues

  if(hide_plots==FALSE){
    layout(matrix(1))
    barplot(out$eigs[1,1:max_PC_viz], main = "PC Loadings")
    abline(h = .05, col = "red")
    abline(h = .01, col = "dark red")
  }



  ## do we want/need pairs plots?

  out$PC_5num <- data.frame(
    "min" = numeric(max_PC_calc),
    "lowerQ" = numeric(max_PC_calc),
    "median" = numeric(max_PC_calc),
    "upperQ" = numeric(max_PC_calc),
    "max" = numeric(max_PC_calc)

  )

  for(ii in seq_along(max_PC_calc)){
    out$PC_5num[ii,] <- fivenum(out$PCA$x[,ii])
  }
  rownames(out$PC_5num) <- paste0("PC_", 1:max_PC_calc, "_summary")


  ## CLSUTERING
  ## add naive Ward's clustering
  out$TREE <- hclust(dist(out$PCA$x[,1:max_PC_calc]),method = "ward.D2")
  ## we probably DO want to go ahead and set this so we can always assume the $tree can have nodes modified.

  out$TREE <- as.dendrogram(out$TREE)

  ## plot dendrogram

  if(hide_plots==FALSE){
    layout(matrix(1))
    plot(mm_ColorLeaves(out$TREE, rep("grey", nr)))
  }

  ## plot scree plot & silhouette plots
  if(hide_plots==FALSE){
    layout(matrix(1:2, ncol = 2))
    mm_ScreePlot(out$PCA$x[,1:max_PC_calc])
    mm_SilPlot(out$PCA$x[,1:max_PC_calc])
    layout(matrix(1))

    }



  class(out) <- "mmDiag"

  return(out)

}





#' Generate Phenotypes
#'
#' Partition sample into clusters, based on information from
#' @param dat Either an Array of shape data, an mmPCA object, or an mmDiag object.
#' @param kgrps A non-negative integer of sub-groups to draw. kgrps=1 will provide results for the whole input dat.
#' @param cuttree_h Optional. Draw clusters by splitting the tree at a given
#' height, h.
#' @param cuttree_k Optional. Draw clsuters by splitting the tree into number of branches, k
#' @param plot_figs Optional. Default = TRUE, plot phenotypes for each set(s) of subgroups.
#' @return If plot_figs=TRUE (Default), plot associated graphs and return a list containing:
#' \itemize{
#'   \item  ALN - an array containing aligned and scaled landmark data, the
#'     output of [mm_ArrayData]
#'   \item PCA - PC scores, eigenvalues, and shape visualizations, the output of
#'     [mm_CalcShapespace]
#'   \item TREE - Dendrogram of PC scores, the output of
#'     [mm_Diagnostics]
#'   \item k_grps - If `kgrps` is specified, a vector defining group membership
#'     (as integer); the results of k-means clustering based on PC scores.
#'   \item cth_grps - If `cth_grps` is specified, a vector defining group
#'     membership (as integer); the results of clustering using
#'     [dendextend::cutree] for a given height.
#'   \item ctk_grps - If `ctk_grps` is specified, a vector defining group
#'     membership (as integer); the results of clustering using
#'   [dendextend::cutree] for a given number of clusters.
#' }
#' @export
#'
mm_Phenotype <- function(dat, kgrps, cuttree_h=NULL, cuttree_k=NULL, plot_figs=TRUE){


  out <- list()
  ## check if array or already PCA

  if(any(class(dat) %in% c("mmPCA", "mmDiag"))){
    out$ALN <- dat$ALN
    out$PCA <- dat$PCA
  } else if(length(dim(dat))==3){
    out$PCA <- prcomp(mm_FlattenArray(dat))
    out$PCA$eigs <- summary(out$PCA)$importance[2:3,]
    out$ALN <- dat
  }

  if(any(class(dat) %in% "mmDiag")){
    out$TREE <- dat$TREE
  } else {

    tt_summary <- mm_Diagnostics(dat,hide_plots = TRUE)
    out$TREE <- tt_summary$TREE
  }


  nr <- nrow(out$PCA$x)

  ## Decide on what if any clusters

  if(!is.null(kgrps)){
    out$k_grps <- kmeans(out$PCA$x, centers = kgrps)
  }

  if(!is.null(cuttree_h)){
    out$cth_grps <- stats::cutree(as.hclust(out$TREE), h = cuttree_h)
  }

  if(!is.null(cuttree_k)){
    out$ctk_grps <- stats::cutree(as.hclust(out$TREE), k = cuttree_k)
  }


  ## plotting

  if(plot_figs){

    ## plot mshp
    mm_PlotArray(out$ALN,lbl = "Full Sample",axis_labels = TRUE)


    if(!is.null(kgrps)){

      uu_kgrps <- unique(out$k_grps$cluster)
      ll_kgrps <- length(uu_kgrps)

      all_cols <- rainbow(ll_kgrps, .4, .8, alpha = .4)
      m_cols <- rainbow(ll_kgrps, .8, .4)

      for(k in uu_kgrps){
        mm_PlotArray(out$ALN[,,out$k_grps$cluster==k],
                     "lbl" = paste0("kgrp_ ", k),
                     MeanCol = m_cols[k],
                     AllCols = all_cols[k],
                     axis_labels = TRUE)
      }
    }

    if(!is.null(cuttree_h)){
      uu_cth_grps <- unique(out$cth_grps)
      ll_cth_grps <- length(uu_cth_grps)

      all_cols <- rainbow(ll_cth_grps, .4, .8, alpha = .4)
      m_cols <- rainbow(ll_cth_grps, .8, .4)

      for(h in uu_cth_grps){
        mm_PlotArray(out$ALN[,,out$cth_grps==h],
                     "lbl" = paste0("hcut_", h),
                     MeanCol = m_cols[h],
                     AllCols = all_cols[h],
                     axis_labels = TRUE
                     )
      }
    }

    if(!is.null(cuttree_k)){
      uu_ctk_grps <- unique(out$ctk_grps)
      ll_ctk_grps <- length(uu_ctk_grps)

      all_cols <- rainbow(ll_ctk_grps, .4, .8, alpha = .4)
      m_cols <- rainbow(ll_ctk_grps, .8, .4)

      for(k in uu_ctk_grps){
        mm_PlotArray(out$ALN[,,out$ctk_grps==k],
                     "lbl" = paste0("kcut_", k),
                     MeanCol = m_cols[k],
                     AllCols = all_cols[k],
                     axis_labels = TRUE)
      }
    }


  }

  return(out)


}

#' Check Imputation
#'
#' Plot Raw (aligned) data along side by side with imputed data.
#'
#'
#' @name mm_CheckImputation
#' @param A1 An aligned array, containing missing data (presumably made with
#'  \code{\link{mm_ArrayData}}`$shape_data_wNA`).
#' @param A2 An aligned and imputed array (presumably made with
#'  \code{\link{mm_ArrayData}}`$Shape_data`).
#' @param ObO One-by-One. If TRUE (default, in interactive sessions), individuals
#'  will be plotted one at a time, requiring the user to advance/exit the
#'  operation. If FALSE, all plots #' will be generated at once to be browsed
#'  or exported from the `Plots` panel.
#' @return A series of plots for each individual in the array. If `ObO=TRUE` user
#' input is required to advance or exit the plotting.
#'
#'
#' @export
#'
mm_CheckImputation <- function(A1, A2, ObO = interactive()){
  if (!identical(dim(A1), dim(A2))){
    stop("Arrays must match")
  }

  on.exit(layout(matrix(1)))

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








# Visualize #####

#' Plot Array
#' Plot individuals and optionally mean form
#' @name mm_PlotArray
#'
#' @param A An array to be plotted
#' @param MeanShape Logical. Should the Mean Shape be calculated and plotted
#' @param AllCols Either a single color for all individuals, or a vector
#'   specifying colors for each individual. If NULL (default) individuals will be plotted in grey
#' @param MeanCol A single color for the mean shape. If Null (default) mean shape
#'   will be plotted in black
#' @param plot_type Should the data be plotted as points or lines.
#' @param lbl A title (main =) for the plot. If NULL (default) the name of the
#'   array will be used.
#' @param yr Y-range, in the form c(0,100)
#' @param axis_labels Should units be printed along the axis. Defaults to FALSE
#'   to maximize the profile shape.
#' @return Plot individual(s) profile(s) in the default graphics device.
#'
#' @export
#'
#'
#'

mm_PlotArray <- function(A,
                         MeanShape = TRUE,
                         AllCols = NULL,
                         MeanCol = NULL,
                         plot_type = c("lines", "points"),
                         lbl = NULL,
                         yr = NULL,
                         axis_labels = FALSE){

  oldpar <- par(no.readonly = TRUE) ## get orig parameters
  on.exit(par(oldpar)) ## ensure old par is reset on exit

  ## default to lines
  if(identical(plot_type, c("lines", "points"))){
    plot_type <- "lines"
  }

  if(length(dim(A)) != 3){
    n <- 1
    if(is.null(yr)){
      y <- range(as.numeric(A[,2]), na.rm = T)
    } else {
      y <- yr
    }
  } else {
    n <- dim(A)[[3]]
    if(is.null(yr)){
      y <- range(as.numeric(A[,2,]), na.rm = T)
    } else {
      y <- yr
    }
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


  if(axis_labels){
    par("mar" = c(2,2,2,1))
  } else {
    par("mar" = c(1,1,2,1))
  }


  if(n == 1){
    mshp <- A

    if(plot_type == "points"){
      plot(
        #mar = c(1,2,1,1),
        mshp, col = MeanCol, cex = 1.5, main = lbl, ylim = y, xlab = "", ylab = "")
    } else {
      plot(
        #mar = c(1,2,1,1),
        mshp, col = MeanCol, type = "l", lwd = 3, main = lbl, ylim = y, xlab = "", ylab = "")
    }


  } else {




  mshp <- apply(A, c(1,2),
                  function(x){
                    mean(x, na.rm= T)
                  })


  plot(
    # mar = c(1,2,1,1),
    mshp, type = "n", main = lbl, ylim = y, xlab = "", ylab = "")

  if(plot_type == "points"){
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
      points(mshp, col = MeanCol, type = "l", lwd = 3, lty = 2)
    }
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
#' @return Returns no values, produces a series of plots.
#'
#' @export

mm_grps_PlotArray <- function(A, grps){

  oldpar <- par(no.readonly = TRUE)
  on.exit(layout(matrix(1)))
  on.exit(par(oldpar),add = TRUE)

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
    mm_PlotArray(A=A[,,grps==q],AllCols = cols[q], MeanCol = mcols[q], plot_type = "lines", lbl = paste("g",q,sep=""), yr = range(A[,2,],na.rm = T))
  }

}




#' Plot Calendar Days
#'

#' Pretty PCA
#'
#' A better PCA plot
#'
#' @param PCA Input data either prcomp or mmPCA.
#' @param xPC The PC to plot on the x axis
#' @param yPC The PC to plot on the y axis
#' @param clas_col A character vector of groupings. Each level will be plotted as a different color.
#' @param legend_cex A scaling factor to be applied specifically to the legend. Set to NULL for scatterplot only.
#'
#' @return Returns no object, plots results of PCA
#' @export
#'
#'
#'
mm_pretty_pca <- function(PCA, xPC=1, yPC=2, clas_col = NULL, legend_cex = .8) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(layout(matrix(1)))
  on.exit(par(oldpar),add = TRUE)

  out <- list()
  if (!class(PCA) %in% c("prcomp", "mmPCA")) {
    stop(paste(deparse(substitute(PCA)), "must be a PCA object"))
  } else if(class(PCA) %in% "mmPCA"){
    out$PCA <- PCA$PCA
  } else {
    out$PCA <- PCA
    out$PCA$eigs <- summary(out$shpPCA)$importance[2:3,]
  }



  x_PC <- out$PCA$x[, xPC]
  y_PC <- out$PCA$x[, yPC]

  n_ind <- nrow(out$PCA$x)

  eigs <- out$PCA$eigs

  sub_lbl <- paste(
    paste0("X: PC ", xPC, "; ", round(eigs[1,xPC],2), "% Variance"),
    paste0("Y: PC ", yPC, "; ", round(eigs[1,yPC],2), "% Variance"), sep = "\n")

  ## Get some values
  x_range <- range(x_PC)
  x_diff <- ceiling(max(abs(x_range)))

  # x_diff <- diff(x_range)/2

  y_range <- range(y_PC)
  y_diff <- ceiling(max(abs(y_range)))
  if(y_diff > x_diff){
    diff <- y_diff
  } else {
    diff <- x_diff
  }
  mean_x <- round(mean(x_PC),2)*100
  mean_y <- round(mean(y_PC),2)*100

  seq <- seq(from = mean_x - diff,
             to = mean_x + diff,
             length.out = 5)

  xax <- data.frame("x" = seq,
                    "y" = rep(0, 5),
                    "labs" = seq)

  yax <- data.frame("x" = rep(0, 5),
                    "y" = seq,
                    "labs" = seq)

  if (is.null(clas_col)) {
    clas_col <- 1

    lbl <- "tmp"

  } else {
    lbl <- deparse(substitute(clas_col))

  }

  if (!inherits(clas_col, "factor")) {
    clas_col <- factor(clas_col)
  }


  if (length(levels(clas_col)) == 1) {

    m_lbl <- paste("PCA of", deparse(substitute(PCA)))

    all_cols <- data.frame("class" = rep(clas_col, n_ind),
                           "col" = rep(hsv(1, 0, .5, .5), n_ind))
    m_pts <- c(mean_x, mean_y)
    m_cols <- "black"
    ggrps <- "No classifier selected"
    cc <- "grey"
  } else {


    m_lbl <- paste(paste("PCA of", deparse(substitute(PCA))), paste("Classified by", lbl, sep = "\n"), sep = "\n")

    ggrps <- levels(clas_col)
    all_cols <- data.frame("class" = clas_col,
                           "col" = character(n_ind))
    cc <- rainbow(length(ggrps),
                  s = .8, v = .6, alpha = .5)
    m_cols <- rainbow(length(ggrps),
                      s = .6, v = .8,
                      alpha = 1)
    m_pts <- matrix(nrow = length(ggrps), ncol = 2)
    for (i in seq_along(ggrps)) {
      all_cols$col[as.integer(all_cols$class) == i] <- cc[i]
      m_pts[i, 1] <- mean(x_PC[as.integer(all_cols$class) == i])
      m_pts[i, 2] <- mean(y_PC[as.integer(all_cols$class) == i])
    }


  }
  par(mar = c(.5,.5,.5,.5))

  if(!is.null(legend_cex)){

  ## Make an empty column and plot legend
  layout(matrix(c(1,2,2), ncol = 3))
  plot(1:5,
       main = "",
       xaxt = "n",
       yaxt = "n",
       bty = "n",
       xlab = "",
       ylab = "",
       type = "n"
  )
  legend("center", legend = ggrps, pch = 22, pt.bg = cc, cex = legend_cex)

  text(1,5, labels = m_lbl, adj = c(0,1))
  text(1,1.2, labels = sub_lbl, adj = c(0,1), cex = .8)

  }


  ## plot individuals without box or annotation

  plot(
    x_PC,
    y_PC,
    pch = 16,
    cex = 2,
    col = all_cols$col,
    main = "",
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    xlab = "",
    ylab = "",
    xlim = range(seq),
    ylim = range(seq)
  )


  ## plot Axes and labels
  points(
    xax,
    type = "l",
    lty = 2,
    lwd = 1,
    col = "black"
  )
  points(
    yax,
    type = "l",
    lty = 2,
    lwd = 1,
    col = "black"
  )
  text(xax,
       labels = xax$labs,
       cex = .8,
       adj = c(0, 1))
  text(yax,
       labels = yax$labs,
       cex = .8,
       adj = c(0, 1))

  ## plot means

  if (length(levels(clas_col)) == 1) {
    points(m_pts[1],
           m_pts[2],
           pch = 23,
           bg = m_cols,
           cex = 1.5)

  } else {
    points(m_pts,
           pch = 23,
           bg = m_cols,
           cex = 1.5)

  }


}





#' Visualize PC axes
#'
#' Plot a scatterplot and vizualize shape change across the X axis.
#'
#' Meant to be a quick diagnostic plot with minimal customization.
#'
#' @param mmPCA Output of \code{mm_CalcShapespace}, containing a PCA object with
#'  PC shapes
#' @param xPC The PC to be plotted on the x axis. If yPC is left null, a
#'   univariate density distribution will be plotted with min/max shapes.
#' @param yPC The PC to be plotted on the y axis.
#' @param yr The y-xis range, in the format c(0,1)
#' @param cols A vector of colors of length n, for use in scatterplot.
#' @param title To be used for the plot
#' @param png_dir A file path to a directory in which to save out PNG figures.
#'   Names will be automatically assigned based on input PC(s).
#' @return Produces a series of plots to visualize PCA analysis. If `png_dir` is
#'   specified, function will save out `.png` files. Otherwise plots will be
#'   displayed in the default plot window.
#' @export



mm_VizShapespace <- function(mmPCA, xPC = 1, yPC = 2, yr = c(0,1.1), cols = NULL, title = "", png_dir = NULL){


  oldpar <- par(no.readonly = TRUE)
  on.exit(layout(matrix(1)))
  on.exit(par(oldpar),add = TRUE)


  if(!any(class(mmPCA) %in% "mmPCA")){
    warnings("First define shapespcae with mm_shapepsace function")
  }


  x_r <- range(mmPCA$PCA$x[,xPC])
  xax <- seq(from = x_r[1], to = x_r[2], length.out = 5)
  x_lab <- paste("PC", xPC, ": ", mmPCA$PCA$eigs[1,xPC]*100, "% Var")

  ## Univariate plot
  if(is.null(yPC)){

    if(!is.null(png)){
      out_path <- file.path(png_dir, paste0("shape_trends_PC_",xPC,"\\_.png"))
      png(filename = out_path, height = 4, width = 10, units = "in",res = 300,family = "cairo")
    }

    layout(matrix(c(2,2,2,2,1,1,3,3,3,3,
                    2,2,2,2,1,1,3,3,3,3), ncol =10, byrow = TRUE))


    ## Plot scatterplot first,

    par("mar" = c(.5,.5,.5,.5))

    plot(density(mmPCA$PCA$x[,xPC]), main = title, xlab = x_lab, ylab = y_lab, as= T, bty = "n", xaxt = "n", yaxt = "n", type = "n")
    abline(h = 0, col = "grey", lty = 2)
    abline(v = 0, col = "grey", lty = 2)

    #text(rep(0, 5), yax, labels = round(yax),cex = 1, col = "red",pos = 4)
    #text(xax, rep(0, 5), labels = round(xax),cex = 1, col = "blue", pos = 1)



    ## Plot Min shape

    par("mar" = c(2,2,2,1))
    plot(mmPCA$Shapes[[xPC]]$min, col = "cyan", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "min"))
    points(mmPCA$Shapes$GrandM, col = "black", type = "l", lty = 2)


    ## Plot Max Shape

    plot(mmPCA$Shapes[[xPC]]$max, col = "blue", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "max"))
    points(mmPCA$Shapes$GrandM, col = "black", type = "l", lty = 2)



  } else {


    if(!is.null(png_dir)){
      out_path <- file.path(png_dir, paste0("shape_trends_PC_",xPC,"_by_",yPC,"\\_.png"))
      png(filename = out_path, height = 12, width = 10, units = "in",res = 300,family = "cairo")
    }

    ## Actually, PCs should be square.

    # y_r <- range(mmPCA$PCA$x[,yPC])
    # yax <- seq(from = y_r[1], to = y_r[2], length.out = 5)

    yax <- seq(from = x_r[1], to = x_r[2], length.out = 5)
    y_lab <- paste("PC", yPC, ": ", mmPCA$PCA$eigs[1,yPC]*100, "% Var")



    layout(matrix(c(
                    0,0,0,0, 5,5,5,5, 0,0,0,0,
                    0,0,0,0, 5,5,5,5, 0,0,0,0,

                    0,0,0,0, 3,3,3,3, 0,0,0,0,
                    1,1,1,1, 3,3,3,3, 4,4,4,4,
                    1,1,1,1, 3,3,3,3, 4,4,4,4,
                    0,0,0,0, 3,3,3,3, 0,0,0,0,

                    0,0,0,0, 2,2,2,2, 0,6,6,6,
                    0,0,0,0, 2,2,2,2, 0,6,6,6
    ), byrow = T, ncol=12, nrow = 8))



    par("mar" = c(2,2,2,1))

    ## xpc min
    plot(mmPCA$Shapes[[xPC]]$min, col = "cyan", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "min"))
    points(mmPCA$Shapes$GrandM, col = "black", type = "l", lty = 2)

    ## ypc min
    plot(mmPCA$Shapes[[yPC]]$min, col = "cyan", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", yPC, "min"))
    points(mmPCA$Shapes$GrandM, col = "black", type = "l", lty = 2)




    ## scatterplot
    par("mar" = c(.5,.5,.5,.5))

    plot(mmPCA$PCA$x[,c(xPC, yPC)], main = title, xlab = x_lab, ylab = y_lab, as= T, bty = "n", xaxt = "n", yaxt = "n", type = "n", xlim = x_r, ylim = x_r, pch = 21, bg = cols)
    abline(h = 0, col = "grey", lty = 2)
    abline(v = 0, col = "grey", lty = 2)




    ##### text(rep(0, 5), yax, labels = round(yax),cex = 1, col = "red",pos = 4)
    ##### text(xax, rep(0, 5), labels = round(xax),cex = 1, col = "blue", pos = 1)

    #### add the points
    points(mmPCA$PCA$x[,c(xPC,yPC)], pch = 16, col = mm_mute_cols("black"))

    ## xpc max

    par("mar" = c(2,2,2,1))

    plot(mmPCA$Shapes[[xPC]]$max, col = "blue", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "max"))
    points(mmPCA$Shapes$GrandM, col = "black", type = "l", lty = 2)

    ## yPC max
    plot(mmPCA$Shapes[[yPC]]$max, col = "blue", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", yPC, "max"))
    points(mmPCA$Shapes$GrandM, col = "black", type = "l", lty = 2)
  }


  ## PLOT Eigs

  par("mar" = c(2,2,2,1))
  barplot(mmPCA$PCA$eigs[1,1:15], ylim = c(0, .25), main = "PC Loadings")
  barplot(mmPCA$PCA$eigs[1,c(xPC, yPC)], col = "red", add= T)
  abline(h = .25, col = "grey", lty = 2)
  abline(h = .05, col = "dark grey", lty = 2)



  if(!is.null(png_dir)){
    dev.off()
  }

}






# Test #####

#' Build, implement, visualize multivariate linear model.
#'
#'Easily evaluate simple model sets (one covariate with up to 2 additional
#'  classifiers/covariates). Helpful for exploratory analysis. For detailed
#'  models or specific combinations of variables, see [geomorph::procD.lm]
#'  for full use of this function.
#'
#' @param shape_data This will be the (multivariate) response variable
#' @param ... Covariate(s)/classifier(s) to build a model set. Individual models
#' are run with interaction effects.
#' @param subgrps Optional. Vector of group membership. Model sets will be run
#' across the whole sample and subgroups. If k is specified, only the full
#' model will be run.
#' @param ff1 An explicit model to test in the format: " coords ~ ...". Names
#' must match those specifed in `...`. Standard lm notation applies.
#' @param univ_series Default (FALSE) will evaluate multiple covariates and their
#'  interaction in a single model. However, it can be helpful to understand the
#'  univariate effects in isolation of interaction/confounding factors. Set
#'  `univ_series=TRUE` to produce a series of model sets, one for each covariate
#'   specified (NOTE: ff1 must also be NULL for this to work).
#'
#' @return A list containing output of one or more multivariate linear models
#'   that can be inspected on their own or interacted with using [mm_VizModel]
#'   or [mm_CompModel].
#'
#' @export


mm_BuildModel <- function(shape_data, ..., subgrps= NULL, ff1 = NULL, univ_series=FALSE){

  input_covs <- list(...)

  chk_lengths <- lapply(input_covs, length)

  nr <- dim(shape_data)[[3]]
  if(!all(chk_lengths==nr)){
    warning("Length of Covariate(s) does not match length of shape data.")
  }

  gdf0 <- geomorph::geomorph.data.frame("coords" = shape_data,  ...)

  all_names <- names(gdf0)

  covariates <- setdiff(all_names, "coords")

  model_sets <- list()

  ## Decide on model
  if(is.null(ff1)){
    ## Exploratory, attempt to test all combinations of covariates
    mod_to_test <- paste("coords ~ ", paste(covariates, collapse = " * "), collapse = "")

    } else {
      ## EXPLICIT: Test a single explicit model once you know best fit
      mod_to_test <- ff1
  }



  ## SPECIFIC model for GROUPS and WHOLE SAMPLE
  if(!is.null(subgrps)){
    ll <- length(levels(as.factor(subgrps)))

    for(gg in 1:ll){
        ## subset data
        sub_gdf <- geomorph::geomorph.data.frame("coords" = gdf0$coords[,,subgrps==gg])
        for(vv in 2:length(gdf0)){
          sub_gdf[[vv]] = gdf0[[vv]][subgrps==gg]
        }
        names(sub_gdf) <- names(gdf0)


        ## apply the model
        model_sets[[gg]] <- list("ff" = mod_to_test)
        model_sets[[gg]]$mod <- geomorph::procD.lm(as.formula(model_sets[[gg]]$ff), data = sub_gdf)
        model_sets[[gg]]$summary <- summary(model_sets[[gg]]$mod)
      }
      names(model_sets) <- paste0("mvlm_gg",1:ll)

      model_sets$mvlm_all <- list("ff" = mod_to_test)
      model_sets$mvlm_all$mod <- geomorph::procD.lm(as.formula(model_sets$mvlm_all$ff), data = gdf0)
      model_sets$mvlm_all$summary <- summary(model_sets$mvlm_all$mod)

      return(model_sets)

    }

  ## SPECIFC model for whole sample ONLY
  if(!is.null(ff1)){
    model_sets$mvlm_all <- list("ff" = mod_to_test)
    model_sets$mvlm_all$mod <- geomorph::procD.lm(as.formula(model_sets$mvlm_all$ff), data = gdf0)
    model_sets$mvlm_all$summary <- summary(model_sets$mvlm_all$mod)

    return(model_sets)
  }


  ## EXPLORE COVARIATES AS A UNIVARIATE SERIES
  if(univ_series){
    for (nn in 1:length(covariates)){
      model_sets[[nn]] <- list("ff" = paste("coords ~ ", covariates[nn], collapse = ""))
      model_sets[[nn]]$mod <- geomorph::procD.lm(as.formula(model_sets[[nn]]$ff), data = gdf0)
      model_sets[[nn]]$summary <- summary(model_sets[[nn]]$mod)
    }

    names(model_sets) <- paste0("mvlm_co", 1:length(covariates))
    return(model_sets)

  }



  ## EXPLORE COVARIATES and INTERACTIONS across whole sample


  for (nn in 1:length(covariates)){
    model_sets[[nn]] <- list("ff" = paste("coords ~ ", paste(covariates[1:nn], collapse = " * "), collapse = ""))
    model_sets[[nn]]$mod <- geomorph::procD.lm(as.formula(model_sets[[nn]]$ff), data = gdf0)
    model_sets[[nn]]$summary <- summary(model_sets[[nn]]$mod)


  }

  names(model_sets) <- paste0("mvlm", 1:length(covariates))

  return(model_sets)



}



#' Visualize Multivariate LM
#'
#' Visualize 2D scatterplot of mvlm including predicted shapes.
#'
#' @param dat Input mvlm, created by [mm_BuildModel] (or by using
#'   [geomorph::procD.lm])
#' @param clas_col A classifier to color the data by. If null (default) all
#'   points will be grey. Otherwise, data will be plotted as rainbow(n) colors.
#' @return A list containing the results of the mvlm, visualizations of shape
#'   trends along the regression line, and the model itself.
#' @export


mm_VizModel <- function(dat, clas_col = NULL){

  oldpar <- par(no.readonly = TRUE)
  on.exit(layout(matrix(1)))
  on.exit(par(oldpar),add = TRUE)

  ## f1 a formula
  ## coords a [p x k x n] matrix of landmark data
  ## covaraites continuous variables to model
  ## classifiers ranked/categorical/binary variables to model
  ## ... a formula to pass to

  out <- list()

  n_ind <- nrow(dat$data)
  png("_mm_tmp__.png")
  ratX <- plot(dat,
               type = "regression",
               predictor = dat$data[[2]],
               reg.type = "RegScore")
  dev.off()
  file.remove("_mm_tmp__.png")

  predsX <- geomorph::shape.predictor(dat$GM$fitted,
                            x = ratX$RegScore,
                            min = min(ratX$RegScore),
                            max = max(ratX$RegScore)
  )

  predsX <- lapply(predsX, function(x){

    nr <- nrow(x)
    nc <- ncol(x)

    x <- matrix(x, nrow=nr, ncol = nc)
  })

  if(is.null(clas_col)){

    clas_col <- hsv(1,.01,.6, alpha = .6)
    all_cols <- data.frame("col" = rep(clas_col, ))

  } else {
    ## make colors better.
    ll <- length(levels(as.factor(clas_col)))


    all_cols <- data.frame("class" = clas_col,
                           "col" = character(n_ind))

    cc <- rainbow(ll,
                  s = .8, v = .6, alpha = .5)
    for (i in seq(ll)) {
      all_cols$col[as.integer(all_cols$class) == i] <- cc[i]

    }
  }

  names(predsX) <- c("min", "max")



  ## plot regMin, scatterplot, regMax

  layout(matrix(c(0,0,0,0,2,2,2,2,2,2,3,3,3,3,
                  0,0,0,0,2,2,2,2,2,2,3,3,3,3,
                  0,0,0,0,2,2,2,2,2,2,0,0,0,0,
                  0,0,0,0,2,2,2,2,2,2,0,0,0,0,
                  1,1,1,1,2,2,2,2,2,2,0,0,0,0,
                  1,1,1,1,2,2,2,2,2,2,0,0,0,0), ncol = 14, byrow = TRUE))

  yset <- c(0,1)

  xlabmod <- names(dat$data)[2] ## this should always be the primary response

  par("mar" = c(2,2,1,1))
  plot(predsX$min, type = "l", lwd = 3, col = hsv(.6,.6,1), xlab = "", ylab = "", ylim = c(0,1))
  plot(ratX$plot_args, pch = 16, col = all_cols$col, xlab = xlabmod, ylab = "RegScore", cex = 3)


  ## add best fit line
  tt_line <- lm(ratX$RegScore ~ ratX$plot_args$x)
  abline(tt_line, col = "black", lty = 2)


  plot(predsX$max, type = "l", lwd = 3, col = hsv(1,.6,1), xlab = "", ylab = "", ylim= c(0,1))



  out$mvlm <- dat
  out$summary <- summary(dat)
  out$Shapes <- list("RegMin" = predsX$min,
                     "RegMax" = predsX$max)


  return(out)
}



#' Compare Complex Model Metrics
#'
#' Compare key figs (Rsq, p-value, etc) across multiple complex models.
#'
#' @param mv_results Input mvlm, created by [mm_BuildModel] (or by using
#'   [geomorph::procD.lm])
#' @param row_labels A character vector to use in output. If NULL (default)
#'   labels from the input data will be used.
#' @param var_labels A character vector to use in output. If NULL (default)
#'   labels from the input data will be used.
#'
#' @param digits Number of decimal places to round to. Default includes 4
#'   decimal places.
#' @return description
#' @export



mm_CompModel_Full <- function(mv_results, row_labels = NULL, var_labels = NULL, digits = 4){


  if(is.null(names(mv_results))){
    stop("Input data must be a list where each model has a unique name.")
  } else {
    mod_labels <- names(mv_results)
  }

  n_mods <- length(mod_labels)

  var_names <- rownames(mv_results[[1]]$summary$table)
  n_vars <- length(var_names)
  out_list <- vector(mode = "list", length = n_vars)

  fill_tab <- data.frame("Rsq" = rep(NA, times = n_mods),
                         "F" = rep(NA, times = n_mods),
                         "Pr(>F)" = rep(NA, times = n_mods)
  )

  out_list <- lapply(out_list, function(x){
    x <- fill_tab
  })


  for(mm in seq_along(mod_labels)){

    for(vv in 1:n_vars){
      tmp_tab <- mv_results[[mod_labels[mm]]]$summary$table[vv,c("Rsq", "F", "Pr(>F)")]
      out_list[[vv]][mm,] <- tmp_tab
    }

  }

  if(is.null(var_labels)){
    var_labels <- var_names
  }

  if(is.null(row_labels)){
    row_labels <- mod_labels
  }

  out_list <- lapply(out_list, function(x){
    rownames(x) <- row_labels
    return(x)
  })

  names(out_list) <- var_labels

  sub_list <- out_list[!names(out_list) %in% c("Residuals", "Total")]

  if(is.numeric(digits)){
    sub_list <- lapply(sub_list, function(x){
      round(x, digits = digits)
    })
  }

  sub_list <- lapply(sub_list, function(x){
    names(x)[names(x)=="Pr(>F)"] <- "p-val"
    return(x)
  })


  return(sub_list)


}


#' Compare Model Metrics
#'
#' Compare key figs (Rsq, p-value, etc) across multiple models.
#'
#' @param mv_results Input mvlm, created by mm_BuildModel (or by using geomorph::procD.lm)
#' @param row_labels A character vector to use in output. If NULL (default) labels from the input data will be used.
#' @param digits Number of decimal places to round to. Default includes 4 decimal places.
#' @return A list containing the results of the mvlm, visualizations of shape
#'   trends along the regression line, and the model itself.
#' @export


mm_CompModel <- function(mv_results, row_labels = NULL, digits=4){


  if(is.null(names(mv_results))){
    stop("Input data must be a list where each model has a unique name.")
  } else {
    mod_labels <- names(mv_results)
  }

  n_mods <- length(mod_labels)

  out_tab <- data.frame("Rsq" = rep(NA, times = n_mods),
                        "F" = rep(NA, times = n_mods),
                        "Pr(>F)" = rep(NA, times = n_mods)
  )


  for(mm in seq_along(mod_labels)){

    tmp_tab <- mv_results[[mod_labels[mm]]]$summary$table[1,c("Rsq", "F", "Pr(>F)")]
    out_tab[mm,] <- tmp_tab


  }


  if(is.null(row_labels)){
    row_labels <- mod_labels
  }

  rownames(out_tab) <- row_labels
  names(out_tab)[names(out_tab)=="Pr(>F)"] <- "p-val"


  if(is.numeric(digits)){
    out_tab <- round(out_tab, digits = digits)
  }

  return(out_tab)


}





