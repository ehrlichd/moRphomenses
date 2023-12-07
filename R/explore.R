
# Explore #####

#' Define Shapespace of aligned dataset.
#'
#' Conduct PCA of shape data and visualize major shape trends.
#'
#' @param dat A 3D array of shape data to be analyzed.
#' @param max_Shapes The maximum amount of PCs to visualize. Default 10.
#'
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

    out$Shapes[[i]] <- geomorph::shape.predictor(
      A = dat,
      x = out$PCA$x[,i],
      min = min(out$PCA$x[,i]),
      max = max(out$PCA$x[,i]))

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

#' Run a suite of diagnostic analysis/plots
#'
#'
#'
#' @param dat A 3D array or a mmPCA object (output of mm_CalcShapespace).
#' @param max_PC_viz Maximum number of PCs to include in visualizations (EG Eigenplots, or shape trends.
#' @param max_PC_calc By default (NULL), all PCs will be included in calculations. However, if fewer PCs are required users may specify an integer, n, to get the first n PCS.
#' @param hide_plots By default (FALSE), helpful visuals are plotted.
#' @export
#'
mm_Diagnostics <- function(dat, max_PC_viz=10, max_PC_calc=NULL, hide_plots = FALSE){


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



  ## DEMOGRAPHICS??

  class(out) <- "mmDiag"

  return(out)

}





#' Generate Phenotypes
#'
#' Partition sample into clusters, based on information from
#' @param dat Either an Array of shape data, an mmPCA object, or an mmDiag object.
#' @param kgrps A non-negative integer of sub-groups to draw. kgrps=1 will provide results for the whole input dat.
#' @param cuttree_h Optional. Draw clusters by splitting the tree at a given height, h.
#' @param cuttree_k Optional. Draw clsuters by splitting the tree into number of branches, k
#' @param plot_figs Optional. Default = TRUE, plot phenotypes for each set(s) of subgroups.
#'
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
    out$cth_grps <- dendextend::cutree(out$TREE, h = cuttree_h)
  }

  if(!is.null(cuttree_k)){
    out$ctk_grps <- dendextend::cutree(out$TREE, k = cuttree_k)
  }


  ## plotting

  if(plot_figs){

    ## plot mshp
    mm_PlotArray(out$ALN,lbl = "Whole Sample")


    if(!is.null(kgrps)){

      uu_kgrps <- unique(out$k_grps$cluster)
      ll_kgrps <- length(uu_kgrps)

      all_cols <- rainbow(ll_kgrps, .4, .8, alpha = .4)
      m_cols <- rainbow(ll_kgrps, .8, .4)

      for(k in uu_kgrps){
        mm_PlotArray(out$ALN[,,out$k_grps$cluster==k], "lbl" = paste0("kmeans_", k), MeanCol = m_cols[k], AllCols = all_cols[k])
      }
    }

    if(!is.null(cuttree_h)){
      uu_cth_grps <- unique(out$cth_grps)
      ll_cth_grps <- length(uu_cth_grps)

      all_cols <- rainbow(ll_cth_grps, .4, .8, alpha = .4)
      m_cols <- rainbow(ll_cth_grps, .8, .4)

      for(h in uu_cth_grps){
        mm_PlotArray(out$ALN[,,out$cth_grps==h], "lbl" = paste0("hcut_", h), MeanCol = m_cols[h], AllCols = all_cols[h])
      }
    }

    if(!is.null(cuttree_k)){
      uu_ctk_grps <- unique(out$ctk_grps)
      ll_ctk_grps <- length(uu_ctk_grps)

      all_cols <- rainbow(ll_ctk_grps, .4, .8, alpha = .4)
      m_cols <- rainbow(ll_ctk_grps, .8, .4)

      for(k in uu_ctk_grps){
        mm_PlotArray(out$ALN[,,out$ctk_grps==k], "lbl" = paste0("kcut_", k), MeanCol = m_cols[k], AllCols = all_cols[k])
      }
    }


  }


}

#' Check Imputation
#'
#' Plot Raw (alinged) data along side by side with imputed data.
#'
#'
#' @name mm_CheckImputation
#' @param A1 An aligned array, containing missing data (presumably made with  \code{\link{mm_ArrayData}}`$shape_data_wNA`).
#' @param A2 An aligned and imputed array (presumably made with \code{\link{mm_ArrayData}}`$Shape_data`).
#' @param ObO One-by-One. If TRUE (default), individuals will be plotted one at a time, requiring the user to advance/exit the operation. If FALSE, all plots will be generated at once to be browsed/exported from  the Plot History
#'
#'
#'
#' @export
#'
## NEEDS TO BE MODIFIED

mm_CheckImputation <- function(A1, A2, ObO = TRUE){
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








# Visualize #####

#' Plot Array
#' Plot individuals and optionally mean form
#' @name mm_PlotArray
#'
#' @param A An array to be plotted
#' @param MeanShape Logical. Should the Mean Shape be calculated and plotted
#' @param AllCols Either a single color for all individuals, or a vector specifying colors for each individual. If NULL (default) individuals will be plotted in grey
#' @param MeanCol A single color for the mean shape. If Null (default) mean shape will be plotted in black
#' @param plot_type Should the data be plotted as points or lines.
#' @param lbl A title (main =) for the plot. If NULL (default) the name of the array will be used.
#' @param yr Y-range, in the form c(0,100)
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
                         yr = NULL){
  ## added lbl and yr arguments
  ## added flexibility for missing data (i think)
  ## can handle groups of n=1

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



  if(n == 1){
    mshp <- A

    if(plot_type == "points"){
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

par("mar" = c(1,1,2,.5))
  plot(
    mar = c(1,2,1,1),
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
    mm_PlotArray(A=A[,,grps==q],AllCols = cols[q], MeanCol = mcols[q], plot_type = "lines", lbl = paste("g",q,sep=""), yr = range(A[,2,],na.rm = T))
  }
  layout(matrix(1))
}


#' Visualize PC axes
#'
#' Plot a scatterplot and vizualize shape change across the X axis.
#'
#' Meant to be a quick diagnostic plot with minimal customization.
#'
#' @param mmPCA Output of \code{mm_CalcShapespace}, containing a PCA object with PC shapes
#' @param xPC The PC to be plotted on the x axis. If yPC is left null, a univariate density distribution will be plotted with min/max shapes.
#' @param yPC The PC to be plotted on the y axis.
#' @param yr The y-xis range, in the format c(0,1)
#' @param cols A vector of colors of length n, for use in scatterplot.
#' @param title To be used for the plot
#' @param png_dir A file path to a directory in which to save out PNG figures. Names will be automatically assigned based on input PC(s).
#'
#' @export



mm_Viz_shapespace <- function(mmPCA, xPC = 1, yPC = 2, yr = c(0,1), cols = NULL, title = "", png_dir = NULL){


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

    layout(matrix(1:3, ncol =3))

    plot(mmPCA$Shapes[[xPC]]$min, col = "cyan", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "min"))
    points(mmPCA$Shapes$GrandM, col = "grey", type = "l")


    plot(density(mmPCA$PCA$x[,xPC]), main = title, xlab = x_lab, ylab = y_lab, as= T, bty = "n", xaxt = "n", yaxt = "n", type = "n")
    abline(h = 0, col = "grey", lty = 2)
    abline(v = 0, col = "grey", lty = 2)

    #text(rep(0, 5), yax, labels = round(yax),cex = 1, col = "red",pos = 4)
    #text(xax, rep(0, 5), labels = round(xax),cex = 1, col = "blue", pos = 1)


    plot(mmPCA$Shapes[[xPC]]$max, col = "blue", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "max"))
    points(mmPCA$Shapes$GrandM, col = "grey", type = "l")


    ## Bivariate plots - big

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


    layout(matrix(c(0,0,0, 5,5,5,0, 0,0,0,
                    0,0,0, 5,5,5,0, 6,6,6,
                    0,0,0, 5,5,5,0, 6,6,6,
                    0,0,0, 5,5,5,0, 0,0,0,

                    1,1,1, 3,3,3,3, 4,4,4,
                    1,1,1, 3,3,3,3, 4,4,4,
                    1,1,1, 3,3,3,3, 4,4,4,
                    1,1,1, 3,3,3,3, 4,4,4,

                    0,0,0, 2,2,2,0, 7,7,0,
                    0,0,0, 2,2,2,0, 7,7,0,
                    0,0,0, 2,2,2,0, 7,7,0,
                    0,0,0, 2,2,2,0, 0,0,0
                    ), byrow = T, ncol=10, nrow = 12))

    par("mar" = c(2,2,2,1))

    ## xpc min
    plot(mmPCA$Shapes[[xPC]]$min, col = "cyan", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", xPC, "min"))
    points(mmPCA$Shapes$GrandM, col = "grey", type = "l")

    ## MODIFY FOR ypc min
    plot(mmPCA$Shapes[[yPC]]$min, col = "cyan", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", yPC, "min"))
    points(mmPCA$Shapes$GrandM, col = "grey", type = "l")




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
    points(mmPCA$Shapes$GrandM, col = "grey", type = "l")

    ## MODIFY for yPC max
    plot(mmPCA$Shapes[[yPC]]$max, col = "blue", bty = "n", xlab = "", ylab = "", type = "l", lwd = 2, ylim = yr, main = paste("PC", yPC, "max"))
    points(mmPCA$Shapes$GrandM, col = "grey", type = "l")
  }


  ## PLOT Summary stuff

  par("mar" = c(0,0,.1,.1))
  barplot(mmPCA$PCA$eigs[1,1:10])
  barplot(mmPCA$PCA$eigs[1,c(xPC, yPC)], col = "red", add= T)


  plot(1:12, xlab = "", ylab = "", as= T, bty = "n", xaxt = "n", yaxt = "n", type = "n")
  text(x=4, y=12, labels = paste0("n:   " , nrow(mmPCA$PCA$x)))
  text(x=5, y=11, labels = c("Eigs.    Cumul"))
  pad <- 1
  for(ii in 1:10){
    pad <- pad + 1
    text(x=7, y=11-pad,
         labels = paste0(
           paste0("PC", ii, ":  "),
           paste(t(round(mmPCA$PCA$eigs[,ii],2)),collapse = "    "))
         )
  }

  if(!is.null(png_dir)){
    dev.off()
  }

  layout(matrix(1))
  par("mar" = c(5.1,4.1,4.1,2.1))

}






# Test #####

#' Build, implement, visualize multivariate linear model.
#'
#'


mm_Model <- function(){


}





