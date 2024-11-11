


library(shiny)
library(moRphomenses)

## load data
example_data <- mm_data

## step 0: align data
aln_11 <- mm_ArrayData(IDs = example_data$ID,
                       DAYS = example_data$CYCLEDAY,
                       VALUE = example_data$E1G,
                       MID = example_data$MIDPOINT,
                       targetLENGTH = 28,
                       targetMID = 16,
                       transformation = "minmax",
                       impute_missing = 3
)


## step 1: pca
shpPCA <- mm_CalcShapespace(aln_11$Shape_data)

# input <- list("x_pc" = 1,
#               "y_pc" = 2)
#

# Define UI for application that draws a histogram
ui <- fixedPage(

  ## 60 px/col ; 78 px/col ; 97 px/col

  # Application title
  titlePanel("mm Explorer"),

    tabsetPanel(
        tabPanel("Shapespace",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                     "Real Shapes",
                     HTML("<br>"),
                     "Individual curves shown in grey.",
                     "The mean shape of the data is located at coordindates (0,0).",
                     "The hypothetical shape of selected coordinates (below) is visualized shown above in green.",
                     "Click anywhere in the scatterplot to select new coordinates.",
                     textOutput("target_coords"),
                     HTML("<br>"),
                     "Hypothetical shapes are visualized at the extrema of each axis.",
                     HTML("<br>"),
                     numericInput("x_pc",
                                  "PC on X axis",
                                  value = 1,
                                  min = 1,
                                  max = ncol(shpPCA$PCA$x)-1,
                                  width = "25%"),

                     numericInput("y_pc",
                                  "PC on Y axis",
                                  value = 2,
                                  min = 2,
                                  max = ncol(shpPCA$PCA$x),
                                  width = "25%")
                     ),

                   mainPanel(width = 9,

                     fixedRow( ## open main pca
                   column(width = 9, ## open LH column
                     fixedRow(
                       column(width = 3, offset = 5,
                       plotOutput("pc_y_max", height = "90px", width = "180px")

                       )),

                     fixedRow(
                       column(width = 3,
                              fixedRow(plotOutput("pad", height = "150px")),
                              fixedRow(plotOutput("pc_x_min", height = "90px", width = "180px"))),
                       column(width = 4, offset = 0,
                              plotOutput("pc_scatter",  height = "240px", width = "240px", click = "pc_click"))


                     ),
                     fixedRow(column(width = 3, offset = 5,
                                     plotOutput("pc_y_min", height = "90px", width = "180px")))
                     ), ## close LH column

                   column( ## open RH column
                     width = 3, offset = 0,
                     fixedRow(plotOutput("target_shp", height = "90px", width = "180px")),

                     fixedRow(plotOutput("pad", height = "50px")),
                     fixedRow(plotOutput("pc_x_max", height = "90px", width = "180px"))
                       ## target shape
                  ) ## close RH column

                  ), ## close main PCA
                 ) ## close main panel

                 ) ## close sidebar layout

                 ), ## close tab panel pca

        tabPanel("subgroups",
                 sidebarLayout(
                   sidebarPanel(width = 3,

                                uiOutput("dyn_k_grps"),


                     radioButtons("cluster_type",
                                  label = "Clustering Algorithm:",
                                  choices = c("Hierarchical", "k-means")
                                  ),
                     checkboxInput("show_diagnostics",
                                   "Show Diagnostic Plots")


                   ), ## close sidebar panel
                   mainPanel(width = 9,

                       column(width = 6,

                              fixedRow(plotOutput("dendro_plot",click = "cut", brush = "inds")),
                              fixedRow(textOutput("debug_text")),
                              fixedRow(plotOutput("sse_plot")),
                              fixedRow(plotOutput("sil_plot"))
                              ),

                       column(width = 3,
                              fixedRow(plotOutput("target_shp2", height = "200px", width = "400px")),
                              fixedRow(plotOutput("tree_shp", height = "200px", width = "400px"))
                       )

                   ) ## close main panel
                 ) ## close layout

                 ),




        tabPanel("Phenotypes",
                plotOutput("phenotypes")
        ), ## close tab panel




      ), ## close tabset

  tags$head(tags$style(
    ".pannel{height:400px; width:600px;}"
  ))

) ## close UI fluidpage

# Define server logic required to draw a histogram
server <- function(input, output) {





  ## PCA SCATTER


  rv_pca_full <- mm_CalcShapespace(aln_11$Shape_data)

  rv_diagnostics <- mm_Diagnostics(rv_pca_full,hide_plots = TRUE)

  all_nodes_xy <- dendextend::get_nodes_xy(rv_diagnostics$TREE)
  all_inds_xy <- all_nodes_xy[all_nodes_xy[,2]==0,]

  rv_pca_select <- reactive({
    rv_pca_full$PCA$x[,c(input$x_pc, input$y_pc)]
  })


  output$pc_scatter <- renderPlot({


    x_r <- range(rv_pca_select()[,1])
    par("mar" = c(.5,.5,.5,.5))

    plot(rv_pca_select()[,c(1,2)], main = "", xlab = "", ylab = "", as= T, bty = "n", xaxt = "n", yaxt = "n", type = "n", xlim = x_r, ylim = x_r, pch = 21, bg = "dark grey")
    abline(h = 0, col = "dark grey", lty = 2, lwd = 2)
    abline(v = 0, col = "dark grey", lty = 2, lwd = 2)

        ##### text(rep(0, 5), yax, labels = round(yax),cex = 1, col = "red",pos = 4)
    ##### text(xax, rep(0, 5), labels = round(xax),cex = 1, col = "blue", pos = 1)

    #### add the points
    points(rv_pca_select(), pch = 16, col = mm_mute_cols("black"))

    ## add target

    abline(h = selected_coords$y, col = "green")
    abline(v = selected_coords$x, col = "green")
    points(selected_coords$x,selected_coords$y, pch = 12, col = "green", cex = 2)



  })


  selected_coords <- reactiveValues(
    "x" = 0,
    "y" = 0
  )
  observeEvent(input$pc_click, {

    selected_coords$x <- input$pc_click$x
    selected_coords$y <- input$pc_click$y
  })


  ## PCA SHAPES

  rv_target <- reactive({

    target_x <- selected_coords$x
    target_y <- selected_coords$y


    out <- mm_coords_to_shape(A = aln_11$Shape_data,
                              PCA = rv_pca_full$PCA,
                              target_PCs = c(input$x_pc, input$y_pc),
                              target_coords = c(target_x, target_y))

    out
  })

  output$target_shp <- renderPlot({
    mm_PlotArray(aln_11$Shape_data,lbl = "Individual shapes",MeanShape = FALSE)
    points(rv_target(), type = "l", col = "green", lwd = 2)

  })

  output$target_shp2 <- renderPlot({
    mm_PlotArray(aln_11$Shape_data,lbl = "Individual shapes",MeanShape = FALSE)
    points(rv_target(), type = "l", col = "green", lwd = 2)

  })


  output$target_coords <- renderText({
    paste(
      "Selected Shape coords",
      paste0("X: ", round(selected_coords$x, 4)),
      paste0("Y: ", round(selected_coords$y, 4)),
      sep = "\n")
  })


  output$pc_x_min <- renderPlot({
    mm_PlotArray(rv_pca_full$Shapes[[input$x_pc]]$min, lbl = paste("Negative X, shape trend", sep = "\n"))
  })

  output$pc_x_max <- renderPlot({
    mm_PlotArray(rv_pca_full$Shapes[[input$x_pc]]$max, lbl = paste("Positive X, shape trend", sep = "\n"))
  })

  output$pc_y_min <- renderPlot({
    mm_PlotArray(rv_pca_full$Shapes[[input$y_pc]]$min, lbl = paste("Negative Y, shape trend", sep = "\n"))
  })

  output$pc_y_max <- renderPlot({
    mm_PlotArray(rv_pca_full$Shapes[[input$y_pc]]$max, lbl = paste("Positive Y, shape trend", sep = "\n"))
  })


  ## CLUSTERING

  output$dyn_k_grps <- renderUI({
    numericInput("k_grps",
                 label = "Subgroups to visualize:",
                 value = brush_ind$max_k, min = 1, max = 25)

  })


  output$sil_plot <- renderPlot({
    mm_SilPlot(rv_pca_full$PCA$x)

  })
  output$sse_plot <- renderPlot({
    mm_ScreePlot(rv_pca_full$PCA$x)

  })
  output$dendro_plot <- renderPlot({

    ## note: hard coded
    rv_diagnostics$TREE <- dendextend::place_labels(rv_diagnostics$TREE, rep("lllll", 75))
    plot(rv_diagnostics$TREE)
    points(brush_ind$in_nodes, col = "orange")
    points(brush_ind$out_nodes, col = "green")
    points(all_inds_xy, col = brush_ind$leaf_cols, pch = 16)

  })

  brush_ind <- reactiveValues(


    "to_keep" = rep(TRUE, 75),
    "in_nodes" = NULL,
    "out_nodes" = all_nodes_xy[all_nodes_xy[,2] > 0,],
    "leaf_cols" = rep("grey", 75),
    "max_k" = 1,
    "grps" = rep(1, 75)
  )



  observeEvent(input$cut, {


    new_groups <- dendextend::cutree(rv_diagnostics$TREE, h = input$cut$y)

    brush_ind$max_k <- max(new_groups)


    all_cols <- rainbow(brush_ind$max_k, 0.8, 0.8)
    new_cols <- character(length(new_groups))

    for(ii in seq_along(all_cols)){
      new_cols[new_groups==ii] <- all_cols[ii]
    }

    brush_ind$leaf_cols <- new_cols[order.dendrogram(rv_diagnostics$TREE)]
    brush_ind$grps <- new_groups

  })


  observeEvent(input$k_grps, {


    new_groups <- dendextend::cutree(rv_diagnostics$TREE, k = input$k_grps)

    brush_ind$max_k <- max(new_groups)


    all_cols <- rainbow(brush_ind$max_k, 0.8, 0.8)
    new_cols <- character(length(new_groups))

    for(ii in seq_along(all_cols)){
      new_cols[new_groups==ii] <- all_cols[ii]
    }

    brush_ind$leaf_cols <- new_cols[order.dendrogram(rv_diagnostics$TREE)]
    brush_ind$grps <- new_groups

  })

  observeEvent(input$inds, {



     select_inds_xy <- (all_inds_xy[,1] < input$inds$xmax) & (all_inds_xy[,1] > input$inds$xmin)
    brush_ind$to_keep <- select_inds_xy[order.dendrogram(rv_diagnostics$TREE)]

    brush_ind$in_nodes <- all_nodes_xy[((all_nodes_xy[,1] < input$inds$xmax) & (all_nodes_xy[,1] > input$inds$xmin)) & all_nodes_xy[,2] > 0,]

    brush_ind$out_nodes <- all_nodes_xy[!((all_nodes_xy[,1] < input$inds$xmax) & (all_nodes_xy[,1] > input$inds$xmin)) & all_nodes_xy[,2] > 0,]

  })

  output$tree_shp <- renderPlot({
    sub_aln <- aln_11$Shape_data[,,brush_ind$to_keep]
    mm_PlotArray(sub_aln,lbl = "Individual shapes",MeanShape = FALSE)
    tmp_mshp <- apply(sub_aln, c(1,2), mean)
    points(tmp_mshp, type = "l", col = "orange", lwd = 2)
  })



  output$phenotypes <- renderPlot({

    layout(matrix(1:brush_ind$max_k, ncol = 1))

    all_cols <- rainbow(brush_ind$max_k, 0.4, 0.8, alpha = .4)
    m_cols <- rainbow(brush_ind$max_k, 0.8, 0.4)
    for(nn in 1:brush_ind$max_k){
      mm_PlotArray(aln_11$Shape_data[,,brush_ind$grps==nn], MeanCol = m_cols[nn], AllCols = all_cols[nn],lbl = paste("Group", nn, "of", brush_ind$max_k))
    }


  })

  output$debug_text <- renderText({

  })

}

# Run the application
shinyApp(ui = ui, server = server)


