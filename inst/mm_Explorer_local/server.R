



library(shiny)
library(dplyr)
library(moRphomenses)


# Define server logic required to draw a histogram
function(input, output, session) {
  # Setup and global params ---------------------

  ##  940 px wid total; but this should be handled by the ui
  # pc_shape_height <- 120
  # pc_shape_width <- 240


  ## try to maintain one single rv list
  rv <- reactiveValues(
    "id_var" = "ID",
    "x_var" = "CYCLEDAY",
    "y_var" = "E1G",
    "aln_var" = "MIDPOINT",
    "mm_dataset" = mm_data,
    "n_inds" = 75,
    "data_ready" = TRUE
  )


  output$dyn_var_assign <- renderUI({
    var_choices <- names(rv$mm_dataset)

    var_select_ui <- list(
      selectInput(
        "assign_id",
        "Choose ID variable",
        choices = var_choices,
        selected = var_choices[1]
      ),
      selectInput(
        "assign_x",
        "Choose Time variable",
        choices = var_choices,
        selected = var_choices[2]
      ),
      selectInput(
        "assign_y",
        "Choose Value variable",
        choices = var_choices,
        selected = var_choices[4]
      ),
      selectInput(
        "assign_aln",
        "Choose Alignment variable",
        choices = var_choices,
        selected = var_choices[3]
      )
    )

    var_select_ui
  })


  observeEvent(input$file_select, {
    filepath <- input$file_select
    rv$mm_dataset <- read.csv(filepath$datapath)
    ## Reset to false
    rv$data_ready <- FALSE
  })

  observeEvent(input$assign_id, {
    rv$id_var <- input$assign_id
    rv$n_inds <-  nrow(count(rv$mm_dataset, !!sym(eval(
      input$assign_id
    ))))
  })

  observeEvent(input$assign_x, {
    rv$x_var <- input$assign_x
  })

  observeEvent(input$assign_y, {
    rv$y_var <- input$assign_y
  })

  observeEvent(input$assign_aln, {
    rv$aln_var <- input$assign_aln
  })


  observeEvent(input$go, {

    rv$data_ready <- TRUE

  })

  # alignment ----------------------------------

  aln_11 <- reactive({
    input$go

    tmp_id <-  pull(rv$mm_dataset, !!sym(eval(rv$id_var)))
    tmp_x <-  pull(rv$mm_dataset, !!sym(eval(rv$x_var)))
    tmp_y <- pull(rv$mm_dataset, !!sym(eval(rv$y_var)))
    tmp_aln <- pull(rv$mm_dataset, !!sym(eval(rv$aln_var)))

    mm_ArrayData(
      IDs = tmp_id,
      DAYS = tmp_x,
      VALUE = tmp_y,
      MID = tmp_aln,
      targetLENGTH = input$tar_length,
      targetMID = input$tar_alignment,
      transformation = "minmax",
      impute_missing = input$knn_n
    )

  })




  ## renderPlot aln_scl ----------------------------

  output$dyn_pick1 <- renderUI({
    ind_list <- unique(pull(rv$mm_dataset, !!sym(eval(rv$id_var))))

    selectInput("ind_pick1", label = "Individual to highlight", choices = ind_list)

  })


  output$aln_scl <- renderPlot({
    if(!rv$data_ready){
      return()
    }
    mm_PlotArray(aln_11()$Shape_data,
                 MeanShape = FALSE,
                 lbl = "")

    points(
      aln_11()$Shape_data[, , input$ind_pick1],
      type = "l",
      col = "green",
      lwd = 2
    )


  })

  output$ind_missing <- renderPlot({
    if(!rv$data_ready){
      return()
    }
    mm_PlotArray(aln_11()$shape_data_wNA,
                 MeanShape = FALSE,
                 lbl = "")
    points(aln_11()$shape_data_wNA[, , input$ind_pick1], pch = 16, col = "dark grey")
    points(aln_11()$Shape_data[, , input$ind_pick1], col = "red")


  })

  ## renderPrint output$aln_summary -------------------

  output$aln_summary <- renderPrint({
    if(!rv$data_ready){
      return()
    }
    tmp <- print_summary(aln_11())
    cat(tmp)

  })

  output$n_miss_summary <- renderTable({
    if(!rv$data_ready){
      return()
    }
    tmp_tab <- data.frame(table(aln_11()$knn_info$nmis))
    names(tmp_tab) <- c("Missing Days", "Frequency")
    tmp_tab

  })



  # shapespace ----------------------------------

  ## CalcShapespcae ----------------------------

  rv_pca_full <- reactive({
    if(!rv$data_ready){
      return()
    }

    mm_CalcShapespace(aln_11()$Shape_data)
  })



  ## Reactive PCs to plot ----------------------

  rv_pca_select <- reactive({
    if(!rv$data_ready){
      return()
    }

    rv_pca_full()$PCA$x[, c(input$x_pc, input$y_pc)]

  })


  ## renderPlot output$pc_scatter ----------------------
  output$pc_scatter <- renderPlot(# height = pc_shape_width, width = pc_shape_width,
    {
      if(!rv$data_ready){
        return()
      }

      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar),add = TRUE)

      x_r <- range(rv_pca_select()[, 1])
      par("mar" = c(.5, .5, .5, .5))

      plot(
        rv_pca_select()[, c(1, 2)],
        main = "",
        xlab = "",
        ylab = "",
        as = T,
        bty = "n",
        xaxt = "n",
        yaxt = "n",
        type = "n",
        xlim = x_r,
        ylim = x_r,
        pch = 21,
        bg = "dark grey"
      )

      abline(
        h = 0,
        col = "dark grey",
        lty = 2,
        lwd = 2
      )
      abline(
        v = 0,
        col = "dark grey",
        lty = 2,
        lwd = 2
      )

      ##### text(rep(0, 5), yax, labels = round(yax),cex = 1, col = "red",pos = 4)
      ##### text(xax, rep(0, 5), labels = round(xax),cex = 1, col = "blue", pos = 1)

#
#       if(input$k_grps > 1){
#         set.seed(13)
#         tmp_grps <- kmeans(rv_pca_select(),centers = brush_ind$max_k)
#         cols <- rainbow(brush_ind$max_k)
#         for(gg in 1:brush_ind$max_k){
#           points(rv_pca_select()[tmp_grps$cluster == gg,],
#                  pch = 16,
#                  col = mm_mute_cols(cols[gg])
#                  )
#         }
#       } else {
#         #### add the points
#         points(rv_pca_select(), pch = 16, col = mm_mute_cols("black"))
#       }

      points(rv_pca_select(), pch = 16, col = mm_mute_cols("black"))

      ## add target
      abline(h = selected_coords$y, col = "green")
      abline(v = selected_coords$x, col = "green")
      points(
        selected_coords$x,
        selected_coords$y,
        pch = 12,
        col = "green",
        cex = 2
      )

    })

  ## reactiveValues selected_coords ---------------------------
  #### INITIALIZE an object that will change interacticely

  selected_coords <- reactiveValues(
    "x" = 0,
    "y" = 0)

  ## observeEvent input$pc_click -------------------------
  observeEvent(input$pc_click, {
    selected_coords$x <- input$pc_click$x
    selected_coords$y <- input$pc_click$y
  })


  ## Interactive PC shape on click

  ### reactive pc_target_shp ------------------

  pc_target_shp <- reactive({
    if(!rv$data_ready){
      return()
    }
    target_x <- selected_coords$x
    target_y <- selected_coords$y

    out <- mm_coords_to_shape(
      A = aln_11()$Shape_data,
      PCA = rv_pca_full()$PCA,
      target_PCs = c(input$x_pc, input$y_pc),
      target_coords = c(target_x, target_y)
    )

    out
  })


  ## renderPlot output$target_shp --------------------------

  output$target_shp <- renderPlot(# height = pc_shape_height, width = pc_shape_width,
    {
      if(!rv$data_ready){
        return()
      }
      mm_PlotArray(aln_11()$Shape_data,
                   lbl = "Individual shapes",
                   MeanShape = FALSE)
      points(
        pc_target_shp(),
        type = "l",
        col = "green",
        lwd = 2
      )

    })

  ## renderText output$target_coords ---------------------

  output$target_coords <- renderText({
    if(!rv$data_ready){
      return()
    }
    paste("Selected Shape coords",
          paste0("X: ", round(selected_coords$x, 4)),
          paste0("Y: ", round(selected_coords$y, 4)),
          sep = "\n")
  })


  ## render PC shapes ------------------------------

  ### renderPlot output$pc_x_min -------------------------
  output$pc_x_min <- renderPlot({
      if(!rv$data_ready){
        return()
      }
      mm_PlotArray(rv_pca_full()$Shapes[[input$x_pc]]$min, lbl = "")
    })


  ### renderPlot output$pc_x_max -------------------------
  output$pc_x_max <- renderPlot({
      if(!rv$data_ready){
        return()
      }
      mm_PlotArray(
        rv_pca_full()$Shapes[[input$x_pc]]$max, lbl = "")
    })


  ### renderPlot output$pc_y_min -------------------------
  output$pc_y_min <- renderPlot(# height = pc_shape_height, width = pc_shape_width,
    {
      if(!rv$data_ready){
        return()
      }
      mm_PlotArray(rv_pca_full()$Shapes[[input$y_pc]]$min, lbl = "")
    })


  ### renderPlot output$pc_y_max -------------------------
  output$pc_y_max <- renderPlot(# height = pc_shape_height, width = pc_shape_width,
    {
      if(!rv$data_ready){
        return()
      }
      mm_PlotArray(rv_pca_full()$Shapes[[input$y_pc]]$max, lbl = "")
    })


  # CLUSTERING --------------------------------------------


  ## Reactive brush_ind  ----------------------------

  #### Reactive object that tracks which indibviduals
  #### are selected from dendro

  ## need to update this so it can initialize without all_nodes_xy or without hard coding n inds

  brush_ind <- reactiveValues(
    "to_keep" = rep(TRUE, 75),
    "in_nodes" = NULL,
    "out_nodes" = NULL,
    "leaf_cols" = rep("grey", 75),
    "max_k" = 1,
    "grps" = rep(1, 75)
  )


  observeEvent(input$inds, {
    if(!rv$data_ready){
      return()
    }
    brush_ind$out_nodes <- all_nodes_xy()[all_nodes_xy()[, 2] > 0, ]
  })
  ## Reactive mm_Diagnostics ------------------------

  #### need to update
  rv_diagnostics <- reactive({
    if(!rv$data_ready){
      return()
    }
    mm_Diagnostics(rv_pca_full(), hide_plots = TRUE)
  })

  all_nodes_xy <- reactive({
    if(!rv$data_ready){
      return()
    }
    dendextend::get_nodes_xy(rv_diagnostics()$TREE)

  })


  all_inds_xy <- reactive({
    if(!rv$data_ready){
      return()
    }
    all_nodes_xy()[all_nodes_xy()[, 2] == 0, ]

  })



  ## RenderUI dyn_k_grps --------------------------------
  output$dyn_k_grps <- renderUI({
    numericInput(
      "k_grps",
      label = "Subgroups to visualize:",
      value = brush_ind$max_k,
      min = 1,
      max = 25
    )

  })


  ## Diagnostics --------------------------------------

  ### Render sil_plot ---------------------------------
  output$sil_plot <- renderPlot({
    if(!rv$data_ready){
      return()
    }

    if(input$show_diagnostics){
      mm_SilPlot(rv_pca_full()$PCA$x)
    }
  })

  ### Render sse_plot ---------------------------------
  output$sse_plot <- renderPlot({
    if(!rv$data_ready){
      return()
    }

    if(input$show_diagnostics){
      mm_ScreePlot(rv_pca_full()$PCA$x)
      }
    })


  output$all_shp1 <- renderPlot({
    if(!rv$data_ready){
      return()
    }
    mm_PlotArray(aln_11()$Shape_data,
                 MeanShape = TRUE,
                 lbl = "Whole Sample")

  })


  ### Render dendro_plot -------------------------------
  output$dendro_plot <- renderPlot({
    if(!rv$data_ready){
      return()
    }
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar),add = TRUE)

    tmp_TREE <-  dendextend::place_labels(rv_diagnostics()$TREE, rep("", 75))
    tmp_inds <- all_inds_xy()

    par("mar" = c(0,2,0,0))
    plot(tmp_TREE)
    points(brush_ind$in_nodes, col = "black", cex = 1.2, pch = 16)
    points(brush_ind$in_nodes, col = "orange", cex = 1.1, pch = 16)
    points(brush_ind$out_nodes, col = "yellow", cex = 1.1, pch = 16)
    points(tmp_inds, col = brush_ind$leaf_cols, pch = 16)

  })



  ## ObserveEvent input$cut -------------------------

  observeEvent(input$cut, {
    new_groups <- dendextend::cutree(rv_diagnostics()$TREE, h = input$cut$y)

    brush_ind$max_k <- max(new_groups)

    all_cols <- rainbow(brush_ind$max_k, 0.8, 0.8)
    new_cols <- character(length(new_groups))

    for (ii in seq_along(all_cols)) {
      new_cols[new_groups == ii] <- all_cols[ii]
    }

    brush_ind$leaf_cols <- new_cols[order.dendrogram(rv_diagnostics()$TREE)]
    brush_ind$grps <- new_groups

  })


  ## ObserveEvent input$k_grps ------------------------------

  observeEvent(input$k_grps, {
    new_groups <- dendextend::cutree(rv_diagnostics()$TREE, k = input$k_grps)

    brush_ind$max_k <- max(new_groups)

    all_cols <- rainbow(brush_ind$max_k, 0.8, 0.8)
    new_cols <- character(length(new_groups))

    for (ii in seq_along(all_cols)) {
      new_cols[new_groups == ii] <- all_cols[ii]
    }

    brush_ind$leaf_cols <- new_cols[order.dendrogram(rv_diagnostics()$TREE)]
    brush_ind$grps <- new_groups

  })

  ## ObserveEvent input$inds --------------------------------
  observeEvent(input$inds, {
    select_inds_xy <- (all_inds_xy()[, 1] < input$inds$xmax) &
      (all_inds_xy()[, 1] > input$inds$xmin)
    brush_ind$to_keep <- select_inds_xy[order.dendrogram(rv_diagnostics()$TREE)]

    brush_ind$in_nodes <- all_nodes_xy()[((all_nodes_xy()[, 1] < input$inds$xmax) &
                                            (all_nodes_xy()[, 1] > input$inds$xmin)) & all_nodes_xy()[, 2] == 0, ]

    brush_ind$out_nodes <- all_nodes_xy()[!((all_nodes_xy()[, 1] < input$inds$xmax) &
                                              (all_nodes_xy()[, 1] > input$inds$xmin)) & all_nodes_xy()[, 2] == 0, ]

  })


  ## Render tree shape from brush ---------------------------------

  output$tree_shp <- renderPlot({
    if(!rv$data_ready){
      return()
    }
    sub_aln <- aln_11()$Shape_data[, , brush_ind$to_keep]
    mm_PlotArray(sub_aln, lbl = "Selected Group", MeanShape = FALSE)
    tmp_mshp <- apply(sub_aln, c(1, 2), mean)
    points(tmp_mshp,
           type = "l",
           col = "orange",
           lwd = 2)

  })



  # Render Phenotypes ------------------------------------------

  output$phenotypes <- renderPlot(height = reactive(200 * brush_ind$max_k),
                                  width = 600,
                                  {
                                    if(!rv$data_ready){
                                      return(NULL)
                                    }

                                    on.exit(layout(matrix(1)))

                                    layout(matrix(1:brush_ind$max_k, ncol = 1))

                                    all_cols <- rainbow(brush_ind$max_k, 0.4, 0.8, alpha = .4)
                                    m_cols <- rainbow(brush_ind$max_k, 0.8, 0.4)
                                    for (nn in 1:brush_ind$max_k) {
                                      mm_PlotArray(
                                        aln_11()$Shape_data[, , brush_ind$grps == nn],
                                        MeanCol = m_cols[nn],
                                        AllCols = all_cols[nn],
                                        lbl = paste("Group", nn, "of", brush_ind$max_k)
                                      )
                                    }

                                  })


  output$pheno_summary <- renderPrint({
    if(!rv$data_ready){
      return()
    }

      pheno_summary <- print_summary(aln_11(), brush_ind$grps)
    lapply(pheno_summary, cat)

  })


  # Debug ------------------------------------
  output$debug_text <- renderText({

  })

}
