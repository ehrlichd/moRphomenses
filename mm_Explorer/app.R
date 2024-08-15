


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

                   mainPanel(

                     fixedRow( ## open main pca
                   column(width = 8, ## open LH column
                     fixedRow(column(width = 2, offset = 4,
                       plotOutput("pc_y_max", height = "100px", width = "200px"))),
                     fixedRow(
                       column(width = 2,
                              fixedRow(plotOutput("pad", height = "150px")),
                              fixedRow(plotOutput("pc_x_min", height = "100px", width = "200px"))),
                       column(width = 4, offset = 1,
                              plotOutput("pc_scatter",  height = "400px", width = "400px", click = "pc_click")),
                       column(width = 2, offset = 2,
                              fixedRow(plotOutput("pad", height = "150px")),
                              fixedRow(plotOutput("pc_x_max", height = "100px", width = "200px")))
                     ),
                     fixedRow(column(width = 2, offset = 4,
                                     plotOutput("pc_y_min", height = "100px", width = "200px")))
                     ), ## close LH column

                   column( ## open RH column
                     width = 4, offset = 0,
                     plotOutput("target_shp", height = "200px", width = "400px")
                       ## target shape
                  ) ## close RH column

                  ), ## close main PCA
                 ) ## close main panel

                 ) ## close sidebar layout

                 ), ## close tab panel pca

        tabPanel("subgroups",
                 sidebarLayout(
                   sidebarPanel(

                     numericInput("k_grps",
                                  label = "Subgroups to visualize:",
                                  value = 1, min = 1, max = 25),
                     radioButtons("cluster_type",
                                  label = "Clustering Algorithm:",
                                  choices = c("Hierarchical", "k-means")
                                  )


                   ), ## close sidebar panel
                   mainPanel(

                       column(width = 5,
                              plotOutput("dendro_plot")),
                       column(width = 3,
                              plotOutput("sse_plot"),
                              plotOutput("sil_plot"))
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
                              PCA = rv_pca_full,
                              target_PCs = c(input$x_pc, input$y_pc),
                              target_coords = c(target_x, target_y))

    out
  })

  output$target_shp <- renderPlot({
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

  output$sil_plot <- renderPlot({
    mm_SilPlot(rv_pca_full$PCA$x)

  })
  output$sse_plot <- renderPlot({
    mm_ScreePlot(rv_pca_full$PCA$x)

  })
  output$dendro_plot <- renderPlot({
    tmp <- mm_Diagnostics(rv_pca_full,hide_plots = TRUE)
    ## note: hard coded
    tmp$TREE <- dendextend::place_labels(tmp$TREE, rep("lllll", 75))
    plot(tmp$TREE)

  })



  output$phenotypes <- renderPlot({

  })

  output$debug_text <- renderText({

    # names(rv_pca_select())

  })

}

# Run the application
shinyApp(ui = ui, server = server)


