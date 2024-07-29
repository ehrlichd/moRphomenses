#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

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

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("mm Explorer"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          ## PCA controls
          numericInput("x_pc",
                       "PC on X axis",
                       value = 1,
                       min = 1,
                       max = ncol(shpPCA$PCA$x)-1),
          numericInput("y_pc",
                       "PC on Y axis",
                       value = 2,
                       min = 2,
                       max = ncol(shpPCA$PCA$x)),
          ## k-means controls
          ## sample info
          ## hypo testing

        ),

        mainPanel(

          tabsetPanel(
            tabPanel("PCA",
                     ##
                     fixedPage(



                       ## top row
                       fixedRow(
                         column(width = 4,
                         ),
                         column(width = 4,
                                fluidRow(class = "shp_row",
                                         plotOutput("pad")),
                                fluidRow(class="shp_row",
                                         plotOutput("pc_y_max")
                                         )

                         ),
                         column(width = 4,
                                plotOutput("target_shape")
                         )
                       ),

                       ## mid row
                       fixedRow(
                         column(width = 4,
                                fixedRow(),
                                fixedRow(class = "shp_row",
                                         plotOutput("pc_x_min")
                                         ),
                                fixedRow()

                         ),
                         column(width = 4,
                                plotOutput("pc_scatter")
                         ),
                         column(width = 4,
                                fixedRow(),
                                fixedRow(class = "shp_row",
                                  plotOutput("pc_x_max")
                                  ),
                                fixedRow()

                         ),
                       ),

                       ## bottom row
                       fixedRow(class = "shp_row",
                         column(width = 4
                         ),
                         column(width = 4,
                                plotOutput("pc_y_min")
                         ),
                         column(width = 4
                         )

                       ),
                       fixedRow(
                         ## interactive text output
                       )
                     ) # close fluid page

                     ), ## close tab panel

            tabPanel("Groups"
                     ), ## close tab panel
            tabPanel("Summary Info"
                     ), ## close tab panel
            tabPanel("Hypothesis"
                     ) ## close tab panel
          ) ## close tabset
          ) ## close main panel

    ), ## close sidebar layout
    tags$head(tags$style(
      ".shp_row{height:50%;}"
    ))
) ## close UI fluidpage

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$pad <- renderPlot({

  })

    output$pc_scatter <- renderPlot({

      x_r <- range(shpPCA$PCA$x[,input$x_pc])
      par("mar" = c(.5,.5,.5,.5))

      plot(shpPCA$PCA$x[,c(input$x_pc, input$y_pc)], main = "", xlab = "", ylab = "", as= T, bty = "n", xaxt = "n", yaxt = "n", type = "n", xlim = x_r, ylim = x_r, pch = 21, bg = "dark grey")
      abline(h = 0, col = "grey", lty = 2)
      abline(v = 0, col = "grey", lty = 2)




      ##### text(rep(0, 5), yax, labels = round(yax),cex = 1, col = "red",pos = 4)
      ##### text(xax, rep(0, 5), labels = round(xax),cex = 1, col = "blue", pos = 1)

      #### add the points
      points(shpPCA$PCA$x[,c(input$x_pc,input$y_pc)], pch = 16, col = mm_mute_cols("black"))


    })

    output$pc_x_min <- renderPlot({
      mm_PlotArray(shpPCA$Shapes[[input$x_pc]]$min, lbl = "PC X Min")
    })

    output$pc_x_max <- renderPlot({
      mm_PlotArray(shpPCA$Shapes[[input$x_pc]]$max, lbl = "PC X Max")
    })

    output$pc_y_min <- renderPlot({
      mm_PlotArray(shpPCA$Shapes[[input$y_pc]]$max, lbl = "PC Y Min")
    })

    output$pc_y_max <- renderPlot({
      mm_PlotArray(shpPCA$Shapes[[input$y_pc]]$max, lbl = "PC Y Max")
    })
    output$target_shape <- renderPlot({


    })



    output$sil_wid <- renderPlot({

    })
    output$group_sse <- renderPlot({

    })
    output$dendro <- renderPlot({

    })
}

# Run the application
shinyApp(ui = ui, server = server)
